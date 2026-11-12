"""Explicitly synthetic Collins DAT/archive fixtures; no measured arrays are vendored."""
import copy
import hashlib
import io
from pathlib import Path
import tarfile

import pytest

from collinsflow import experimental as ex
from collinsflow.data import validate_dataset
from collinsflow.provenance import atomic_json, file_hash, read_json


ROOT = Path(__file__).resolve().parents[1]
MEMBER = "HEPdata_tables/Collins.dat"


def synthetic_text():
    lines = ["Collins asymmetry (COMPASS deuteron 2022 data (~70%))", ""]
    for charge in ("+", "-"):
        for projection in ("x", "z", "pT"):
            lines += [f"h^{charge}; {projection}-dependence", "",
                      f"{projection}-bins <x> <y> <z> <p_T> <W> <Q^2> Asym er_stat er_tot"]
            for index in range(2):
                means = {"x": "0.150", "y": "0.330", "z": "0.420", "pT": "0.450"}
                low, high, mean = {
                    "x": [("0.100", "0.200", "0.150"), ("0.200", "0.300", "0.250")],
                    "z": [("0.400", "0.500", "0.420"), ("0.500", "0.600", "0.520")],
                    "pT": [("0.400", "0.600", "0.450"), ("0.600", "0.800", "0.650")],
                }[projection][index]
                means[projection] = mean
                # Equal rounded stat/total fields are intentionally retained without decomposition.
                lines.append(f"[{low};{high}] {means['x']} {means['y']} {means['z']} {means['pT']} "
                             "6.200 2.400 -0.0170 0.0040 0.0040")
            lines.append("")
    return "\n".join(lines) + "\n"


def fixture(tmp_path, text=None, extra=()):
    payloads = [(MEMBER, (synthetic_text() if text is None else text).encode()),
                ("HEPdata_tables/Sivers.dat", b"Sivers explicitly synthetic and not a Collins input\n"), *extra]
    path = tmp_path / "synthetic.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, payload in payloads:
            entry = tarfile.TarInfo(name)
            entry.size = len(payload)
            archive.addfile(entry, io.BytesIO(payload))
    source = {"id": "synthetic_COMPASS", "publication": "explicitly synthetic parser fixture",
        "url": "https://example.invalid/synthetic-COMPASS", "version": "synthetic-v1",
        "sha256": file_hash(path), "filename": path.name,
        "adapter": {"kind": "compass2022_dat", "contract": "COMPASS2022-DATA-001", "member": MEMBER,
            "member_sha256": {name: hashlib.sha256(payload).hexdigest() for name, payload in payloads[:2]},
            "tables": [{"view_id": f"synthetic_h_{charge}_{projection}", "hadron": "h" + charge,
                        "projection": projection, "expected_rows": 2}
                       for charge in ("+", "-") for projection in ("x", "z", "pT")]}}
    return path, source


def candidate(views):
    return {"schema_version": 1, "scope": "experimental_candidate", "production_ready": False,
            "views": views}


def test_partial_hadron_identity_errors_units_and_convention_are_preserved(tmp_path):
    path, source = fixture(tmp_path)
    views = ex.parse_compass2022(path, source)
    assert ex.validate_candidate(candidate(views)) == {
        "views": 6, "points_across_overlapping_views": 12, "production_ready": False}
    assert [v["hadron"] for v in views] == ["h+"] * 3 + ["h-"] * 3
    for dataset in views:
        assert dataset["role"] == "review_required" and dataset["target_material"] == "6LiD"
        assert dataset["sample"] == "partial_2022_approximately_70_percent"
        assert dataset["convention"] == "published_COMPASS_sin(phi_h+phi_S-pi)_DNN_f_Pt_removed"
        assert dataset["depolarization_removed"] is True
        assert dataset["additional_depolarization_or_sign_conversion_applied"] is False
        assert dataset["overlap_with"] == ["COMPASS_deuteron_2022_complete"]
        scales = dataset["reported_scale_uncertainties"]
        assert [s["relative_sigma"] for s in scales] == [0.03, 0.02]
        assert all(s["inclusion_in_er_tot"] == "not_established" and
                   s["applied_to_likelihood"] is False for s in scales)
        item = dataset["points"][0]
        assert item["value"] == -0.017 and item["statistical_error"] == item["reported_total_error"] == 0.004
        assert item["systematic_components"] == []
        assert item["systematic_error_status"] == "not_decomposed_in_source"
        assert item["source_decimal_fields"]["Asym"] == "-0.0170"
        assert item["source_decimal_fields"]["er_tot"] == "0.0040"
        assert item["source_locator"]["archive_member"] == MEMBER
        assert file_hash(path) == source["sha256"]
    with pytest.raises(ValueError, match="scope"):
        validate_dataset(candidate(views))
    with pytest.raises(ValueError, match="overlap"):
        ex.choose_views(candidate(views), [v["id"] for v in views[:2]])
    complete = copy.deepcopy(views[0])
    complete["id"] = "synthetic_complete_2022"
    for item in complete["points"]:
        item["id"] += ":complete"
    with pytest.raises(ValueError, match="overlap"):
        ex.choose_views(candidate([views[0], complete]), [views[0]["id"], complete["id"]])


@pytest.mark.parametrize("before,after,message", [
    ("Collins asymmetry", "Sivers asymmetry", "source identity"),
    ("(~70%)", "(100%)", "source identity"),
    ("h^+; x-dependence", "pi^+; x-dependence", "section heading"),
    ("er_stat er_tot", "er_stat er_sys", "column identity"),
    ("-0.0170 0.0040 0.0040", "-0.0170 0.0040 0.0039", "uncertainty"),
    ("-0.0170 0.0040 0.0040", "-0.0170 0.0000 0.0040", "uncertainty"),
    ("-0.0170 0.0040 0.0040", "nan 0.0040 0.0040", "numeric cell"),
    ("[0.100;0.200] 0.150", "[0.100;0.200] 0.250", "mean outside"),
    ("[0.200;0.300]", "[0.210;0.300]", "noncontiguous"),
    ("6.200 2.400", "4.900 2.400", "published cuts"),
    ("6.200 2.400", "6.200 0.900", "published cuts"),
    ("6.200 2.400 -0.0170", "6.200 -0.0170", "numerical row"),
])
def test_changed_observable_sample_schema_or_bad_numerical_row_fails(tmp_path, before, after, message):
    text = synthetic_text()
    assert before in text
    path, source = fixture(tmp_path, text.replace(before, after, 1))
    with pytest.raises(ValueError, match=message):
        ex.parse_compass2022(path, source)


def test_projection_sections_row_counts_and_member_identity_are_exact(tmp_path):
    text = synthetic_text()
    path, source = fixture(tmp_path, text[:text.index("h^-; pT-dependence")])
    with pytest.raises(ValueError, match="sections changed"):
        ex.parse_compass2022(path, source)
    path, source = fixture(tmp_path)
    source["adapter"]["tables"][0]["expected_rows"] = 3
    with pytest.raises(ValueError, match="row count"):
        ex.parse_compass2022(path, source)
    source["adapter"]["tables"][0]["expected_rows"] = 2
    source["adapter"]["member"] = "HEPdata_tables/Sivers.dat"
    with pytest.raises(ValueError, match="member identity"):
        ex.parse_compass2022(path, source)


@pytest.mark.parametrize("extra", [[("../outside.dat", b"synthetic")],
                                    [("HEPdata_tables/unexpected.dat", b"synthetic")],
                                    [(MEMBER, b"synthetic duplicate")]])
def test_archive_traversal_unknown_or_duplicate_members_fail_without_extraction(tmp_path, extra):
    path, source = fixture(tmp_path, extra=extra)
    with pytest.raises(ValueError, match="unsafe|unexpected or duplicate"):
        ex.parse_compass2022(path, source)
    assert not (tmp_path / "outside.dat").exists()


def test_archive_member_hash_and_links_are_checked(tmp_path):
    path, source = fixture(tmp_path)
    source["adapter"]["member_sha256"][MEMBER] = "0" * 64
    with pytest.raises(ValueError, match="checksum"):
        ex.parse_compass2022(path, source)
    with tarfile.open(path, "w:gz") as archive:
        entry = tarfile.TarInfo(MEMBER)
        entry.type, entry.linkname = tarfile.SYMTYPE, "outside.dat"
        archive.addfile(entry)
    with pytest.raises(ValueError, match="unsafe"):
        ex.parse_compass2022(path, source)


@pytest.mark.parametrize("change", ["missing", "nan", "too_small", "boolean", "extra_systematic"])
def test_candidate_validation_preserves_stat_total_representation(tmp_path, change):
    path, source = fixture(tmp_path)
    data = candidate(ex.parse_compass2022(path, source))
    item = data["views"][0]["points"][0]
    if change == "missing":
        item.pop("reported_total_error")
    elif change == "extra_systematic":
        item["systematic_components"] = [{"name": "invented", "sigma": 0.001}]
    else:
        item["reported_total_error"] = {"nan": float("nan"), "too_small": 0.001, "boolean": True}[change]
    with pytest.raises(ValueError, match="statistical/total"):
        ex.validate_candidate(data)


@pytest.mark.parametrize("field", ["reported_total_error", "source_decimal_fields", "convention"])
def test_rehashed_compass_candidate_still_requires_exact_source_replay(tmp_path, monkeypatch, field):
    path, source = fixture(tmp_path)
    cache, output = tmp_path / "raw", tmp_path / "out"
    cache.mkdir()
    (cache / path.name).write_bytes(path.read_bytes())
    manifest = tmp_path / "manifest.json"
    atomic_json(manifest, {"sources": [source]})
    def decisions(data, paths):
        return {**ex.validate_candidate(data), "status": "SYNTHETIC_PARSER_DECISION_ONLY",
                "source_manifest_sha256": data["manifest_sha256"]}
    monkeypatch.setattr(ex, "_source_decisions", decisions)
    ex.run_intake(manifest, cache, output, offline=True)
    assert ex.verify_intake(manifest, cache, output)["verified"] is True
    data = read_json(output / "candidates.json")
    item = data["views"][0]["points"][0]
    if field == "reported_total_error":
        item[field] *= 2
    elif field == "source_decimal_fields":
        item[field]["er_tot"] = "0.0080"
    else:
        data["views"][0][field] = "unreviewed sign convention"
    atomic_json(output / "candidates.json", data)
    report = read_json(output / "intake-report.json")
    report["candidate_file_sha256"] = file_hash(output / "candidates.json")
    atomic_json(output / "intake-report.json", report)
    with pytest.raises(ValueError, match="locked source replay"):
        ex.verify_intake(manifest, cache, output)


def test_source_manifest_and_scientific_contract_bind_partial_hadron_tables():
    manifest = read_json(ROOT / "config/experimental-sources.json")
    source = next(s for s in manifest["sources"] if s["id"] == "COMPASS2022_partial_public_tables")
    assert source["sha256"] == "4287a629578739bbe00675ef4e563d2af7437eb45b4c84005649498b8a12af76"
    assert len(source["adapter"]["tables"]) == 6
    assert sum(d["expected_rows"] for d in source["adapter"]["tables"]) == 52
    contract = read_json(ROOT / "config/compass-2022-data-contract.json")
    assert contract["id"] == source["adapter"]["contract"]
    assert contract["scope"] == "experimental_candidate" and contract["production_ready"] is False
    assert contract["hadron_species"] == ["h+", "h-"]
    assert len(contract["numeric_field_order"]) == 11
    assert "deuteron_nuclear_prescription_not_certified" in contract["production_blockers"]
    assert any(s["id"] == "COMPASS2022_complete_preliminary_tables"
               for s in manifest["known_unavailable_sources"])
