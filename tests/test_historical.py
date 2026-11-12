"""Explicitly synthetic source-layout and admission tests; no experimental arrays."""
import copy
import hashlib
import io
import tarfile

import pytest
from openpyxl import Workbook

from collinsflow import historical as hist
from collinsflow.data import validate_dataset
from collinsflow.provenance import file_hash


def synthetic_tex():
    kine = [r"\begin{table*}", "x&Q2&y&z&P&W&Wp&a&b&c&d" + r"\\"]
    results = [r"\begin{table*}", "&x&Collins Moment&Collins Moment&Sivers Moment&Sivers Moment" + r"\\"]
    for i in range(1, 5):
        kine.append(f"{i / 10}&{i + 1}&0.6&0.4&0.2&3&2&0&0&0&0" + r"\\")
        results.append(f"n&{i / 10}&" + r"0.01\pm0.02\pm0.003 (0.004)&-0.02\pm0.03\pm0.004 (0.005)&0\pm1\pm0 (0)&0\pm1\pm0 (0)\\")
    kine.extend([r"\label{table:kine}", r"\end{table*}"])
    results.extend([r"\label{table:results}", r"\end{table*}"])
    body = "\n".join(kine + results)
    commented_duplicate = "\n".join("%" + line for line in results)
    return body + "\n" + commented_duplicate


def remap_contract():
    return {"HallA_remap": {"row_id_field": "Q2", "fields": {
        "Q2_GeV2": "x", "x": "y", "y": "z", "z": "pT", "pT_GeV": "value", "central_value": "stat_u",
        "statistical_error": "syst_u", "experimental_systematic": "fit", "model_systematic": "usr"}},
        "SIDIS_projection_cut": {"upper_GeV": .8}}


def legacy_neutron_row():
    return {"Q2": 1, "x": 2., "y": .1, "z": .6, "pT": .4, "value": .2,
            "stat_u": .01, "syst_u": .02, "fit": .003, "usr": .004}


def write_xlsx(path, rows):
    book = Workbook()
    sheet = book.active
    sheet.title = "synthetic"
    headers = list(rows[0])
    sheet.append(headers)
    for row in rows:
        sheet.append([row[k] for k in headers])
    book.save(path)
    book.close()
    return headers


def synthetic_source(path, headers, raw_rows, expected_points, projection="pt"):
    return {"id": "synthetic_history", "publication": "synthetic parser fixture", "version": "synthetic-v1",
        "url": "https://example.invalid/synthetic", "sha256": file_hash(path), "adapter": {
            "kind": "sidis", "experiment": "hermes", "target": "proton", "hadron": "pi+", "projection": projection,
            "source_obs": "AUTcollins", "input_unit": "fraction", "expected_points": expected_points,
            "expected_raw_rows": raw_rows, "expected_headers": headers, "sheet": "synthetic",
            "uncertainty_layout": "stat_exp", "depolarization_status": "raw_DNN_included_primary_verified"}}


def test_primary_tables_ignore_commented_duplicate_and_keep_neutron_harmonic():
    result = hist.parse_halla_tex(synthetic_tex())
    assert len(result) == 8
    row = result[("pi+", .1)]
    assert {k: row[k] for k in ("Q2_GeV2", "y", "z", "pT_GeV", "central_value", "statistical_error")} == {
        "Q2_GeV2": 2., "y": .6, "z": .4, "pT_GeV": .2, "central_value": .01, "statistical_error": .02}
    assert row["experimental_systematic"] == .003 and row["model_systematic"] == .004
    assert result[("pi-", .1)]["central_value"] == -.02


def test_primary_harmonic_and_incomplete_table_fail():
    with pytest.raises(ValueError, match="harmonic"):
        hist.parse_halla_tex(synthetic_tex().replace("Collins Moment", "Sivers Moment"))
    with pytest.raises(ValueError, match="incomplete"):
        hist.parse_halla_tex("\n".join(line for line in synthetic_tex().splitlines() if not line.startswith("n&0.4")))


@pytest.mark.parametrize("unsafe", ["../escape.tex", "/absolute.tex", "symlink"])
def test_primary_archive_does_not_extract_unsafe_members(tmp_path, unsafe):
    content = synthetic_tex().encode()
    path = tmp_path / "source.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        item = tarfile.TarInfo(unsafe)
        if unsafe == "symlink":
            item.type, item.linkname = tarfile.SYMTYPE, "../outside"
        else:
            item.size = len(content)
        archive.addfile(item, io.BytesIO(content) if item.isfile() else None)
    spec = {"sha256": file_hash(path), "member": unsafe, "member_sha256": hashlib.sha256(content).hexdigest()}
    with pytest.raises(ValueError, match="unsafe"):
        hist.halla_primary(path, spec)


def test_primary_archive_hash_and_member_are_bound(tmp_path):
    content = synthetic_tex().encode()
    path = tmp_path / "source.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        item = tarfile.TarInfo("synthetic.tex")
        item.size = len(content)
        archive.addfile(item, io.BytesIO(content))
    spec = {"sha256": file_hash(path), "member": "synthetic.tex", "member_sha256": hashlib.sha256(content).hexdigest()}
    assert len(hist.halla_primary(path, spec)) == 8
    spec["member_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="changed primary TeX"):
        hist.halla_primary(path, spec)


def test_legacy_remap_requires_exact_primary_values_and_preserves_components():
    adapter = {"uncertainty_layout": "stat_exp_model_remap", "hadron": "pi+"}
    row = legacy_neutron_row()
    primary = hist.parse_halla_tex(synthetic_tex())
    corrected = hist.sidis_point(row, adapter, remap_contract(), primary)
    assert corrected["averages"] == {"x": .1, "y": .6, "z": .4, "Q2_GeV2": 2., "pT_GeV": .2}
    assert corrected["central_value"] == .01 and corrected["statistical_error"] == .02
    assert corrected["systematic_components"] == {"experimental": .003, "model": .004}
    assert row == legacy_neutron_row()
    row["fit"] = .004
    with pytest.raises(ValueError, match="disagrees with primary"):
        hist.sidis_point(row, adapter, remap_contract(), primary)
    row = legacy_neutron_row()
    row["Q2"] = True
    with pytest.raises(ValueError, match="row-id"):
        hist.sidis_point(row, adapter, remap_contract(), primary)


def test_strict_pT_cut_keeps_source_locators_and_does_not_cut_other_projections(tmp_path):
    rows = [{"Q2": 2., "x": .1, "y": .6, "z": .4, "pT": p, "obs": "AUTcollins", "value": .01,
             "stat_u": .02, "syst_u": .003, "target": "proton", "hadron": "pi+", "col": "HERMES",
             "dependence": "pt"} for p in (.79, .8, .81)]
    path = tmp_path / "synthetic.xlsx"
    headers = write_xlsx(path, rows)
    spec = synthetic_source(path, headers, 3, 1)
    view = hist.parse_historical_xlsx(path, spec, remap_contract(), {})
    assert len(view["points"]) == 1 and [v["reported_pT_GeV"] for v in view["excluded"]] == [.8, .81]
    assert [v["locator"]["excel_row"] for v in view["excluded"]] == [3, 4]
    for row in rows:
        row["dependence"] = "x"
    headers = write_xlsx(path, rows)
    view = hist.parse_historical_xlsx(path, synthetic_source(path, headers, 3, 3, "x"), remap_contract(), {})
    assert len(view["points"]) == 3 and not view["excluded"]


def test_sia_percent_conversion_does_not_invent_mean_P_from_mean_qT():
    row = {"z1": .2, "z2": .3, "Q": 10., "qT": 1.5, "S2/1+C2": .7, "value": 1.2, "stat_u": .3, "sys_u": .4}
    adapter = {"input_unit": "percent", "error_columns": {"statistical": "stat_u", "reported_systematic": "sys_u"},
               "angular_column": "S2/1+C2"}
    p = hist.sia_point(row, adapter)
    assert (p["central_value"], p["statistical_error"], p["systematic_components"]["reported_systematic"]) == pytest.approx((.012, .003, .004))
    assert p["averages"]["qT_GeV"] == 1.5 and "pT_GeV" not in p["averages"]
    assert p["angular_source"]["average_type"] == "publication_review_pending"
    row["syt_u"] = row.pop("sys_u")
    with pytest.raises(ValueError, match="declared.*column"):
        hist.sia_point(row, adapter)
    adapter["error_columns"]["reported_systematic"] = "syt_u"
    assert hist.sia_point(row, adapter)["systematic_components"]["reported_systematic"] == .004


def test_source_harmonic_and_checksum_are_not_inferred(tmp_path):
    row = {"Q2": 2., "x": .1, "y": .6, "z": .4, "pT": .3, "obs": "AUTsivers", "value": .01,
           "stat_u": .02, "syst_u": .003, "target": "proton", "hadron": "pi+", "col": "HERMES", "dependence": "pt"}
    path = tmp_path / "synthetic.xlsx"
    headers = write_xlsx(path, [row])
    spec = synthetic_source(path, headers, 1, 1)
    with pytest.raises(ValueError, match="harmonic"):
        hist.parse_historical_xlsx(path, spec, remap_contract(), {})
    spec["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checksum"):
        hist.parse_historical_xlsx(path, spec, remap_contract(), {})


@pytest.mark.parametrize("key", ["production_ready", "physical_extraction", "historical_fit_ready"])
def test_historical_candidates_cannot_claim_admission(key):
    data = {"schema_version": 1, "scope": hist.SCOPE, "production_ready": False,
            "physical_extraction": False, "historical_fit_ready": False, "views": []}
    with pytest.raises(ValueError, match="scope"):
        validate_dataset(data)
    changed = copy.deepcopy(data)
    changed[key] = True
    with pytest.raises(ValueError, match="self-certify"):
        hist.validate_historical_candidate(changed, {})
