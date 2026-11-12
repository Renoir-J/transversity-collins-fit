"""Synthetic parser and admission regressions; no experimental arrays are vendored."""
import copy
import io
import tarfile

import pytest
import yaml

from collinsflow import experimental as ex
from collinsflow.data import validate_dataset
from collinsflow.provenance import atomic_json, file_hash, read_json


def source(**kwargs):
    return {"id": "synthetic_source", "publication": "synthetic parser fixture", "sha256": "0" * 64,
            "url": "https://example.invalid/synthetic", "version": "synthetic-v1", **kwargs}


def candidate():
    dataset = ex.view("synthetic_view", source(), "SIDIS", "pi+", "test", "same_sample", "test")
    dataset["points"] = [ex.point("synthetic:0", 0.01, 0.002, 0.001, bins={"x": [0.1, 0.2]},
                                  averages={"x": 0.15}, locator={"synthetic_row": 1})]
    return {"schema_version": 1, "scope": "experimental_candidate", "views": [dataset]}


@pytest.mark.parametrize("unit,cell,expected", [
    ("percent", r"1.25 \pm 0.20 \pm 0.10", (0.0125, 0.002, 0.001)),
    ("fraction", r"0.0125 \pm 0.0020 \pm 0.0010", (0.0125, 0.002, 0.001)),
])
def test_unit_conversion_preserves_stat_and_sys_separately(unit, cell, expected):
    assert ex.measurement(cell, unit) == pytest.approx(expected)


@pytest.mark.parametrize("cell", ["1.0 (2)", "NaN", "1+2", "1.0 or True", "1.0 +/- 0.1"])
def test_ambiguous_numeric_cells_are_not_guessed(cell):
    with pytest.raises(ValueError):
        ex.scalar(cell)


def test_math_annotation_is_read_once_and_merged_cells_are_rejected(tmp_path):
    path = tmp_path / "synthetic.html"
    path.write_text('<figure id="synthetic"><table><tr><td><math alttext="1.25">'
                    '<mn>1.25</mn><annotation>1.25</annotation></math></td></tr></table></figure>')
    assert ex.html_rows(path, "synthetic") == [["1.25"]]
    path.write_text('<figure id="synthetic"><table><tr><td colspan="2">1</td></tr></table></figure>')
    with pytest.raises(ValueError, match="merged"):
        ex.html_rows(path, "synthetic")


def test_candidates_are_not_fit_data_and_cannot_self_certify():
    data = candidate()
    assert ex.validate_candidate(data)["production_ready"] is False
    with pytest.raises(ValueError, match="scope"):
        validate_dataset(data)
    data["views"][0]["production_admission"]["status"] = "PASS"
    with pytest.raises(ValueError, match="certify"):
        ex.validate_candidate(data)


def test_stat_and_systematic_errors_are_validated():
    data = candidate()
    data["views"][0]["points"][0]["statistical_error"] = 0
    with pytest.raises(ValueError, match="statistical"):
        ex.validate_candidate(data)
    data = candidate()
    data["views"][0]["points"][0]["systematic_components"][0]["sigma"] = -0.01
    with pytest.raises(ValueError, match="systematic"):
        ex.validate_candidate(data)


def test_overlapping_projection_and_double_ratio_selection_stops():
    data = candidate()
    other = copy.deepcopy(data["views"][0])
    other["id"] = "synthetic_alternate_projection"
    other["points"][0]["id"] = "synthetic:1"
    data["views"].append(other)
    with pytest.raises(ValueError, match="overlap"):
        ex.choose_views(data, [v["id"] for v in data["views"]])
    assert len(ex.choose_views(data, ["synthetic_view"])) == 1


def test_locked_cache_and_seed_are_checked_without_network(tmp_path):
    seed = tmp_path / "seed.txt"
    seed.write_text("explicitly synthetic source\n")
    spec = source(filename="synthetic.txt", sha256=file_hash(seed))
    cached = ex.acquire(spec, tmp_path / "cache", seed, offline=True)
    assert cached.read_bytes() == seed.read_bytes()
    cached.write_text("changed\n")
    with pytest.raises(ValueError, match="checksum"):
        ex.acquire(spec, tmp_path / "cache", offline=True)
    spec["filename"] = "../escape.txt"
    with pytest.raises(ValueError, match="basename"):
        ex.acquire(spec, tmp_path / "cache", offline=True)


def test_besiii_exchange_fold_is_a_union_not_a_single_box(tmp_path):
    header = '<tr>' + ''.join(f'<th>c{i}</th>' for i in range(7)) + '</tr>'
    zrow = ['[0.2,0.3][0.3,0.5]', '0.32', '0.32', '0.4', '0.6',
            r'0.02 \pm 0.003 \pm 0.001', r'0.01 \pm 0.002 \pm 0.001']
    ptrow = ['[0.2,0.3]', '0.25', '0.32', '0.33', '0.6',
             r'0.03 \pm 0.003 \pm 0.001', r'0.02 \pm 0.002 \pm 0.001']
    def row(values):
        return '<tr>' + ''.join(f'<td>{cell}</td>' for cell in values) + '</tr>'
    path = tmp_path / "synthetic.html"
    path.write_text('<figure id="S0.T1"><table>' + header + row(zrow) * 6 + header + row(ptrow) * 5 + '</table></figure>')
    views = ex.parse_besiii(path, source(version="1507.06824v2"))
    folded = views[0]["points"][0]
    assert folded["bins"] == {}
    assert folded["bin_union"] == [{"z1": [0.2, 0.3], "z2": [0.3, 0.5]},
                                   {"z1": [0.3, 0.5], "z2": [0.2, 0.3]}]
    with pytest.raises(ValueError, match="final"):
        ex.parse_besiii(path, source(version="1507.06824v1"))


def test_hermes_source_correction_is_explicit_and_requires_original_value(monkeypatch):
    record = {"obs": "AUTcollins", "hadron": "pi+", "target": "p", "col": "HERMES", "depol": False,
              "Q2": 2.0, "x": 0.1, "y": 0.4, "z": 0.6, "pT": 0.3, "zmin": 0.049, "zmax": 0.7,
              "value": 0.02, "stat_u": 0.003, "syst_u": 0.001}
    monkeypatch.setattr(ex, "xlsx_records", lambda *args: [dict(record)])
    adapter = {"view_id": "synthetic_HERMES", "hadron": "pi+", "dimension": 1, "projection": "z",
               "depolarization_column": "depol", "sheet": "synthetic", "bins": {"z": ["zmin", "zmax"]},
               "expected_rows": 1, "corrections": [{"source_excel_row": 2, "column": "zmin", "old_value": 0.049,
                                                    "new_value": 0.49, "reference": "synthetic fixture"}]}
    result = ex.parse_hermes("unused", source(adapter=adapter))[0]
    assert result["points"][0]["bins"]["z"] == [0.49, 0.7]
    assert result["points"][0]["source_corrections"][0]["old_value"] == 0.049
    assert record["zmin"] == 0.049
    record["zmin"] = 0.49
    with pytest.raises(ValueError, match="no longer matches"):
        ex.parse_hermes("unused", source(adapter=adapter))


def test_real_source_manifest_locks_versions_and_corrections():
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    manifest = read_json(root / "config/experimental-sources.json")
    sources = {s["id"]: s for s in manifest["sources"]}
    assert len(sources) == len(manifest["sources"]) == 23
    assert sources["arxiv_1507.06824v2"]["version"] == "1507.06824v2"
    corrections = sources["HERMES_JAM_3709"]["adapter"]["corrections"]
    assert [c["source_excel_row"] for c in corrections] == [14, 15, 16, 17]
    assert all(c["old_value"] == 0.049 and c["new_value"] == 0.49 for c in corrections)
    assert sources["HERMES_JAM_3710"]["adapter"]["depolarization_column"] == "depole"
    for identifier, tables, rows in [("STAR200_HEPData_v1", 44, 328), ("STAR510_HEPData_v1", 12, 76)]:
        definitions = sources[identifier]["adapter"]["tables"]
        assert len(definitions) == tables
        assert sum(d["expected_rows"] for d in definitions) == rows


def test_missing_locked_source_does_not_create_a_candidate_report(tmp_path):
    manifest = tmp_path / "manifest.json"
    atomic_json(manifest, {"sources": [source(filename="missing.html")]})
    with pytest.raises(FileNotFoundError):
        ex.run_intake(manifest, tmp_path / "cache", tmp_path / "outputs", offline=True)
    assert not (tmp_path / "outputs/candidates.json").exists()


def synthetic_archive(path, members):
    with tarfile.open(path, "w:gz") as archive:
        for name, payload in members:
            entry = tarfile.TarInfo(name)
            entry.size = len(payload)
            archive.addfile(entry, io.BytesIO(payload))
    return path


def test_star_plot_uncertainties_do_not_become_physical_bins_or_means(tmp_path):
    filename = r"figure2_$\pi^+$.yaml"
    table = {"independent_variables": [{"header": {"name": "$x_T$"},
                                      "values": [{"low": 0.09, "high": 0.11}]}],
             "dependent_variables": [{"header": {"name": r"$A_{UT}^{sin(\phi_{S}-\phi_{H})}$"},
                 "qualifiers": [{"name": "synthetic sample", "value": "positive hemisphere"}],
                 "values": [{"value": 0.01, "errors": [{"label": "stat", "symerror": 0.002},
                                                          {"label": "sys", "symerror": 0.001}]}]}]}
    metadata = [{"hepdata_doi": "synthetic-doi", "data_license": {"name": "synthetic-license"}},
                {"data_file": filename, "table_doi": "synthetic-table-doi", "name": "Synthetic table",
                 "description": "Horizontal box represents systematic uncertainty in xT."}]
    adapter = {"doi": "synthetic-doi", "energy_GeV": 510, "tables": [
        {"data_file": filename, "view_id": "synthetic_STAR", "hadron": "pi+", "expected_rows": 1,
         "coordinate": "xT", "sample": "synthetic sample"}]}
    def write():
        return synthetic_archive(tmp_path / "synthetic.tar.gz", [
            ("submission.yaml", yaml.safe_dump_all(metadata).encode()),
            (filename, yaml.safe_dump(table).encode())])
    views = ex.parse_star_hepdata(write(), source(adapter=adapter))
    data = {"schema_version": 1, "scope": "experimental_candidate", "views": views}
    assert ex.validate_candidate(data)["production_ready"] is False
    item = views[0]["points"][0]
    assert item["bins"] == {} and item["averages"] == {}
    assert item["source_intervals"]["xT"] == {
        "low": 0.09, "high": 0.11, "meaning": "reported_plot_systematic_interval_not_physical_binning"}
    assert item["source_coordinate_record"] == {"low": 0.09, "high": 0.11}
    assert views[0]["source_qualifiers"] == table["dependent_variables"][0]["qualifiers"]
    table["dependent_variables"][0]["header"]["name"] = r"$A_{UT}^{sin(\phi_{S}-2\phi_{H})}$"
    with pytest.raises(ValueError, match="Collins harmonic"):
        ex.parse_star_hepdata(write(), source(adapter=adapter))


@pytest.mark.parametrize("name", ["../outside.yaml", "/absolute.yaml"])
def test_archive_traversal_is_rejected_without_extracting(tmp_path, name):
    path = synthetic_archive(tmp_path / "synthetic.tar.gz", [(name, b"synthetic: true\n")])
    with pytest.raises(ValueError, match="unsafe"):
        ex.hepdata_documents(path)
    assert not (tmp_path / "outside.yaml").exists()


def test_archive_links_and_duplicate_filenames_are_rejected(tmp_path):
    path = tmp_path / "synthetic.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        entry = tarfile.TarInfo("table.yaml")
        entry.type = tarfile.SYMTYPE
        entry.linkname = "outside.yaml"
        archive.addfile(entry)
    with pytest.raises(ValueError, match="unsafe"):
        ex.hepdata_documents(path)
    synthetic_archive(path, [("a/table.yaml", b"synthetic: true"), ("b/table.yaml", b"synthetic: false")])
    with pytest.raises(ValueError, match="duplicate"):
        ex.hepdata_documents(path)


def test_bin_unions_and_plot_intervals_are_checked():
    data = candidate()
    item = data["views"][0]["points"][0]
    item["bin_union"] = [{"z1": [0.3, 0.2]}]
    with pytest.raises(ValueError, match="experimental bin"):
        ex.validate_candidate(data)
    item["bin_union"] = []
    item["source_intervals"] = {"xT": {"low": float("nan"), "high": 0.2,
        "meaning": "reported_plot_systematic_interval_not_physical_binning"}}
    with pytest.raises(ValueError, match="plot-systematic"):
        ex.validate_candidate(data)


@pytest.fixture
def cached_intake(tmp_path, monkeypatch):
    manifest, cache, output = tmp_path / "manifest.json", tmp_path / "raw", tmp_path / "output"
    cache.mkdir()
    raw, lattice = cache / "synthetic-measurement.json", cache / "synthetic-lattice.json"
    atomic_json(raw, {"value": 0.01, "stat": 0.002, "sys": 0.001})
    atomic_json(lattice, [{"value": 0.936, "stat_u": 0.025, "Q": 2}])
    atomic_json(manifest, {"sources": [
        source(id="arxiv_2007.07755v1", filename=raw.name, sha256=file_hash(raw),
               version="explicitly-synthetic-v1", adapter={"kind": "hermes_xlsx"}),
        source(id="JAM_lattice_2000", filename=lattice.name, sha256=file_hash(lattice),
               version="explicitly-synthetic-v1")]})

    def synthetic_parser(path, lock):
        numbers = read_json(path)
        dataset = ex.view("synthetic_view", lock, "SIDIS", "pi+", "synthetic_Collins",
                          "synthetic_sample", "synthetic_sin(phi_h+phi_S)", role="primary_candidate")
        dataset["points"] = [ex.point("synthetic:0", numbers["value"], numbers["stat"], numbers["sys"],
            bins={"x": [0.1, 0.2]}, averages={"x": 0.15}, locator={"explicitly_synthetic_row": 1})]
        return [dataset]

    monkeypatch.setattr(ex, "parse_hermes", synthetic_parser)
    monkeypatch.setattr(ex, "xlsx_records", lambda path: read_json(path))
    monkeypatch.setattr(ex, "crosscheck_besiii", lambda *args: [])
    monkeypatch.setattr(ex, "crosscheck_hermes_grid", lambda *args: [
        {"view": "synthetic_view", "points": 1, "passed": True}])
    monkeypatch.setattr(ex, "source_fingerprint", lambda *args: "synthetic-code-fingerprint")
    ex.run_intake(manifest, cache, output, offline=True)
    return manifest, cache, output, raw


def test_cached_intake_rechecks_artifacts_raw_sources_and_parser(cached_intake, monkeypatch):
    manifest, cache, output, raw = cached_intake
    frozen = (output / "candidates.json").read_bytes()
    original_raw = raw.read_bytes()
    assert ex.verify_intake(manifest, cache, output)["production_ready"] is False
    (output / "candidates.json").write_bytes(frozen + b" ")
    with pytest.raises(ValueError, match="artifact was changed"):
        ex.verify_intake(manifest, cache, output)
    (output / "candidates.json").write_bytes(frozen)
    raw.write_text("changed source\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        ex.verify_intake(manifest, cache, output)
    raw.write_bytes(original_raw)
    monkeypatch.setattr(ex, "source_fingerprint", lambda *args: "changed-synthetic-parser")
    with pytest.raises(ValueError, match="stale experimental parser"):
        ex.verify_intake(manifest, cache, output)


@pytest.mark.parametrize("mutation", [
    "value", "stat", "systematic", "mean", "bin", "harmonic", "locator", "overlap", "role",
    "source_reference", "production_claim"])
def test_cached_intake_rejects_rehashed_candidate_changes(cached_intake, mutation):
    manifest, cache, output, _ = cached_intake
    path, report_path = output / "candidates.json", output / "intake-report.json"
    data, report = read_json(path), read_json(report_path)
    dataset, item = data["views"][0], data["views"][0]["points"][0]
    if mutation == "value":
        item["value"] += 0.001
    elif mutation == "stat":
        item["statistical_error"] *= 2
    elif mutation == "systematic":
        item["systematic_components"][0]["sigma"] *= 2
    elif mutation == "mean":
        item["averages"]["x"] = 0.16
    elif mutation == "bin":
        item["bins"]["x"][0] = 0.11
    elif mutation == "harmonic":
        dataset["convention"] = "synthetic_sin(phi_h-phi_S)"
    elif mutation == "locator":
        item["source_locator"]["explicitly_synthetic_row"] += 1
    elif mutation == "overlap":
        dataset["overlap_group"] = "synthetic_forged_disjoint_sample"
    elif mutation == "role":
        dataset["role"] = "alternative"
        report["primary_candidate_views"] = []
    elif mutation == "source_reference":
        data["sources"][0]["version"] = "synthetic_forged_v2"
    elif mutation == "production_claim":
        data["production_ready"] = True
    atomic_json(path, data)
    report["candidate_file_sha256"] = file_hash(path)
    atomic_json(report_path, report)
    with pytest.raises(ValueError, match="candidate differs from locked source replay"):
        ex.verify_intake(manifest, cache, output)


@pytest.mark.parametrize("mutation", [
    "counts", "primary_decision", "crosscheck_count", "crosscheck_status", "lattice", "STAR",
    "source_digest", "status", "production_claim"])
def test_cached_intake_rejects_changed_source_decisions(cached_intake, mutation):
    manifest, cache, output, _ = cached_intake
    path = output / "intake-report.json"
    report = read_json(path)
    if mutation == "counts":
        report["views"] += 1
    elif mutation == "primary_decision":
        report["primary_candidate_views"] = []
    elif mutation == "crosscheck_count":
        report["HERMES_3D_primary_binning"][0]["points"] += 1
    elif mutation == "crosscheck_status":
        report["HERMES_3D_primary_binning"][0]["passed"] = False
    elif mutation == "lattice":
        report["lattice_JAM_crosscheck"]["value"] += 0.1
    elif mutation == "STAR":
        report["STAR_source_review"]["hepdata_views"] += 1
    elif mutation == "source_digest":
        report["source_manifest_digest"] = "0" * 64
    elif mutation == "status":
        report["status"] = "PASS"
    elif mutation == "production_claim":
        report["physical_extraction"] = True
    atomic_json(path, report)
    with pytest.raises(ValueError, match="report differs from locked source replay"):
        ex.verify_intake(manifest, cache, output)


@pytest.mark.parametrize("artifact,field", [
    ("candidate", "production_ready"), ("candidate", "schema_version"),
    ("report", "production_ready"), ("report", "passed"), ("report", "views")])
def test_cached_intake_rejects_scalar_type_aliases(cached_intake, artifact, field):
    manifest, cache, output, _ = cached_intake
    path, report_path = output / "candidates.json", output / "intake-report.json"
    data, report = read_json(path), read_json(report_path)
    target = data if artifact == "candidate" else report
    if field == "production_ready":
        target[field] = 0
    elif field == "schema_version":
        target[field] = True
    elif field == "passed":
        target["HERMES_3D_primary_binning"][0][field] = 1
    else:
        target[field] = float(target[field])
    atomic_json(path, data)
    report["candidate_file_sha256"] = file_hash(path)
    atomic_json(report_path, report)
    with pytest.raises(ValueError, match="differs from locked source replay"):
        ex.verify_intake(manifest, cache, output)


@pytest.mark.parametrize("mutation", ["raw", "manifest", "parser", "candidate", "report"])
def test_cached_intake_rejects_snapshot_changes_during_replay(cached_intake, monkeypatch, mutation):
    manifest, cache, output, raw = cached_intake
    parser = ex.parse_hermes

    def change_after_reading(path, lock):
        result = parser(path, lock)
        if mutation == "raw":
            raw.write_text("changed explicitly synthetic raw source\n")
        elif mutation == "manifest":
            manifest.write_bytes(manifest.read_bytes() + b" ")
        elif mutation == "parser":
            monkeypatch.setattr(ex, "source_fingerprint", lambda *args: "changed-synthetic-parser")
        else:
            name = "candidates.json" if mutation == "candidate" else "intake-report.json"
            target = output / name
            target.write_bytes(target.read_bytes() + b" ")
        return result

    monkeypatch.setattr(ex, "parse_hermes", change_after_reading)
    with pytest.raises(ValueError, match="changed during|checksum mismatch"):
        ex.verify_intake(manifest, cache, output)
