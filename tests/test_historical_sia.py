"""Synthetic BaBar-shaped text/rows only; no experimental values or admission."""
import copy
from decimal import Decimal
import hashlib
import io
import tarfile

import pytest
from openpyxl import Workbook

from collinsflow import historical as hist
from collinsflow.provenance import file_hash


def synthetic_babar_tex():
    lines = [r"\def\aul{\ensuremath{{ A_{0}^{UL}}}}", r"\def\auc{\ensuremath{{ A_{0}^{UC}}}}"]
    for ratio, macro, measurements in (("ul", "aul", ((-1.25, .08, .03), (2.5, .12, .04))),
                                        ("uc", "auc", ((-.55, .04, .02), (.9, .06, .02)))):
        header = [r"\multirow{2}{*}{$z_1$}", r"\multirow{2}{*}{$\langle z_1\rangle$}",
            r"$\langle p_{t1}\rangle$", r"\multirow{2}{*}{$z_2$}", r"\multirow{2}{*}{$\langle z_2\rangle$}",
            r"$\langle p_{t2}\rangle$", r"$\langle p_{t0}\rangle$",
            r"\multirow{2}{*}{$\left\langle\frac{s_{th}}{1+c_{th}}\right\rangle$}",
            r"\multirow{2}{*}{$\left\langle\frac{s_{2}}{1+c_{2}}\right\rangle$}",
            "\\multirow{2}{*}{\\"+macro+r"Th (\%)}", "", "\\multirow{2}{*}{\\"+macro+r" (\%)}"]
        lines.extend([r"\begin{table*}", "Synthetic amplitudes in percent; statistical and systematic errors.",
                      " & ".join(header)+r"\\"])
        for i, (a, b, z1, z2, p, angular) in enumerate((
                ("[0.2,0.3]", "[0.3,0.4]", .25, .35, .1, .62),
                ("[0.3,0.4]", "[0.2,0.3]", .35, .25, .2, .71))):
            measured = " $\\pm$ ".join(map(str, measurements[i]))
            row = [f"${a}$", str(z1), "7.1", f"${b}$", str(z2), "7.2", str(p), ".99", str(angular),
                   r"99 $\pm$ 1 $\pm$ 2", "", measured]
            lines.append(" & ".join(row)+r"\\")
        lines.extend(["\\label{tab:totZ"+ratio+"}", r"\end{table*}"])
    lines.extend([r"\begin{table*}", r"Synthetic caption mentions $p_{t0}$ but is not the RF0 header.",
        r"$p_{t1}$ & $\langle p_{t1}\rangle$ & RF12 header\\",
        r"$[0.,0.2]$ & .1 & .24 & $[0.,0.2]$ & .1 & .36 & .99 & 99 $\pm$ 1 $\pm$ 2 & & 99 $\pm$ 1 $\pm$ 2\\",
        r"$p_{t0}$ & $\langle p_{t0}\rangle$ & \multirow{2}{*}{$\langle z_1\rangle$} & & &",
        r"\multirow{2}{*}{$\langle z_2\rangle$} & \multirow{2}{*}{$\left\langle\frac{\sin^2\theta_{2}}{1+\cos^2\theta_{2}}\right\rangle$}",
        r"& \multirow{2}{*}{\aul} & & \multirow{2}{*}{\auc}\\",
        r"$[0.,0.2]$ & .1 & .24 & & & .36 & .62 & -.8 $\pm$ .1 $\pm$ .03 & & -.3 $\pm$ .05 $\pm$ .02\\",
        r"$[0.2,0.4]$ & .3 & .26 & & & .34 & .68 & 1.8 $\pm$ .2 $\pm$ .04 & & .7 $\pm$ .1 $\pm$ .02\\",
        r"\label{tab:totPt}", r"\end{table*}"])
    return "\n".join(lines)


def primary_fixture(tmp_path, text=None):
    content = (synthetic_babar_tex() if text is None else text).encode()
    path = tmp_path/"synthetic-primary.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        item = tarfile.TarInfo("synthetic.tex")
        item.size = len(content)
        archive.addfile(item, io.BytesIO(content))
    spec = {"id": "synthetic_babar", "publication": "synthetic-v1", "sha256": file_hash(path),
        "member": "synthetic.tex", "member_sha256": hashlib.sha256(content).hexdigest(),
        "measurement_unit": "percent", "frame": "hadron2",
        "table_rows": {"tab:totZul": 2, "tab:totZuc": 2, "tab:totPt": 2}}
    return path, spec


def secondary_fixture(tmp_path, systematic="sys_u", mutation=None):
    rows = [dict(z1min=.2, z1max=.3, z2min=.3, z2max=.4, z1=.25, z2=.35, pT=.1,
                 Q=7., **{"S2/1+C2": .62}, value=-1.25, stat_u=.08, **{systematic: .03}),
            dict(z1min=.3, z1max=.4, z2min=.2, z2max=.3, z1=.35, z2=.25, pT=.2,
                 Q=7., **{"S2/1+C2": .71}, value=2.5, stat_u=.12, **{systematic: .04})]
    for row in rows:
        row.update(col="babar", obs="synthetic-UL-RF0", hadron1="pi", hadron2="pi")
    if mutation:
        mutation(rows)
    path = tmp_path/"synthetic-legacy.xlsx"
    book = Workbook()
    sheet = book.active
    headers = list(rows[0])
    sheet.append(headers)
    for row in rows:
        sheet.append([row[k] for k in headers])
    book.save(path)
    book.close()
    spec = {"id": "synthetic_legacy", "publication": "synthetic", "version": "synthetic-v1",
        "url": "https://example.invalid/synthetic", "sha256": file_hash(path), "adapter": {
            "kind": "sia", "experiment": "babar", "ratio": "UL", "projection": "z1_z2", "frame": "hadron2",
            "primary_crosscheck_id": "synthetic_babar", "input_unit": "percent", "source_obs": "synthetic-UL-RF0",
            "sheet": sheet.title, "expected_headers": headers, "expected_raw_rows": 2, "expected_points": 2,
            "error_columns": {"statistical": "stat_u", "reported_systematic": systematic}, "angular_column": "S2/1+C2"}}
    return path, spec


def synthetic_contract():
    return {"id": "SYNTHETIC-SOURCE-ONLY", "BaBar_primary": {"publication": "synthetic-v1",
        "legacy_publication": "synthetic", "averaging_measure": "declared synthetic bin-average measure; not an experiment"},
        "expected": {"views": 1, "points": 2, "SIDIS_views": 0, "SIDIS_points": 0, "SIA_views": 1, "SIA_points": 2}}


def candidate_view(tmp_path, systematic="sys_u", mutation=None):
    archive, primary_spec = primary_fixture(tmp_path)
    primary = hist.babar_primary(archive, primary_spec)
    path, source = secondary_fixture(tmp_path, systematic, mutation)
    view = hist.parse_historical_xlsx(path, source, synthetic_contract(), {}, babar_tables={primary["id"]: primary})
    return view, path, source, primary


def test_rf0_selection_keeps_signed_amplitudes_and_ignores_rf12_and_caption_tokens(tmp_path):
    path, spec = primary_fixture(tmp_path)
    primary = hist.babar_primary(path, spec)
    assert set(primary["views"]) == {(r, p) for r in ("UL", "UC") for p in ("z1_z2", "pT")}
    ul = primary["views"][("UL", "z1_z2")]["points"][0]
    assert ul["fields"]["value"] == Decimal("-1.25")
    assert ul["fields"]["pT"] == Decimal("0.1") and ul["fields"]["S2/1+C2"] == Decimal("0.62")
    assert primary["views"][("UL", "pT")]["points"][0]["fields"]["value"] == Decimal("-.8")
    assert primary["views"][("UC", "pT")]["points"][0]["fields"]["value"] == Decimal("-.3")
    assert sum(len(p["fields"]) for v in primary["views"].values() for p in v["points"]) == 80


@pytest.mark.parametrize("systematic", ["sys_u", "syt_u"])
def test_primary_crosscheck_preserves_existing_numbers_and_records_limited_provenance(tmp_path, systematic):
    view, path, source, primary = candidate_view(tmp_path, systematic)
    assert file_hash(path) == source["sha256"]
    assert [p["central_value"] for p in view["points"]] == [-.0125, .025]
    assert [p["averages"]["Q2_GeV2"] for p in view["points"]] == [49., 49.]
    assert view["primary_numeric_crosscheck"]["numeric_cells"] == 22
    assert systematic in view["primary_numeric_crosscheck"]["compared_columns"]
    assert "Q" not in view["primary_numeric_crosscheck"]["compared_columns"]
    for point in view["points"]:
        assert point["angular_source"]["average_type"] == "mean_of_angular_ratio"
        assert point["angular_source"]["measure_status"] == "joint_acceptance_unresolved"
        assert point["primary_locator"]["publication"] == primary["publication"]
        assert point["primary_locator"]["source_line"] > 0
        assert point["uncertainty_source_columns"]["reported_systematic"] == systematic
    assert view["points"][0]["primary_bins"] == {"z1": [.2, .3], "z2": [.3, .4]}
    assert view["points"][0]["angular_source"]["value"] == .62
    assert any("acceptance" in b for b in view["blockers"])
    assert view["production_ready"] is view["historical_fit_ready"] is False


@pytest.mark.parametrize("column", ["value", "stat_u", "sys_u", "z1", "z2", "pT", "S2/1+C2", "z1min", "z2max"])
def test_every_declared_primary_field_requires_exact_agreement(tmp_path, column):
    with pytest.raises(ValueError, match="disagrees with primary"):
        candidate_view(tmp_path, mutation=lambda rows: rows[0].update({column: rows[0][column]+.001}))


def test_extra_decimal_is_not_accepted_by_rounding_to_printed_precision(tmp_path):
    with pytest.raises(ValueError, match="disagrees with primary"):
        candidate_view(tmp_path, mutation=lambda rows: rows[0].update(value=-1.2500000001))


@pytest.mark.parametrize("field,value", [("frame", "thrust"), ("ratio", "UC"), ("projection", "pT")])
def test_equal_numeric_arrays_cannot_authorize_wrong_primary_identity(tmp_path, field, value):
    archive, spec = primary_fixture(tmp_path)
    primary = hist.babar_primary(archive, spec)
    primary["views"][("UL", "z1_z2")][field] = value
    path, source = secondary_fixture(tmp_path)
    with pytest.raises(ValueError, match="frame/ratio/projection"):
        hist.parse_historical_xlsx(path, source, synthetic_contract(), {}, babar_tables={spec["id"]: primary})


@pytest.mark.parametrize("old,new", [
    (r"A_{0}^{UL}", r"A_{12}^{UL}"),
    (r"\aul (\%)", r"\aulTh (\%)"),
    (r"\frac{s_{2}}{1+c_{2}}", r"\frac{s_{th}}{1+c_{th}}"),
    (r"\sin^2\theta_{2}", r"\sin^2\theta_{th}"),
    (r"\aul (\%)", r"\aul (fraction)"),
])
def test_primary_header_and_macro_changes_cannot_select_other_observables(tmp_path, old, new):
    text = synthetic_babar_tex()
    assert old in text
    path, spec = primary_fixture(tmp_path, text.replace(old, new))
    with pytest.raises(ValueError, match="primary BaBar"):
        hist.babar_primary(path, spec)


def test_duplicate_momentum_header_is_ambiguous_even_with_unchanged_rows(tmp_path):
    text = synthetic_babar_tex().replace(r"$p_{t0}$ &", "$p_{t0}$ & duplicate header\n$p_{t0}$ &", 1)
    path, spec = primary_fixture(tmp_path, text)
    with pytest.raises(ValueError, match="ambiguous"):
        hist.babar_primary(path, spec)


@pytest.mark.parametrize("old,new", [("-.8 $\\pm$ .1", "-.8 $\\pm$ 0"),
                                     (".24 & & & .36", "1.24 & & & .36"),
                                     (".36 & .62", ".36 & 1.62")])
def test_primary_uncertainty_and_mean_domains_are_checked(tmp_path, old, new):
    text = synthetic_babar_tex()
    assert old in text
    path, spec = primary_fixture(tmp_path, text.replace(old, new))
    with pytest.raises(ValueError, match="primary BaBar"):
        hist.babar_primary(path, spec)


def test_primary_source_and_percentage_contract_are_required(tmp_path):
    archive, spec = primary_fixture(tmp_path)
    primary = hist.babar_primary(archive, spec)
    path, source = secondary_fixture(tmp_path)
    with pytest.raises(ValueError, match="missing or mismatched"):
        hist.parse_historical_xlsx(path, source, synthetic_contract(), {})
    source["adapter"]["input_unit"] = "fraction"
    with pytest.raises(ValueError, match="percentage"):
        hist.parse_historical_xlsx(path, source, synthetic_contract(), {}, babar_tables={spec["id"]: primary})


def test_primary_metadata_is_source_evidence_and_does_not_certify_acceptance(tmp_path):
    view, _, _, _ = candidate_view(tmp_path)
    candidate = {"schema_version": 1, "scope": hist.SCOPE, "contract_id": synthetic_contract()["id"],
                 "production_ready": False, "physical_extraction": False, "historical_fit_ready": False, "views": [view]}
    assert hist.validate_historical_candidate(candidate, synthetic_contract())["points"] == 2
    altered = copy.deepcopy(candidate)
    altered["views"][0]["points"][0]["angular_source"]["average_type"] = "ratio_of_angular_means"
    with pytest.raises(ValueError, match="unresolved acceptance"):
        hist.validate_historical_candidate(altered, synthetic_contract())
    altered = copy.deepcopy(candidate)
    altered["contract_id"] = "STALE-SYNTHETIC-CONTRACT"
    with pytest.raises(ValueError, match="contract identity"):
        hist.validate_historical_candidate(altered, synthetic_contract())
