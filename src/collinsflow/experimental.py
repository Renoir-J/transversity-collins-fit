"""Read pinned experimental sources into candidates, without constructing a fit likelihood.

Reported errors and overlapping projections remain separate. Unresolved covariance,
acceptance or theory conventions prevent admission to numerical fitting.
"""
from __future__ import annotations

import math
import re
import shutil
import tempfile
import tarfile
from itertools import product
from pathlib import Path, PurePosixPath

import requests
import yaml
from bs4 import BeautifulSoup
from openpyxl import load_workbook

from .provenance import atomic_json, digest, file_hash, read_json, source_fingerprint, utc_now

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def scalar(text):
    text = text.strip().replace("\N{MINUS SIGN}", "-")
    if not re.fullmatch(NUMBER, text):
        raise ValueError(f"ambiguous numeric cell: {text!r}")
    value = float(text)
    if not math.isfinite(value):
        raise ValueError("nonfinite experimental number")
    return value


def interval(text):
    match = re.fullmatch(rf"\[\s*({NUMBER})\s*,\s*({NUMBER})\s*\]", text)
    if not match:
        raise ValueError(f"ambiguous bin cell: {text!r}")
    low, high = map(float, match.groups())
    if low >= high:
        raise ValueError("nonpositive bin width")
    return [low, high]


def measurement(text, unit):
    pieces = re.split(r"\s*(?:\\pm|\N{PLUS-MINUS SIGN})\s*", text)
    if len(pieces) != 3:
        raise ValueError("expected central value, statistical error and systematic error")
    multiplier = {"fraction": 1.0, "percent": 0.01}[unit]
    value, stat, systematic = [scalar(part) * multiplier for part in pieces]
    if stat <= 0 or systematic < 0:
        raise ValueError("invalid uncertainty")
    return value, stat, systematic


def html_rows(path, table_id):
    soup = BeautifulSoup(Path(path).read_text(encoding="utf-8"), "html.parser")
    container = soup.find(id=table_id)
    if container is None or container.find("table") is None:
        raise ValueError(f"missing source table {table_id}")
    table = container.find("table")
    for node in table.find_all("math"):
        if not node.get("alttext"):
            raise ValueError("math cell lacks an unambiguous TeX annotation")
        node.replace_with(node["alttext"])
    rows = []
    for row in table.find_all("tr"):
        cells = row.find_all(["td", "th"], recursive=False)
        if any(int(c.get("colspan", 1)) != 1 or int(c.get("rowspan", 1)) != 1 for c in cells):
            raise ValueError("merged table cells need an explicit source-specific adapter")
        rows.append([c.get_text(" ", strip=True) for c in cells])
    return rows


def xlsx_records(path, expected_sheet=None):
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = workbook[expected_sheet] if expected_sheet else workbook.active
        rows = sheet.iter_rows(values_only=True)
        header = next(rows)
        names = [name.strip() for name in header if name is not None]
        if len(names) != len(set(names)):
            raise ValueError("duplicate column headers")
        return [{name.strip(): value for name, value in zip(header, row) if name is not None}
                for row in rows if any(value is not None for value in row)]
    finally:
        workbook.close()


def acquire(source, cache, seed=None, offline=False):
    """Acquire exactly the locked bytes. A content change requires source review."""
    cache = Path(cache)
    if Path(source["filename"]).name != source["filename"]:
        raise ValueError("source filename must be a basename")
    target = cache / source["filename"]
    cache.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if file_hash(target) != source["sha256"]:
            raise ValueError(f"cached source checksum mismatch: {target.name}")
        return target
    if seed and Path(seed).is_file():
        if file_hash(seed) != source["sha256"]:
            raise ValueError("seed source checksum mismatch")
        shutil.copyfile(seed, target)
        return target
    if offline:
        raise FileNotFoundError(f"locked source not cached: {source['id']}")
    if not source["url"].startswith("https://"):
        raise ValueError("experimental sources must use HTTPS")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=cache, delete=False) as stream:
            temporary = Path(stream.name)
            with requests.get(source["url"], timeout=45, stream=True) as response:
                response.raise_for_status()
                size = 0
                for chunk in response.iter_content(65536):
                    size += len(chunk)
                    if size > 20_000_000:
                        raise ValueError("source exceeds acquisition size limit")
                    stream.write(chunk)
        if file_hash(temporary) != source["sha256"]:
            raise ValueError(f"upstream bytes differ from locked source: {source['id']}")
        temporary.replace(target)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return target


def view(identifier, source, process, hadron, observable, overlap, convention,
         role="alternative", blockers=None, **metadata):
    return {"id": identifier, "source_id": source["id"], "publication": source["publication"],
            "source_url": source["url"], "source_sha256": source["sha256"], "version": source["version"],
            "process": process, "hadron": hadron, "observable": observable, "unit": "fraction",
            "overlap_group": overlap, "convention": convention, "role": role,
            "production_admission": {"status": "BLOCKED", "reasons": blockers or
                ["systematic_covariance_unresolved", "process_adapter_not_certified"]},
            "points": [], **metadata}


def point(identifier, value, stat, systematic, *, bins=None, averages=None, locator=None, **metadata):
    return {"id": identifier, "value": value, "statistical_error": stat,
            "systematic_components": [{"name": "reported_point_systematic", "sigma": systematic,
                                        "correlation": "unresolved"}],
            "bins": bins or {}, "averages": averages or {}, "source_locator": locator, **metadata}


def parse_hermes(path, source):
    config = source["adapter"]
    records = xlsx_records(path, config["sheet"])
    result = view(config["view_id"], source, "SIDIS", config["hadron"], "Collins_SFA",
                  "HERMES_2002_2005_" + config["hadron"], "sin(phi_h+phi_S)_epsilon_prefactor_removed",
                  role="primary_candidate" if config["dimension"] == 3 else "alternative_projection",
                  prediction_prescription="at_reported_mean_kinematics",
                  depolarization_removed=True, target="proton", beam_energy_GeV=27.6,
                  statistical_independence="64 separate kinematic bins for one moment/species" if
                      config["dimension"] == 3 else "do not combine overlapping 1D projections",
                  relative_systematics=[{"name": "target_polarization", "sigma": 0.073,
                      "correlation_group": "HERMES_2002_2005_target_polarization", "kind": "multiplicative"}],
                  replaces_event_sample="HERMES_2010_same_2002_2005_sample",
                  cuts={"Q2_GeV2_min": 1, "W2_GeV2_min": 10, "x": [0.023, 0.6],
                        "z": [0.2, 1.2 if config.get("projection") == "z" else 0.7], "pT_GeV_max": 2})
    for index, record in enumerate(records):
        record = dict(record)
        changes = []
        for correction in config.get("corrections", []):
            if correction["source_excel_row"] == index + 2:
                column = correction["column"]
                if record[column] != correction["old_value"]:
                    raise ValueError("reviewed correction no longer matches the locked source cell")
                record[column] = correction["new_value"]
                changes.append(correction)
        if (record.get("obs") != "AUTcollins" or record.get("hadron") != config["hadron"] or
                record.get("target") != "p" or record.get("col") != "HERMES" or
                record.get(config["depolarization_column"]) is not False):
            raise ValueError("HERMES table metadata disagrees with reviewed mapping")
        bins = {key: [float(record[low]), float(record[high])]
                for key, (low, high) in config["bins"].items()}
        result["points"].append(point(f"{result['id']}:{index}", float(record["value"]),
            float(record["stat_u"]), float(record["syst_u"]), bins=bins,
            averages={key: float(record[column]) for key, column in
                      {"Q2_GeV2": "Q2", "x": "x", "y": "y", "z": "z", "pT_GeV": "pT"}.items()},
            locator={"sheet": config["sheet"], "row": index + 2, "columns": ["value", "stat_u", "syst_u"]},
            source_corrections=changes,
            kinematic_flags=["high_z_semi_exclusive"] if record["z"] > 0.7 else []))
    if len(records) != config["expected_rows"]:
        raise ValueError("unexpected HERMES table row count")
    if config["dimension"] == 3:
        boxes = [p["bins"] for p in result["points"]]
        for i, left in enumerate(boxes):
            for right in boxes[i+1:]:
                if all(max(left[k][0], right[k][0]) < min(left[k][1], right[k][1]) for k in left):
                    raise ValueError("overlapping HERMES 3D bins")
    return [result]


def crosscheck_hermes_grid(views, primary_paper):
    rows = html_rows(primary_paper, "S3.T5")
    if len(rows) != 5 or any(len(row) != 3 for row in rows):
        raise ValueError("HERMES primary binning table schema changed")
    axes = []
    for column in range(3):
        axis = []
        for row in rows[1:]:
            numbers = re.findall(NUMBER, row[column])
            if len(numbers) != 2:
                raise ValueError("ambiguous primary-paper HERMES bin")
            axis.append(tuple(map(float, numbers)))
        axes.append(axis)
    expected = set(product(*axes))
    reports = []
    for dataset in views:
        if not dataset["id"].startswith("HERMES2020_") or not dataset["id"].endswith("_3D"):
            continue
        keys = ["x", "z", "pT_GeV"]
        actual = {tuple(tuple(item["bins"][key]) for key in keys) for item in dataset["points"]}
        if actual != expected or len(dataset["points"]) != 64:
            raise ValueError("HERMES 3D source does not reproduce all 64 primary-paper bins")
        for item in dataset["points"]:
            if any(not bounds[0] <= item["averages"][key] <= bounds[1] for key, bounds in item["bins"].items()):
                raise ValueError("HERMES mean lies outside its reviewed kinematic bin")
        reports.append({"view": dataset["id"], "points": 64, "primary_binning": "2007.07755v1 Table 5",
                        "corrected_source_cells": sum(len(p["source_corrections"]) for p in dataset["points"]),
                        "passed": True})
    return reports


def parse_belle(path, source):
    outputs = []
    for definition in source["adapter"]["tables"]:
        rows = html_rows(path, definition["table_id"])
        if len(rows) - 1 != definition["expected_rows"]:
            raise ValueError("unexpected Belle table row count")
        if "%" not in rows[0][-1] or not all(len(row) == definition["columns"] for row in rows):
            raise ValueError("Belle column/unit schema changed")
        boundary = definition["one_dimensional_rows"]
        for dimensionality, selected in [("1D", rows[1:boundary+1]), ("2D", rows[boundary+1:])]:
            for observable, column in definition["observables"].items():
                identifier = f"Belle2019_{definition['family']}_{dimensionality}_{observable}"
                result = view(identifier, source, "SIA", definition["hadron"], observable,
                    "Belle2019_980fb_pion_pair_sample", "thrust_axis_cos(phi1+phi2)_double_ratio",
                    role="primary_candidate" if definition["family"] == "charged_z" and
                         dimensionality == "2D" and observable == "A12_UL" else "alternative_projection",
                    prediction_prescription="bin_integrated_ratio_with_published_angular_acceptance",
                    sqrt_s_GeV=10.58, angular_factor_definition="mean_of_sin2_over_1_plus_cos2",
                    table_id=definition["table_id"],
                    cuts={"visible_energy_GeV_min": 7, "thrust_min": 0.8,
                          "thrust_polar_angle_rad": [1.34, 2.03], "meson_axis_opening_angle_rad_max": 0.3},
                    cross_publication_overlap="earlier Belle samples require replacement/overlap review")
                for index, row in enumerate(selected):
                    if definition["family"] == "charged_pt":
                        columns = {"pT1_GeV": (0, 1), "z1": (2, 3), "pT2_GeV": (4, 5), "z2": (6, 7)}
                    else:
                        columns = {"z1": (0, 1), "pT1_GeV": (2, 3), "z2": (4, 5), "pT2_GeV": (6, 7)}
                    bins = {key: interval(row[b]) for key, (b, _) in columns.items()}
                    means = {key: scalar(row[m]) for key, (_, m) in columns.items()}
                    means["angular_factor"] = scalar(row[8])
                    value, stat, systematic = measurement(row[column], "percent")
                    source_row = index + (1 if dimensionality == "1D" else boundary+1)
                    result["points"].append(point(f"{identifier}:{index}", value, stat, systematic,
                        bins=bins, averages=means, locator={"table": definition["table_id"],
                        "row_after_header": source_row, "column_index_1based": column+1},
                        kinematic_flags=["low_z_below_0.2"] if bins["z1"][0] < 0.2 or bins["z2"][0] < 0.2 else []))
                outputs.append(result)
    return outputs


def parse_besiii(path, source):
    rows = html_rows(path, "S0.T1")
    if len(rows) != 13 or any(len(row) != 7 for row in rows):
        raise ValueError("BESIII final-revision table schema changed")
    if "%" in rows[0][-1] or source["version"] != "1507.06824v2":
        raise ValueError("BESIII adapter requires the final fraction-valued v2 table")
    outputs = []
    for projection, data_rows in [("z", rows[1:7]), ("pT", rows[8:])]:
        for label, column in [("UL", 5), ("UC", 6)]:
            identifier = f"BESIII3650_{projection}_{label}"
            result = view(identifier, source, "SIA", "pi_charged_pair", f"A0_{label}",
                "BESIII3650_dipion_sample", "hadron_frame_cos(2phi0)_double_ratio",
                role="primary_candidate" if projection == "z" and label == "UL" else "alternative_projection",
                sqrt_s_GeV=3.65, angular_factor_definition="ratio_of_mean_sin2_to_mean_1_plus_cos2",
                prediction_prescription="bin_integrated_ratio_with_exchange_fold_and_angular_acceptance",
                low_scale_power_corrections="require sensitivity study")
            for index, row in enumerate(data_rows):
                if projection == "z":
                    matches = re.findall(r"\[[^\[\]]+\]", row[0])
                    if len(matches) != 2:
                        raise ValueError("BESIII exchanged z bins are ambiguous")
                    z1, z2 = map(interval, matches)
                    union = [{"z1": z1, "z2": z2}]
                    if z1 != z2:
                        union.append({"z1": z2, "z2": z1})
                    bins = {}
                    means = dict(zip(["z1", "z2", "pT_GeV", "angular_factor"], map(scalar, row[1:5])))
                else:
                    union, bins = [], {"pT_GeV": interval(row[0])}
                    means = dict(zip(["pT_GeV", "z1", "z2", "angular_factor"], map(scalar, row[1:5])))
                value, stat, systematic = measurement(row[column], "fraction")
                result["points"].append(point(f"{identifier}:{index}", value, stat, systematic,
                    bins=bins, averages=means, bin_union=union,
                    locator={"table": "S0.T1", "row_after_header": index+1 if projection == "z" else index+8,
                             "column_index_1based": column+1}))
            outputs.append(result)
    return outputs


def parse_star510(path, source):
    rows = html_rows(path, "as1_S5.T1")
    if len(rows) != 15 or any(len(row) != 6 for row in rows) or rows[0][1] != r"A_{\mathrm{UT}}":
        raise ValueError("STAR510 supplemental table schema changed")
    identifier = "STAR510_pi_plus_xF_positive_jet_pT_supplement"
    result = view(identifier, source, "STAR_510", "pi+", "AUT_Collins", "STAR510_pion_in_jet_sample",
        "published_STAR_sin(phiS-phiH)", role="partial_source_only",
        blockers=["jet_pT_bin_edges_not_in_supplement", "unresolved_systematic_correlations",
                  "different_published_binning_from_HEPData_figure2", "process_adapter_not_certified"],
        sqrt_s_GeV=510, prediction_prescription="bin_integrated_polarized_over_unpolarized_cross_sections",
        cuts={"jet_R": 0.5, "jet_axis": "published anti-kt jet axis", "eta": [0, 0.9],
              "z": [0.1, 0.8], "DeltaR_h_min": 0.02,
              "jT_max_GeV": "min((0.0093+0.2982*z)*jet_pT_GeV,2.5)"},
        polarization_scale_uncertainty="requires official table review")
    for index, row in enumerate(rows[1:]):
        mean, value, stat, trigger, cross, resolution = map(scalar, row)
        item = point(f"{identifier}:{index}", value, stat, 0, averages={"jet_pT_GeV": mean},
                     locator={"table": "as1_S5.T1", "row_after_header": index+1}, binning_status="mean_only")
        item["systematic_components"] = [{"name": name, "sigma": sigma, "correlation": "unresolved"}
            for name, sigma in [("trigger", trigger), ("cross", cross), ("resolution", resolution)]]
        result["points"].append(item)
    return [result]


def hepdata_documents(path):
    """Read safe YAML members in memory; TeX backslashes in filenames are not paths."""
    documents, size = {}, 0
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts or member.issym() or member.islnk():
                raise ValueError("unsafe HEPData archive member")
            if not member.isfile():
                continue
            size += member.size
            if member.size > 2_000_000 or size > 20_000_000:
                raise ValueError("HEPData archive exceeds decompression limits")
            if name.suffix not in {".yaml", ".yml"}:
                continue
            if name.name in documents:
                raise ValueError("duplicate HEPData table filename")
            with archive.extractfile(member) as stream:
                documents[name.name] = list(yaml.safe_load_all(stream))
    if "submission.yaml" not in documents:
        raise ValueError("missing HEPData submission metadata")
    return documents


def parse_star_hepdata(path, source):
    documents = hepdata_documents(path)
    config = source["adapter"]
    metadata = documents["submission.yaml"]
    if metadata[0]["hepdata_doi"] != config["doi"]:
        raise ValueError("HEPData submission DOI differs from source lock")
    by_file = {item["data_file"]: item for item in metadata[1:] if "data_file" in item}
    outputs = []
    for definition in config["tables"]:
        name = definition["data_file"]
        meta = by_file[name]
        table_docs = documents[name]
        if len(table_docs) != 1:
            raise ValueError("unexpected multiple YAML table documents")
        table = table_docs[0]
        if len(table["independent_variables"]) != 1 or len(table["dependent_variables"]) != 1:
            raise ValueError("STAR adapter requires one measured coordinate and asymmetry")
        independent, dependent = table["independent_variables"][0], table["dependent_variables"][0]
        header = dependent["header"]["name"].replace(r"\sin", "sin")
        if header != r"$A_{UT}^{sin(\phi_{S}-\phi_{H})}$":
            raise ValueError("HEPData table is not the pion Collins harmonic")
        if len(dependent["values"]) != definition["expected_rows"] or len(independent["values"]) != len(dependent["values"]):
            raise ValueError("STAR source row count changed")
        identifier = definition["view_id"]
        result = view(identifier, source, f"STAR_{config['energy_GeV']}", definition["hadron"], "AUT_Collins",
            f"STAR{config['energy_GeV']}_pion_in_jet_sample", "published_STAR_sin(phiS-phiH)",
            role="review_required",
            blockers=["physical_bin_edges_not_encoded_by_plot_systematic_intervals",
                      "systematic_and_projection_covariance_unresolved", "process_adapter_not_certified"],
            sqrt_s_GeV=config["energy_GeV"], table_id=meta["table_doi"],
            table_name=meta["name"], table_description=meta["description"],
            source_qualifiers=dependent.get("qualifiers", []),
            coordinate_header=independent["header"], sample_composition=definition["sample"],
            data_license=metadata[0].get("data_license"),
            prediction_prescription="requires physical bin/acceptance review; plot x errors are not integration bins")
        for index, (xvalue, yvalue) in enumerate(zip(independent["values"], dependent["values"])):
            errors = {}
            for error in yvalue["errors"]:
                if set(error) != {"label", "symerror"} or error["label"] in errors:
                    raise ValueError("unsupported or duplicate STAR error component")
                errors[error["label"]] = float(error["symerror"])
            if set(errors) != {"stat", "sys"}:
                raise ValueError("STAR stat/sys schema changed")
            averages, source_interval = {}, {}
            coordinate = definition["coordinate"]
            if "value" in xvalue:
                averages[coordinate] = float(xvalue["value"])
            if {"low", "high"} <= xvalue.keys():
                if xvalue["low"] >= xvalue["high"]:
                    raise ValueError("invalid STAR horizontal systematic interval")
                source_interval[coordinate] = {"low": float(xvalue["low"]), "high": float(xvalue["high"]),
                    "meaning": "reported_plot_systematic_interval_not_physical_binning"}
            if not averages and not source_interval:
                raise ValueError("missing STAR measured coordinate")
            result["points"].append(point(f"{identifier}:{index}", float(yvalue["value"]), errors["stat"],
                errors["sys"], averages=averages, source_intervals=source_interval,
                source_coordinate_record=xvalue, binning_status="physical_edges_require_review",
                locator={"archive_table": name, "row_after_header": index+1, "table_doi": meta["table_doi"]}))
        outputs.append(result)
    return outputs


def validate_candidate(candidate):
    if candidate.get("scope") != "experimental_candidate" or candidate.get("schema_version") != 1:
        raise ValueError("experimental candidate scope/schema required")
    identifiers, point_ids = set(), set()
    for dataset in candidate["views"]:
        if dataset["id"] in identifiers:
            raise ValueError("duplicate view ID")
        identifiers.add(dataset["id"])
        if dataset["production_admission"]["status"] != "BLOCKED" or not dataset["production_admission"]["reasons"]:
            raise ValueError("candidate sources cannot certify their own production admission")
        if not dataset["points"] or dataset["unit"] != "fraction":
            raise ValueError("missing data or noncanonical unit")
        for item in dataset["points"]:
            if item["id"] in point_ids:
                raise ValueError("duplicate experimental point ID")
            point_ids.add(item["id"])
            if not math.isfinite(item["value"]) or not math.isfinite(item["statistical_error"]) or item["statistical_error"] <= 0:
                raise ValueError("invalid central value/statistical uncertainty")
            for component in item["systematic_components"]:
                if component["sigma"] < 0 or not math.isfinite(component["sigma"]):
                    raise ValueError("invalid systematic uncertainty")
            for box in [item["bins"], *item.get("bin_union", [])]:
                for bounds in box.values():
                    if len(bounds) != 2 or not all(math.isfinite(v) for v in bounds) or bounds[0] >= bounds[1]:
                        raise ValueError("invalid experimental bin")
            for bounds in item.get("source_intervals", {}).values():
                if (bounds.get("meaning") != "reported_plot_systematic_interval_not_physical_binning"
                        or not all(math.isfinite(bounds[key]) for key in ("low", "high"))
                        or bounds["low"] >= bounds["high"]):
                    raise ValueError("invalid source plot-systematic interval")
            if not all(math.isfinite(v) for v in item["averages"].values()) or not item["source_locator"]:
                raise ValueError("invalid mean kinematics or source locator")
    return {"views": len(identifiers), "points_across_overlapping_views": len(point_ids),
            "production_ready": False}


def choose_views(candidate, identifiers):
    validate_candidate(candidate)
    by_id = {dataset["id"]: dataset for dataset in candidate["views"]}
    if len(identifiers) != len(set(identifiers)) or any(key not in by_id for key in identifiers):
        raise ValueError("unknown or duplicate selected view")
    chosen = [by_id[key] for key in identifiers]
    overlap = [dataset["overlap_group"] for dataset in chosen]
    if len(overlap) != len(set(overlap)):
        raise ValueError("selected views overlap; joint covariance or a disjoint selection is required")
    return chosen


def crosscheck_besiii(views, paths):
    reports = []
    index = {dataset["id"]: dataset for dataset in views}
    for number, key in [(3000, "BESIII3650_z_UL"), (3001, "BESIII3650_z_UC"),
                        (3002, "BESIII3650_pT_UL"), (3003, "BESIII3650_pT_UC")]:
        rows = xlsx_records(paths[f"JAM_SIA_{number}"])
        points = index[key]["points"]
        if len(rows) != len(points):
            raise ValueError("BESIII primary/JAM row-count mismatch")
        errors = [abs(item[field] - float(row[column])) for item, row in zip(points, rows)
                  for field, column in [("value", "value"), ("statistical_error", "stat_u")]]
        errors += [abs(item["systematic_components"][0]["sigma"] - float(row["sys_u"]))
                   for item, row in zip(points, rows)]
        errors += [abs(item["averages"][key] - float(row[column])) for item, row in zip(points, rows)
                   for key, column in [("z1", "z1"), ("z2", "z2"), ("pT_GeV", "pT"),
                                       ("angular_factor", "S2/1+C2")]]
        if any(row["Q"] != 3.65 for row in rows):
            raise ValueError("BESIII source energy mismatch")
        error = max(errors)
        if error > 1e-14:
            raise ValueError("BESIII final-paper values/errors disagree with pinned JAM transcription")
        reports.append({"view": key, "points": len(points), "max_abs_error": error, "passed": True})
    return reports


def _locked_sources(manifest):
    sources = manifest["sources"]
    if len({s["id"] for s in sources}) != len(sources) or len({s["filename"] for s in sources}) != len(sources):
        raise ValueError("duplicate locked source ID or filename")
    if any(not re.fullmatch(r"[a-f0-9]{64}", s["sha256"]) for s in sources):
        raise ValueError("invalid locked source SHA256")
    return sources


def _candidate_from_sources(manifest_sha256, sources, paths):
    views = []
    adapters = {"hermes_xlsx": parse_hermes, "belle_html": parse_belle,
                "besiii_html": parse_besiii, "star510_html": parse_star510,
                "star_hepdata": parse_star_hepdata}
    for source in sources:
        kind = source.get("adapter", {}).get("kind")
        if kind:
            views.extend(adapters[kind](paths[source["id"]], source))
    return {"schema_version": 1, "scope": "experimental_candidate", "production_ready": False,
            "manifest_sha256": manifest_sha256, "views": views,
            "sources": [{"id": s["id"], "url": s["url"], "sha256": s["sha256"],
                         "version": s["version"]} for s in sources]}


def _source_decisions(candidate, paths):
    views = candidate["views"]
    summary = validate_candidate(candidate)
    summary["primary_candidate_views"] = [v["id"] for v in views if v["role"] == "primary_candidate"]
    choose_views(candidate, summary["primary_candidate_views"])
    summary["BESIII_primary_vs_JAM"] = crosscheck_besiii(views, paths)
    summary["HERMES_3D_primary_binning"] = crosscheck_hermes_grid(views, paths["arxiv_2007.07755v1"])
    summary["STAR_source_review"] = {"physical_bin_edges_inferred_from_plot_errors": False,
        "combined_and_component_year_samples_kept_as_alternatives": True,
        "hepdata_views": sum("source_intervals" in v["points"][0] for v in views)}
    lattice = xlsx_records(paths["JAM_lattice_2000"])
    if len(lattice) != 1 or lattice[0]["value"] != 0.936 or lattice[0]["stat_u"] != 0.025 or lattice[0]["Q"] != 2:
        raise ValueError("lattice source differs from the selected original-paper constraint")
    summary["lattice_JAM_crosscheck"] = {"passed": True, "value": 0.936, "sigma": 0.025, "mu_GeV": 2}
    summary["source_manifest_digest"] = digest(candidate["sources"])
    summary["source_manifest_sha256"] = candidate["manifest_sha256"]
    summary["status"] = "CANDIDATES_VALIDATED_PRODUCTION_BLOCKED"
    return summary


def _check_source_snapshot(manifest_path, manifest_sha256, fingerprint, sources, cache):
    if file_hash(manifest_path) != manifest_sha256:
        raise ValueError("experimental source manifest changed during intake verification")
    if source_fingerprint(Path(__file__).resolve().parents[2]) != fingerprint:
        raise ValueError("experimental parser/configuration changed during intake verification")
    for source in sources:
        acquire(source, cache, offline=True)


def run_intake(manifest_path, cache, output, *, seed_root=None, offline=False):
    manifest_sha256 = file_hash(manifest_path)
    manifest = read_json(manifest_path)
    sources = _locked_sources(manifest)
    fingerprint = source_fingerprint(Path(__file__).resolve().parents[2])
    paths = {}
    for source in sources:
        seed = Path(seed_root) / source["seed_path"] if seed_root else None
        paths[source["id"]] = acquire(source, cache, seed, offline)
    candidate = _candidate_from_sources(manifest_sha256, sources, paths)
    summary = _source_decisions(candidate, paths)
    _check_source_snapshot(manifest_path, manifest_sha256, fingerprint, sources, cache)
    summary["source_fingerprint"] = fingerprint
    summary["completed_at"] = utc_now()
    output = Path(output)
    atomic_json(output / "candidates.json", candidate)
    summary["candidate_file_sha256"] = file_hash(output / "candidates.json")
    atomic_json(output / "intake-report.json", summary)
    return summary


def verify_intake(manifest_path, cache, output):
    manifest_sha256 = file_hash(manifest_path)
    manifest = read_json(manifest_path)
    sources = _locked_sources(manifest)
    fingerprint = source_fingerprint(Path(__file__).resolve().parents[2])
    output = Path(output)
    report_path, candidate_path = output / "intake-report.json", output / "candidates.json"
    report_sha256, candidate_sha256 = file_hash(report_path), file_hash(candidate_path)
    report = read_json(report_path)
    if report.get("source_manifest_sha256") != manifest_sha256:
        raise ValueError("stale experimental source manifest")
    if report.get("source_fingerprint") != fingerprint:
        raise ValueError("stale experimental parser/configuration")
    if report.get("candidate_file_sha256") != candidate_sha256:
        raise ValueError("experimental candidate artifact was changed")
    paths = {source["id"]: acquire(source, cache, offline=True) for source in sources}
    candidate = read_json(candidate_path)
    validate_candidate(candidate)
    reconstructed = _candidate_from_sources(manifest_sha256, sources, paths)
    if digest(candidate) != digest(reconstructed):
        raise ValueError("experimental candidate differs from locked source replay")
    expected = _source_decisions(reconstructed, paths)
    expected.update(source_fingerprint=fingerprint, candidate_file_sha256=candidate_sha256)
    if (set(report) != set(expected) | {"completed_at"}
            or digest({key: value for key, value in report.items() if key != "completed_at"}) != digest(expected)
            or not isinstance(report["completed_at"], str) or not report["completed_at"]):
        raise ValueError("experimental intake report differs from locked source replay")
    _check_source_snapshot(manifest_path, manifest_sha256, fingerprint, sources, cache)
    if file_hash(candidate_path) != candidate_sha256 or file_hash(report_path) != report_sha256:
        raise ValueError("experimental cached artifacts changed during source replay")
    return {"verified": True, "production_ready": False, "views": expected["views"],
            "points_across_overlapping_views": expected["points_across_overlapping_views"]}
