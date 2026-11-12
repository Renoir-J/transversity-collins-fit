"""Pinned historical candidates with explicit source repairs; no fit admission."""
from __future__ import annotations

import hashlib
import math
import re
import tarfile
from pathlib import Path, PurePosixPath

from .experimental import NUMBER, acquire, scalar, xlsx_records
from .provenance import atomic_json, file_hash, read_json, source_fingerprint, utc_now

SCOPE = "historical_benchmark_candidate"


def number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"finite numeric {name} required")
    if (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(f"invalid {name}")
    return float(value)


def tex_table(text, label):
    # Comments include an obsolete duplicate table in the pinned source.
    # Reading TeX as text does not execute its macros or file commands.
    clean = re.sub(r"(?<!\\)%[^\n]*", "", text)
    matches = [m for m in re.finditer(r"\\begin\{table\*?\}.*?\\end\{table\*?\}", clean, re.S)
               if "\\label{" + label + "}" in m.group()]
    if len(matches) != 1:
        raise ValueError(f"missing or ambiguous primary table {label}")
    match = matches[0]
    start = clean[:match.start()].count("\n") + 1
    return [(start + i, line) for i, line in enumerate(match.group().splitlines())]


def tex_number(cell):
    return scalar(cell.replace("$", "").replace("~", "").strip())


def tex_measurement(cell):
    cell = cell.replace("$", "").replace("~", " ").strip()
    pattern = rf"({NUMBER})\s*\\pm\s*({NUMBER})\s*\\pm\s*({NUMBER})\s*\(\s*({NUMBER})\s*\)"
    match = re.fullmatch(pattern, cell)
    if not match:
        raise ValueError("ambiguous primary Collins measurement")
    value, stat, experimental, model = map(scalar, match.groups())
    if stat <= 0 or min(experimental, model) < 0:
        raise ValueError("invalid primary uncertainty")
    return {"central_value": value, "statistical_error": stat,
            "experimental_systematic": experimental, "model_systematic": model}


def parse_halla_tex(text):
    """Read Tables I/II of the pinned final source, keeping neutron/3He distinct."""
    kinematics = {}
    for line_number, line in tex_table(text, "table:kine"):
        cells = line.split("&")
        if len(cells) != 11 or not re.match(r"\s*\d", cells[0]):
            continue
        values = [tex_number(v) for v in cells[:7]]
        x, q2, y, z, momentum, w, wp = values
        if not 0 < x < 1 or x in kinematics:
            raise ValueError("invalid/duplicate primary x")
        kinematics[x] = dict(zip(("x", "Q2_GeV2", "y", "z", "pT_GeV", "W_GeV", "Wprime_GeV"), values))
        kinematics[x]["source_line"] = line_number
    table = tex_table(text, "table:results")
    header = next((line.split("&") for _, line in table if "Collins Moment" in line), None)
    if header is None or len(header) != 6 or any("Collins Moment" not in header[i] for i in (2, 3)) or any(
            "Sivers Moment" not in header[i] for i in (4, 5)):
        raise ValueError("primary harmonic columns changed")
    result = {}
    for line_number, line in table:
        cells = line.split("&")
        if len(cells) != 6 or cells[0].strip() != "n":
            continue
        x = tex_number(cells[1])
        if x not in kinematics:
            raise ValueError("primary neutron/kinematic rows disagree")
        for charge, column in (("pi+", 2), ("pi-", 3)):
            key = (charge, x)
            if key in result:
                raise ValueError("duplicate primary neutron row")
            result[key] = {**kinematics[x], **tex_measurement(cells[column]),
                           "measurement_source_line": line_number}
    if len(kinematics) != 4 or len(result) != 8:
        raise ValueError("incomplete primary Hall A tables")
    return result


def halla_primary(path, spec):
    if file_hash(path) != spec["sha256"]:
        raise ValueError("changed primary archive")
    with tarfile.open(path, mode="r:*") as archive:
        members = archive.getmembers()
        if len(members) > 256 or sum(m.size for m in members) > 32_000_000:
            raise ValueError("oversized primary archive")
        for member in members:
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts or member.issym() or member.islnk():
                raise ValueError("unsafe primary archive member")
        chosen = [m for m in members if m.name == spec["member"] and m.isfile()]
        if len(chosen) != 1 or chosen[0].size > 2_000_000:
            raise ValueError("missing/oversized primary TeX member")
        data = archive.extractfile(chosen[0]).read()
    if hashlib.sha256(data).hexdigest() != spec["member_sha256"]:
        raise ValueError("changed primary TeX source")
    return parse_halla_tex(data.decode("utf-8"))


def row_bounds(row, averages):
    result = {}
    for name, key in (("x", "x"), ("z", "z"), ("pT", "pT_GeV"), ("z1", "z1"), ("z2", "z2")):
        low, high = row.get(name + "min"), row.get(name + "max")
        if low is None and high is None:
            continue
        a, b = number(low, name + "min"), number(high, name + "max")
        if a >= b or (key in averages and not a <= averages[key] <= b):
            raise ValueError("invalid published bounds or mean outside bin")
        result[key] = [a, b]
    return result


def sidis_point(row, adapter, contract, primary):
    corrected = adapter["uncertainty_layout"] == "stat_exp_model_remap"
    remap = contract["HallA_remap"]["fields"] if corrected else {
        "Q2_GeV2": "Q2", "x": "x", "y": "y", "z": "z", "pT_GeV": "pT", "central_value": "value"}
    values = {key: number(row[field], key) for key, field in remap.items()}
    if corrected:
        row_id = row[contract["HallA_remap"]["row_id_field"]]
        if isinstance(row_id, bool) or row_id not in (1, 2, 3, 4):
            raise ValueError("Hall A remap requires the pinned legacy row-id layout")
        key = (adapter["hadron"], values["x"])
        reference = primary.get(key)
        if reference is None or any(values[k] != reference[k] for k in remap):
            raise ValueError("Hall A remap disagrees with primary numerical tables")
        stat = values["statistical_error"]
        systematic = {"experimental": values["experimental_systematic"], "model": values["model_systematic"]}
    elif adapter["uncertainty_layout"] == "stat_exp":
        stat = number(row["stat_u"], "statistical error", positive=True)
        systematic = {"experimental": number(row["syst_u"], "experimental systematic", nonnegative=True)}
    else:
        stat, systematic = None, {}
    averages = {key: values[key] for key in ("x", "y", "z", "Q2_GeV2", "pT_GeV")}
    if not all(0 < averages[k] < 1 for k in ("x", "y", "z")) or averages["Q2_GeV2"] <= 0 or averages["pT_GeV"] < 0:
        raise ValueError("invalid SIDIS mean kinematics")
    if stat is not None and (stat <= 0 or any(v < 0 for v in systematic.values())):
        raise ValueError("invalid remapped statistical/systematic uncertainty")
    result = {"central_value": values["central_value"], "statistical_error": stat,
              "systematic_components": systematic, "averages": averages, "bins": row_bounds(row, averages)}
    if stat is None:
        result["reported_error_u"] = number(row["error_u"], "reported error_u", positive=True)
        result["uncertainty_status"] = "reported_error_u_composition_unresolved"
    if corrected:
        result["source_remap"] = {"kind": "reviewed_legacy_column_remap", "fields": remap, "source_bin_id": row_id,
            "primary_tables": ["table:kine", "table:results"], "primary_lines": [reference["source_line"], reference["measurement_source_line"]]}
    return result


def sia_point(row, adapter):
    if adapter["input_unit"] != "percent":
        raise ValueError("historical SIA input must use its declared percent unit")
    averages = {"z1": number(row["z1"], "z1"), "z2": number(row["z2"], "z2"),
                "Q2_GeV2": number(row["Q"], "Q", positive=True) ** 2}
    if not 0 < averages["z1"] <= 1 or not 0 < averages["z2"] <= 1:
        raise ValueError("invalid SIA fractions")
    for source_key, key in (("pT", "pT_GeV"), ("qT", "qT_GeV")):
        if row.get(source_key) is not None:
            averages[key] = number(row[source_key], source_key, nonnegative=True)
    columns = adapter["error_columns"]
    if any(name not in row for name in columns.values()):
        raise ValueError("declared historical SIA error column missing")
    angular_column = adapter["angular_column"]
    angular = number(row[angular_column], "angular source value", nonnegative=True)
    if angular > 1:
        raise ValueError("invalid angular source value")
    return {"central_value": .01 * number(row["value"], "asymmetry"),
            "statistical_error": .01 * number(row[columns["statistical"]], "statistical error", positive=True),
            "systematic_components": {"reported_systematic": .01 * number(row[columns["reported_systematic"]], "reported systematic error", nonnegative=True)},
            "uncertainty_source_columns": columns,
            "averages": averages, "bins": row_bounds(row, averages),
            "angular_source": {"column": angular_column, "value": angular, "average_type": "publication_review_pending"}}


def parse_historical_xlsx(path, source, contract, primary):
    if file_hash(path) != source["sha256"]:
        raise ValueError("historical source checksum mismatch")
    config = source["adapter"]
    records = xlsx_records(path, config["sheet"])
    if len(records) != config["expected_raw_rows"] or list(records[0]) != config["expected_headers"]:
        raise ValueError("historical source row/header schema changed")
    sidis = config["kind"] == "sidis"
    if config["kind"] not in ("sidis", "sia") or (sidis and config["input_unit"] != "fraction"):
        raise ValueError("unknown historical adapter or input unit")
    blockers = ["historical candidate only; no approved likelihood or production admission",
                "historical projection overlap/covariance prescription is not approved for updated primary fits"]
    if config.get("uncertainty_layout") != "stat_exp_model_remap":
        blockers.append("secondary transcription primary numerical crosscheck pending")
    if config.get("uncertainty_layout") == "reported_unspecified":
        blockers.append("source error_u statistical/systematic composition unresolved")
    if config.get("depolarization_status") == "publication_review_pending":
        blockers.append("publication depolarization convention review pending")
    if not sidis:
        blockers.extend(["publication angular-average type review pending", "historical momentum integration/acceptance adapter pending"])
    if config.get("uncertainty_layout") == "stat_exp_model_remap":
        blockers.append("original historical fit inclusion of Hall A model-systematic errors unresolved")
    dataset = {"id": source["id"], "source_id": source["id"], "publication": source["publication"],
        "source_sha256": source["sha256"], "source_version": source["version"], "source_url": source["url"],
        "process": "SIDIS" if sidis else "SIA", "projection": config["projection"],
        "input_units": config["input_unit"], "output_units": "fraction", "source_observable": config["source_obs"],
        "overlap_group": config["experiment"] + "_" + source["publication"],
        "production_ready": False, "historical_fit_ready": False, "blockers": blockers, "points": [], "excluded": []}
    if sidis:
        dataset.update(target=config["target"], hadron=config["hadron"],
                       depolarization_removed=False if config["depolarization_status"] == "raw_DNN_included_primary_verified" else None,
                       prediction_prescription="reported mean kinematics")
        if config["target"] == "neutron":
            dataset.update(target_material="3He", target_result="published extracted neutron; no second nuclear extraction")
        if config["experiment"] == "hermes":
            dataset["replaced_by_in_updated_fit"] = "HERMES2020_same_2002_2005_event_sample"
    else:
        dataset.update(ratio=config["ratio"], frame="hadron2",
                       prediction_prescription="reported P mean" if config["projection"] == "pT" else "P-integrated source view; theory adapter pending")
    for i, record in enumerate(records, start=2):
        if record["obs"] != config["source_obs"] or str(record["col"]).lower() != config["experiment"]:
            raise ValueError("historical experiment/harmonic identity changed")
        if sidis and (record["target"] != config["target"] or record["hadron"] != config["hadron"] or
                      str(record["dependence"]).lower() != config["projection"]):
            raise ValueError("historical target/species/projection mismatch")
        if not sidis and (record["hadron1"] != "pi" or record["hadron2"] != "pi"):
            raise ValueError("historical SIA is pion-pair only")
        point = sidis_point(record, config, contract, primary) if sidis else sia_point(record, config)
        point.update(id=f"{source['id']}:r{i}", locator={"sheet": config["sheet"], "excel_row": i})
        if sidis and config["projection"] == "pt" and point["averages"]["pT_GeV"] >= contract["SIDIS_projection_cut"]["upper_GeV"]:
            dataset["excluded"].append({"locator": point["locator"], "reported_pT_GeV": point["averages"]["pT_GeV"],
                "reason": "fails published strict pT<0.8 GeV projection cut"})
        else:
            dataset["points"].append(point)
    if len(dataset["points"]) != config["expected_points"]:
        raise ValueError("historical point selection disagrees with reviewed count")
    return dataset


def validate_historical_candidate(candidate, contract):
    if candidate.get("schema_version") != 1 or candidate.get("scope") != SCOPE or any(
            candidate.get(k) is not False for k in ("production_ready", "physical_extraction", "historical_fit_ready")):
        raise ValueError("historical candidates cannot self-certify fit admission")
    views = candidate["views"]
    if len({v["id"] for v in views}) != len(views):
        raise ValueError("duplicate historical view")
    ids = set()
    for view in views:
        if view["production_ready"] is not False or view["historical_fit_ready"] is not False or not view["blockers"]:
            raise ValueError("historical view cannot self-certify")
        for p in view["points"]:
            if p["id"] in ids:
                raise ValueError("duplicate historical point")
            ids.add(p["id"])
            number(p["central_value"], "central value")
            if p["statistical_error"] is None:
                if p.get("uncertainty_status") != "reported_error_u_composition_unresolved":
                    raise ValueError("missing uncertainty provenance")
                number(p["reported_error_u"], "reported error", positive=True)
            else:
                number(p["statistical_error"], "statistical error", positive=True)
            for value in p["systematic_components"].values():
                number(value, "systematic error", nonnegative=True)
    counts = {"views": len(views), "points": len(ids)}
    for process in ("SIDIS", "SIA"):
        chosen = [v for v in views if v["process"] == process]
        counts[process + "_views"] = len(chosen)
        counts[process + "_points"] = sum(len(v["points"]) for v in chosen)
    if counts != contract["expected"]:
        raise ValueError("historical totals do not match the source contract")
    return counts


def run_historical_intake(manifest_path, cache, output, *, seed_root=None, offline=False):
    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent.parent
    manifest = read_json(manifest_path)
    contract_path = root / manifest["contract"]
    contract = read_json(contract_path)
    if manifest.get("scope") != SCOPE or manifest.get("production_ready") is not False or manifest.get("physical_extraction") is not False:
        raise ValueError("invalid historical manifest scope")
    fingerprint = source_fingerprint(root)
    hashes = {str(p): file_hash(p) for p in (manifest_path, contract_path)}
    paths = {}
    for source in [manifest["primary_crosscheck"], *manifest["sources"]]:
        seed = Path(seed_root) / source["seed_path"] if seed_root else None
        if seed is not None and not seed.is_file():
            seed = Path(seed_root) / "data/raw/historical" / source["filename"]
        path = acquire(source, cache, seed=seed if seed is not None and seed.is_file() else None, offline=offline)
        paths[source["id"]] = path
        hashes[str(path.resolve())] = file_hash(path)
    primary = halla_primary(paths[manifest["primary_crosscheck"]["id"]], manifest["primary_crosscheck"])
    views = [parse_historical_xlsx(paths[s["id"]], s, contract, primary) for s in manifest["sources"]]
    candidate = {"schema_version": 1, "scope": SCOPE, "production_ready": False, "physical_extraction": False,
        "historical_fit_ready": False, "contract_id": contract["id"], "views": views,
        "manifest_sha256": file_hash(manifest_path), "contract_sha256": file_hash(contract_path), "source_fingerprint": fingerprint}
    counts = validate_historical_candidate(candidate, contract)
    if source_fingerprint(root) != fingerprint or any(file_hash(path) != value for path, value in hashes.items()):
        raise ValueError("historical source/parser changed during intake")
    output = Path(output)
    atomic_json(output / "candidates.json", candidate)
    report = {"status": "CANDIDATES_VALIDATED", "scope": SCOPE, "production_ready": False, "physical_extraction": False,
        "historical_fit_ready": False, "counts": counts, "source_hashes": hashes, "source_fingerprint": fingerprint,
        "manifest": str(manifest_path), "cache_root": str(Path(cache).resolve()),
        "candidate_sha256": file_hash(output / "candidates.json"), "primary_HallA_points_crosschecked": len(primary),
        "excluded_source_rows": sum(len(v["excluded"]) for v in views), "limitations": contract["limitations"], "completed_at": utc_now()}
    atomic_json(output / "report.json", report)
    return report


def verify_historical_intake(manifest_path, cache, output):
    """Revalidate cached outputs, locked inputs and current parser/configuration."""
    manifest_path, cache, output = Path(manifest_path).resolve(), Path(cache).resolve(), Path(output)
    root = manifest_path.parent.parent
    manifest = read_json(manifest_path)
    contract_path = root / manifest["contract"]
    contract = read_json(contract_path)
    report, candidate = read_json(output / "report.json"), read_json(output / "candidates.json")
    if report.get("status") != "CANDIDATES_VALIDATED" or report.get("scope") != SCOPE or any(
            report.get(k) is not False for k in ("production_ready", "physical_extraction", "historical_fit_ready")):
        raise ValueError("cached historical report cannot admit fitting")
    current = source_fingerprint(root)
    if candidate["source_fingerprint"] != current or report["source_fingerprint"] != current:
        raise ValueError("stale historical parser/configuration")
    if report["candidate_sha256"] != file_hash(output / "candidates.json"):
        raise ValueError("changed historical candidates")
    if candidate["manifest_sha256"] != file_hash(manifest_path) or candidate["contract_sha256"] != file_hash(contract_path):
        raise ValueError("changed historical source specification")
    sources = [manifest["primary_crosscheck"], *manifest["sources"]]
    expected = {str(manifest_path): file_hash(manifest_path), str(contract_path): file_hash(contract_path),
                **{str(cache / source["filename"]): source["sha256"] for source in sources}}
    if report["manifest"] != str(manifest_path) or report["cache_root"] != str(cache) or report["source_hashes"] != expected:
        raise ValueError("missing or changed historical input bindings")
    if any(file_hash(path) != value for path, value in expected.items()):
        raise ValueError("changed locked historical input")
    counts = validate_historical_candidate(candidate, contract)
    if report["counts"] != counts or report["primary_HallA_points_crosschecked"] != 8 or report["excluded_source_rows"] != 14:
        raise ValueError("changed historical intake decisions")
    return {"scope": SCOPE, "production_ready": False, "physical_extraction": False,
            "historical_fit_ready": False, "counts": counts}
