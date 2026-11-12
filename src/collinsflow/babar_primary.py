"""Read reviewed BaBar RF0 table fields as data, without fitting or executing TeX."""
from __future__ import annotations

from decimal import Decimal
import math
import re

from .experimental import NUMBER


def compact(text):
    return re.sub(r"\s+", "", text)


def decimal_cell(cell):
    value = cell.replace("$", "").strip()
    if not re.fullmatch(NUMBER, value):
        raise ValueError("ambiguous primary BaBar numeric cell")
    result = Decimal(value)
    if not result.is_finite():
        raise ValueError("nonfinite primary BaBar number")
    return result


def bin_cell(cell):
    match = re.fullmatch(rf"\[\s*({NUMBER})\s*,\s*({NUMBER})\s*\]", cell.replace("$", "").strip())
    if not match:
        raise ValueError("ambiguous primary BaBar finite bin")
    low, high = map(Decimal, match.groups())
    if not 0 <= low < high:
        raise ValueError("invalid primary BaBar bin")
    return low, high


def amplitude_cell(cell):
    pieces = cell.replace("$", "").split(r"\pm")
    if len(pieces) != 3:
        raise ValueError("primary BaBar amplitude needs stat and systematic errors")
    value, stat, systematic = map(decimal_cell, pieces)
    if stat <= 0 or systematic < 0:
        raise ValueError("invalid primary BaBar uncertainty")
    return {"value": value, "statistical": stat, "reported_systematic": systematic}


def data_rows(lines, columns):
    result = []
    for line_number, line in lines:
        if not re.match(r"\s*\$\[", line):
            continue
        line = line.strip()
        if not line.endswith(r"\\"):
            raise ValueError("incomplete primary BaBar row")
        cells = line[:-2].split("&")
        if len(cells) != columns:
            raise ValueError("primary BaBar row column count changed")
        result.append((line_number, cells))
    return result


def checked_fields(fields, bins, line, table):
    if (not all(0 < fields[k] <= 1 for k in ("z1", "z2")) or fields["pT"] < 0 or
            not 0 <= fields["S2/1+C2"] <= 1):
        raise ValueError("invalid primary BaBar mean kinematics")
    for key, (low, high) in bins.items():
        if not low <= fields[key] <= high:
            raise ValueError("primary BaBar mean outside its bin")
        fields[key+"min"], fields[key+"max"] = low, high
    return {"fields": fields, "bins": {k: list(map(float, v)) for k, v in bins.items()},
            "source_line": line, "table": table}


def parse_tables(text, tables, spec):
    """Select RF0 using source identity and column meaning, not numeric agreement."""
    labels = {"tab:totZul", "tab:totZuc", "tab:totPt"}
    if (set(tables) != labels or set(spec["table_rows"]) != labels or
            spec["measurement_unit"] != "percent" or spec["frame"] != "hadron2"):
        raise ValueError("unreviewed primary BaBar table/frame/unit contract")
    clean = re.sub(r"(?<!\\)%[^\n]*", "", text)
    for macro, ratio in (("aul", "UL"), ("auc", "UC")):
        definitions = [compact(line) for line in clean.splitlines()
                       if re.match(r"\s*\\def\\"+macro+r"\s*\{", line)]
        expected = "\\def\\"+macro+"{\\ensuremath{{A_{0}^{"+ratio+"}}}}"
        if definitions != [expected]:
            raise ValueError("primary BaBar RF0 amplitude macro changed")
    result = {}
    for ratio, label in (("UL", "tab:totZul"), ("UC", "tab:totZuc")):
        lines = tables[label]
        headers = [compact(line).removesuffix(r"\\").split("&") for _, line in lines
                   if compact(line).startswith(r"\multirow{2}{*}{$z_1$}&")]
        macro = "aul" if ratio == "UL" else "auc"
        if (len(headers) != 1 or len(headers[0]) != 12 or
                r"\langlez_1\rangle" not in headers[0][1] or
                r"\langlez_2\rangle" not in headers[0][4] or
                headers[0][6] != r"$\langlep_{t0}\rangle$" or
                r"\left\langle\frac{s_{2}}{1+c_{2}}\right\rangle" not in headers[0][8] or
                headers[0][9] != "\\multirow{2}{*}{\\"+macro+r"Th(\%)}" or
                headers[0][10] or headers[0][11] != "\\multirow{2}{*}{\\"+macro+r"(\%)}" or
                "percent" not in " ".join(line for _, line in lines)):
            raise ValueError("primary BaBar z-table RF0/header/unit identity changed")
        rows = data_rows(lines, 12)
        if len(rows) != spec["table_rows"][label]:
            raise ValueError("incomplete primary BaBar z table")
        points = []
        seen = set()
        for line, cells in rows:
            bins = {"z1": bin_cell(cells[0]), "z2": bin_cell(cells[3])}
            key = tuple(bins.items())
            if key in seen:
                raise ValueError("duplicate primary BaBar z bin")
            seen.add(key)
            fields = {k: decimal_cell(cells[i]) for k, i in (("z1", 1), ("z2", 4), ("pT", 6), ("S2/1+C2", 8))}
            if cells[10].strip():
                raise ValueError("primary BaBar z-table separator changed")
            points.append(checked_fields({**fields, **amplitude_cell(cells[11])}, bins, line, label))
        result[(ratio, "z1_z2")] = {"frame": "hadron2", "ratio": ratio, "projection": "z1_z2", "points": points}
    lines = tables["tab:totPt"]
    starts = [i for i, (_, line) in enumerate(lines) if re.match(r"\s*\$p_\{t0\}\$\s*&", line)]
    if len(starts) != 1:
        raise ValueError("missing or ambiguous primary BaBar RF0 momentum header")
    section = lines[starts[0]:]
    first_row = next((i for i, (_, line) in enumerate(section) if re.match(r"\s*\$\[", line)), None)
    if first_row is None:
        raise ValueError("missing primary BaBar RF0 momentum rows")
    header = compact("".join(line for _, line in section[:first_row])).split(r"\\", 1)[0].split("&")
    if (len(header) != 10 or header[0] != r"$p_{t0}$" or header[1] != r"$\langlep_{t0}\rangle$" or
            r"\langlez_1\rangle" not in header[2] or any(header[i] for i in (3, 4, 8)) or
            r"\langlez_2\rangle" not in header[5] or
            r"\left\langle\frac{\sin^2\theta_{2}}{1+\cos^2\theta_{2}}\right\rangle" not in header[6] or
            header[7] != r"\multirow{2}{*}{\aul}" or header[9] != r"\multirow{2}{*}{\auc}"):
        raise ValueError("primary BaBar momentum-table RF0/header identity changed")
    rows = data_rows(section, 10)
    if len(rows) != spec["table_rows"]["tab:totPt"]:
        raise ValueError("incomplete primary BaBar momentum table")
    for ratio, index in (("UL", 7), ("UC", 9)):
        points = []
        seen = set()
        for line, cells in rows:
            bounds = bin_cell(cells[0])
            if bounds in seen or any(cells[i].strip() for i in (3, 4, 8)):
                raise ValueError("duplicate bin or changed primary BaBar momentum separators")
            seen.add(bounds)
            fields = {k: decimal_cell(cells[i]) for k, i in (("z1", 2), ("z2", 5), ("pT", 1), ("S2/1+C2", 6))}
            points.append(checked_fields({**fields, **amplitude_cell(cells[index])}, {"pT": bounds}, line, "tab:totPt"))
        result[(ratio, "pT")] = {"frame": "hadron2", "ratio": ratio, "projection": "pT", "points": points}
    return {"id": spec["id"], "publication": spec["publication"], "sha256": spec["sha256"],
            "member": spec["member"], "member_sha256": spec["member_sha256"], "views": result}


def compare_row(record, reference, adapter):
    """Require exact published decimals before leaving existing normalization intact."""
    if adapter["input_unit"] != "percent":
        raise ValueError("primary BaBar percentage unit required")
    columns = {"statistical": adapter["error_columns"]["statistical"],
               "reported_systematic": adapter["error_columns"]["reported_systematic"],
               "S2/1+C2": adapter["angular_column"]}
    for key, expected in reference["fields"].items():
        value = record.get(columns.get(key, key))
        if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or
                Decimal(str(value)) != expected):
            raise ValueError("secondary BaBar row disagrees with primary numerical table: "+key)
    return len(reference["fields"])
