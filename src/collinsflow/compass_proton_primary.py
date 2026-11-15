"""Partial primary-text checks; Sigma and unequal legacy means stay unresolved."""
from __future__ import annotations

from decimal import Decimal
import hashlib
import math
from pathlib import PurePosixPath
import re
import tarfile

from .experimental import NUMBER
from .provenance import file_hash

MEANS = ("Q2", "y", "pT", "x", "z", "W")
COORDINATES = {"x": "x", "pt": "pT", "z": "z"}
QUANTITIES = {"x": "x", "pt": "pT_GeV", "z": "z"}
HEADER_SUFFIX = "Range Asym +/- Sigma | <Q^2> (GeV/c)^2 <y> <p_T> (GeV/c) <x> <z> <W> (GeV/c^2) |"
HEADERS = {prefix+" "+HEADER_SUFFIX: projection
           for projection, prefix in (("x", "x"), ("pt", "p_T (GeV/c)"), ("z", "z"))}
SECTION = re.compile(r"(Collins|Sivers) Asymmetry: (Positive|Negative) Pions\s+- full range")
ROW = re.compile(rf"\s*({NUMBER})\s+\[\s*({NUMBER})\s*,\s*({NUMBER}|inf\.)\s*\)"
                 rf"\s*({NUMBER})\s*\+/-\s*({NUMBER})\s*\|([^|]+)\|\s*")


def primary_text(path, spec):
    """Read one bounded, locked regular member without extracting the archive."""
    if file_hash(path) != spec["sha256"]:
        raise ValueError("changed primary COMPASS proton archive")
    with tarfile.open(path, mode="r:*") as archive:
        members = archive.getmembers()
        if len(members) > 256 or sum(m.size for m in members) > 32_000_000:
            raise ValueError("oversized primary COMPASS proton archive")
        names = set()
        for member in members:
            name = PurePosixPath(member.name)
            if (name.is_absolute() or ".." in name.parts or "\\" in member.name or
                    not (member.isdir() or member.isfile()) or member.name in names):
                raise ValueError("unsafe or duplicate primary COMPASS proton archive member")
            names.add(member.name)
        chosen = [m for m in members if m.name == spec["member"] and m.isfile()]
        if len(chosen) != 1 or chosen[0].size > 2_000_000:
            raise ValueError("missing/oversized primary COMPASS proton member")
        data = archive.extractfile(chosen[0]).read()
    if hashlib.sha256(data).hexdigest() != spec["member_sha256"]:
        raise ValueError("changed primary COMPASS proton member")
    return data.decode("utf-8")


def decimal_number(literal):
    if re.fullmatch(NUMBER, literal) is None:
        raise ValueError("invalid primary COMPASS proton number")
    result = Decimal(literal)
    if not result.is_finite():
        raise ValueError("nonfinite primary COMPASS proton number")
    return result


def source_row(match, line, projection):
    coordinate, low, high, value, sigma, means = match.groups()
    literals = dict(zip(MEANS, means.split(), strict=True))
    numbers = {key: decimal_number(literal) for key, literal in literals.items()}
    c, lower, a, error = map(decimal_number, (coordinate, low, value, sigma))
    upper = None if high == "inf." else decimal_number(high)
    numeric = [c, lower, a, error, *numbers.values()] + ([] if upper is None else [upper])
    if (any(not math.isfinite(float(v)) for v in numeric) or
            lower < 0 or not lower <= c or (upper is not None and not c < upper) or
            (upper is None and projection != "pt") or error <= 0 or
            numbers["Q2"] <= 0 or numbers["W"] <= 0 or numbers["pT"] < 0 or
            not all(0 < numbers[key] < 1 for key in ("x", "y", "z")) or
            (projection in ("x", "z") and (c <= 0 or upper > 1))):
        raise ValueError("invalid primary COMPASS proton bin/mean/uncertainty")
    return {"source_line": line, "coordinate": c, "value": a, "Sigma": error,
            "coordinate_literal": coordinate, "value_literal": value, "Sigma_literal": sigma,
            "means": literals, "low": lower, "high": upper,
            "bin": {"quantity": QUANTITIES[projection], "low": float(lower),
                    "high": None if upper is None else float(upper), "upper_unbounded": upper is None,
                    "low_literal": low, "high_literal": high, "interval": "left_closed_right_open"}}


def parse_tables(text, selection):
    if (selection["selected_harmonic"] != "Collins" or
            selection["harmonics_in_member"] != ["Collins", "Sivers"] or
            selection["hadrons"] != ["pi+", "pi-"] or set(selection["projections"]) != set(COORDINATES)):
        raise ValueError("unreviewed primary COMPASS proton selection")
    groups = {}
    section = current = None
    waiting_header = False
    for line_number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        match = SECTION.fullmatch(stripped)
        if match:
            if waiting_header or (current is not None and not groups[current]):
                raise ValueError("incomplete primary COMPASS proton section")
            section = (match[1], "pi+" if match[2] == "Positive" else "pi-")
            waiting_header = True
            current = None
            continue
        header = HEADERS.get(" ".join(stripped.split()))
        if header is not None:
            if not waiting_header or section is None:
                raise ValueError("ambiguous primary COMPASS proton header")
            current = (*section, header)
            if current in groups:
                raise ValueError("duplicate primary COMPASS proton section")
            groups[current] = []
            waiting_header = False
            continue
        match = ROW.fullmatch(line)
        if match is None or waiting_header or current is None:
            raise ValueError("changed primary COMPASS proton section/header/row")
        row = source_row(match, line_number, current[2])
        previous = groups[current][-1] if groups[current] else None
        if previous and (row["coordinate"] <= previous["coordinate"] or previous["high"] is None or
                         row["low"] < previous["high"]):
            raise ValueError("duplicate/unordered primary COMPASS proton coordinate or bin")
        groups[current].append(row)
    expected = {(kind, hadron, projection): count for kind in selection["harmonics_in_member"]
                for hadron in selection["hadrons"] for projection, count in selection["projections"].items()}
    if (waiting_header or set(groups) != set(expected) or len(groups) != selection["source_groups"] or
            any(len(groups[key]) != count for key, count in expected.items())):
        raise ValueError("incomplete primary COMPASS proton groups/rows")
    return {(hadron, projection): {"harmonic": "Collins", "hadron": hadron, "projection": projection,
            "points": rows} for (kind, hadron, projection), rows in groups.items() if kind == "Collins"}


def read_primary(path, spec, contract):
    if contract is None:
        raise ValueError("reviewed primary COMPASS proton contract required")
    source = contract["primary_source"]
    if (spec["publication"] != contract["publication"] or
            any(spec[key] != source[field] for key, field in
                (("sha256", "archive_sha256"), ("member", "member"), ("member_sha256", "member_sha256")))):
        raise ValueError("mismatched primary COMPASS proton source identity")
    views = parse_tables(primary_text(path, spec), contract["selection"])
    return {**{key: spec[key] for key in ("id", "publication", "sha256", "member", "member_sha256")},
            "contract_id": contract["id"], "harmonic": "Collins", "target": "proton", "views": views}


def legacy_decimal(record, key):
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("invalid legacy COMPASS proton numeric field: "+key)
    return Decimal(str(value))


def compare_row(record, reference, projection):
    """Check three exact cells; unequal means are evidence, never a rounded PASS."""
    for key, field in ((COORDINATES[projection], "coordinate"), ("value", "value"), ("error_u", "Sigma")):
        if legacy_decimal(record, key) != reference[field]:
            raise ValueError("secondary COMPASS proton row disagrees with primary numerical table: "+key)
    return {key: {"legacy_decimal": str(legacy_decimal(record, key)), "primary_literal": literal,
                  "exact_decimal_equal": legacy_decimal(record, key) == decimal_number(literal),
                  "unequal_mean_verified": False}
            for key, literal in reference["means"].items()}
