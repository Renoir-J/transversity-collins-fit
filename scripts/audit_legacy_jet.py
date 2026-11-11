"""Probe two legacy jet-adapter edge cases with synthetic backends, not physics data.

Only the AST of RESIDUALS._get_theory is executed. No upstream imports, grids,
pickled replicas, network requests, or real cross-section calculations are used.
"""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.integrate import fixed_quad

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "external/jam3d_dev_lib/obslib/AN_pp/residuals.py"


def audit():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "RESIDUALS")
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "_get_theory")
    namespace = {"np": np, "fixed_quad": fixed_quad, "conf": {"unpol PDF FF fixed": True}}
    calls = {"denominator": 0}

    def denominator(*args, **kwargs):
        calls["denominator"] += 1
        return 2.0

    namespace["AN_theory_hj"] = SimpleNamespace(get_sig=denominator, get_sigST=lambda *a, **k: 1.0)
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(SOURCE), "exec"), namespace)
    compute = namespace["_get_theory"]
    common = {"obs": ["AUTCollins"], "col": ["STAR"], "target": ["p"],
              "hadron": ["pi+"], "jperp": [0.3], "z": [0.25], "rs": [510.0]}

    # Both dimensions are fixed: the upstream adapter has no matching branch.
    fixed = SimpleNamespace(tabs={1: dict(common, pT=[20.0], eta=[0.4])}, storage={})
    try:
        compute(fixed, (1, 0))
        fixed_result = {"reproduced": False}
    except (KeyError, UnboundLocalError) as error:
        fixed_result = {"reproduced": True, "exception": type(error).__name__, "detail": str(error)}

    # Integrate pT, hold eta fixed; the 4-field membership test cannot hit a 5-field key.
    integrated = SimpleNamespace(tabs={1: dict(common, pTmin=[10.0], pTmax=[20.0], eta=[0.4])}, storage={})
    first = compute(integrated, (1, 0))
    after_first = calls["denominator"]
    second = compute(integrated, (1, 0))
    cache_result = {"reproduced": calls["denominator"] == 2 * after_first and after_first > 0,
                    "calls_after_first_evaluation": after_first,
                    "calls_after_two_identical_evaluations": calls["denominator"],
                    "calls_expected_with_cache": after_first,
                    "outputs_equal": bool(first == second),
                    "stored_key_lengths": sorted({len(k) for k in integrated.storage})}
    return {"scope": "synthetic adapter-control-flow probes; no validation of published fits",
            "source": SOURCE.relative_to(ROOT).as_posix(),
            "fixed_pT_fixed_eta_missing_branch": fixed_result,
            "denominator_cache_membership_mismatch": cache_result}


if __name__ == "__main__":
    result = audit()
    destination = ROOT / "artifacts/legacy-jet-audit.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not all(result[key]["reproduced"] for key in
               ["fixed_pT_fixed_eta_missing_branch", "denominator_cache_membership_mismatch"]):
        raise SystemExit("Upstream behavior changed; review this audit before reusing the findings.")
