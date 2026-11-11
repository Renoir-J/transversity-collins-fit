"""Apply explicitly configured project administration using normal Git credentials.

Run manually. Credentials stay in memory and are never printed or written to disk.
The connected GitHub app lacks access to this account's new private repository.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "Renoir-J/transversity-collins-fit"
BASE = "https://api.github.com/repos/" + REPOSITORY


def authenticated_session():
    env = os.environ.copy()
    env.update(GIT_TERMINAL_PROMPT="0", GCM_INTERACTIVE="Never")
    result = subprocess.run(["git", "credential", "fill"], input="protocol=https\nhost=github.com\n\n",
                            capture_output=True, text=True, check=True, timeout=20, env=env)
    fields = dict(line.split("=", 1) for line in result.stdout.splitlines() if "=" in line)
    session = requests.Session()
    session.headers.update(Authorization="Bearer " + fields["password"], Accept="application/vnd.github+json")
    identity = session.get("https://api.github.com/user", timeout=30)
    identity.raise_for_status()
    if identity.json()["login"] != "Renoir-J":
        raise ValueError("authenticated GitHub account differs from project owner")
    return session


def request(session, method, route, **kwargs):
    response = session.request(method, BASE + route, timeout=30, **kwargs)
    response.raise_for_status()
    return response.json() if response.content else None


def bootstrap(session):
    settings = request(session, "PATCH", "", json={"has_issues": True, "has_wiki": False,
        "allow_squash_merge": True, "allow_merge_commit": False, "allow_rebase_merge": False,
        "delete_branch_on_merge": True})
    milestone_specs = [
        ("M0 — Reproducible workflow", "Julia engine, Linux/Snakemake closure, provenance and CI."),
        ("M1 — Experimental data", "Reviewed public source acquisition, bin conventions, covariance and overlap."),
        ("M2 — Physics validation", "KPSY reconstruction, evolution, bin-integrated adapters and independent benchmarks."),
        ("M3 — Global fits and products", "B0–B4, V1, uncertainty/coverage, holdouts and publication artifacts."),
        ("M4 — NLO matching", "FeynFacet hard functions with matching, subtraction and jet definitions.")]
    milestones = {m["title"]: m for m in request(session, "GET", "/milestones?state=all&per_page=100")}
    for title, description in milestone_specs:
        if title not in milestones:
            milestones[title] = request(session, "POST", "/milestones", json={"title": title, "description": description})
    label_specs = {"workflow": "1d76db", "physics": "5319e7", "data": "0075ca", "validation": "0e8a16",
                   "blocked:science": "b60205", "status:todo": "fbca04", "status:done": "0e8a16", "tracking": "c5def5"}
    labels = {label["name"] for label in request(session, "GET", "/labels?per_page=100")}
    for name, color in label_specs.items():
        if name not in labels:
            request(session, "POST", "/labels", json={"name": name, "color": color})
    issues = request(session, "GET", "/issues?state=all&per_page=100")
    existing = {issue["title"]: issue for issue in issues}
    specs = json.loads((ROOT / "config/project-tasks.json").read_text(encoding="utf-8"))
    created = []
    for spec in specs:
        issue = existing.get(spec["title"])
        if issue is None:
            issue = request(session, "POST", "/issues", json={"title": spec["title"], "body": spec["body"],
                "labels": spec["labels"], "assignees": ["Renoir-J"],
                "milestone": milestones[milestone_specs[spec["milestone"]][0]]["number"]})
        created.append({"number": issue["number"], "title": issue["title"], "url": issue["html_url"], "state": issue["state"]})
    state = {"repository": REPOSITORY, "private": settings["private"], "url": settings["html_url"],
             "merge_strategy": "squash_only", "delete_branch_on_merge": settings["delete_branch_on_merge"],
             "milestones": [{"number": m["number"], "title": m["title"], "url": m["html_url"]} for m in milestones.values()],
             "issues": created, "branch_protection": "NOT_YET_ATTEMPTED"}
    path = ROOT / "config/github-state.json"
    path.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"repository": REPOSITORY, "issues": len(created), "milestones": len(milestones)}))


def protect(session):
    policy = {"required_status_checks": {"strict": True, "contexts": ["validation"]},
              "enforce_admins": True, "required_pull_request_reviews": {"required_approving_review_count": 0,
                  "dismiss_stale_reviews": True, "require_code_owner_reviews": False}, "restrictions": None,
              "required_linear_history": True, "allow_force_pushes": False, "allow_deletions": False,
              "required_conversation_resolution": True}
    response = session.put(BASE + "/branches/main/protection", json=policy, timeout=30)
    state_path = ROOT / "config/github-state.json"
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if response.status_code == 200:
        state["branch_protection"] = {"enforced": True, "policy": policy}
    else:
        state["branch_protection"] = {"enforced": False, "http_status": response.status_code,
            "message": response.json().get("message"), "desired_policy": policy}
    state_path.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(state["branch_protection"], ensure_ascii=False))
    if response.status_code not in (200, 403):
        response.raise_for_status()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["bootstrap", "protect"])
    args = parser.parse_args()
    session = authenticated_session()
    bootstrap(session) if args.action == "bootstrap" else protect(session)


if __name__ == "__main__":
    main()

