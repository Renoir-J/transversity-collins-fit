"""Cache public primary-source papers for this research audit; never executes them."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import requests
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
PAPERS = {
    "1505.05589": "1505.05589v1",
    "2002.08384": "2002.08384v2",
    "2205.00999": "2205.00999v2",
    "2205.11800": "2205.11800v2",
    "2507.16355": "2507.16355v2",
    "2610.31745": "2610.31745v1",
    "2401.00309": "2401.00309v1",
    "1909.01857": "1909.01857v1",
    "1707.00913": "1707.00913v1",
    "2506.21959": "2506.21959v1",
    "2609.35958": "2609.35958v1",
    "2603.26990": "2603.26990v1",
    "2610.16733": "2610.16733v1",
    "2610.36224": "2610.36224v1",
    "2603.22471": "2603.22471v1",
    "2512.19799": "2512.19799v2",
    "1909.00485": "1909.00485v2",
    "2007.07755": "2007.07755v1",
    "1507.06824": "1507.06824v1",
}


def fetch(item):
    key, version = item
    url = f"https://arxiv.org/html/{version}"
    record = {"arxiv": key, "version": version, "url": url}
    try:
        target = ROOT / "references" / "papers"
        target.mkdir(parents=True, exist_ok=True)
        cached = target / f"{version}.html"
        if cached.exists():
            payload = cached.read_bytes()
        else:
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            payload = response.content
            cached.write_bytes(payload)
        soup = BeautifulSoup(payload, "html.parser")
        article = soup.find("article") or soup.body or soup
        # Paragraphs, headings, displays and tables, without duplicate table cells.
        nodes = article.select("h1,h2,h3,h4,h5,p,table,figure")
        chunks = []
        for node in nodes:
            if any(p.name in {"p", "table", "figure"} for p in node.parents if p is not article):
                continue
            chunks.append(node.get_text(" ", strip=True))
        text = "\n\n".join(chunks)
        (target / f"{version}.txt").write_text(text, encoding="utf-8")
        record.update(sha256=hashlib.sha256(payload).hexdigest(),
                      bytes=len(payload), text_chars=len(text), status="cached")
    except Exception as error:
        record.update(status="unavailable", error=str(error))
    return record


if __name__ == "__main__":
    records = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        for future in as_completed([executor.submit(fetch, p) for p in PAPERS.items()]):
            item = future.result()
            records.append(item)
            print(item["arxiv"], item["status"], item.get("text_chars", item.get("error")), flush=True)
    out = ROOT / "references" / "literature-cache.json"
    out.write_text(json.dumps({"retrieved_at": datetime.now(timezone.utc).isoformat(),
                               "papers": sorted(records, key=lambda x: x["arxiv"])},
                              indent=2), encoding="utf-8")
