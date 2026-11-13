"""Pin and parse source expressions, without importing legacy numerical programs."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

import requests

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / 'config/star-born-contract.json'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def source_expression(value):
    # Only serialize syntax. The independent Julia reader accepts arithmetic
    # and literal index lookups; it never evals Python or Wolfram programs.
    return ast.unparse(value)


def legacy_records(text):
    functions = {node.name: node for node in ast.parse(text).body if isinstance(node, ast.FunctionDef)}
    coefficients = {}
    for name, target in [('get_Hupol', 'Hupol'), ('get_HTffa', 'HTffa')]:
        records = {}
        for node in functions[name].body:
            if isinstance(node, ast.Assign) and len(node.targets) == 1:
                key = node.targets[0]
                if isinstance(key, ast.Subscript) and isinstance(key.value, ast.Name) and key.value.id == target:
                    assert isinstance(key.slice, ast.Constant) and type(key.slice.value) is int
                    records[str(key.slice.value)] = source_expression(node.value)
        indices = range(1, 15) if name == 'get_Hupol' else range(13)
        assert set(records) == {str(i) for i in indices}
        coefficients[name] = records
    contractions = {}
    for name, accumulator in [('_get_dsig', 'upol'), ('_get_dsigST', 'ffcs')]:
        aliases, terms = {}, []
        for node in functions[name].body:
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and isinstance(node.value, ast.Subscript):
                parent = node.value
                while isinstance(parent, ast.Subscript):
                    parent = parent.value
                if isinstance(parent, ast.Name) and parent.id in ('f', 'ft', 'd', 'h', 'Hxxpz'):
                    aliases[node.targets[0].id] = source_expression(node.value)
            if isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name) and node.target.id == accumulator:
                assert isinstance(node.op, ast.Add)
                terms.append(source_expression(node.value))
        assert aliases and terms
        contractions[name] = {'aliases': aliases, 'terms': terms}
    return {'coefficients': coefficients, 'contractions': contractions}


def feyn_records(unpolarized, transfer):
    section = re.search(r'"UU"\s*->\s*<\|(.*?)\|>,\s*"LL"\s*->', unpolarized, re.S)
    assert section is not None
    pairs = re.findall(r'"([^"\n]+)"\s*->\s*(.*?)(?=,\s*"[^"\n]+"\s*->|\s*$)', section.group(1), re.S)
    UU = {key: value.strip() for key, value in pairs}
    assert len(UU) == len(pairs) == 29
    pattern = r'(\d+)\s*->\s*<\|"Process"\s*->\s*"([^"]+)"\s*,\s*"Born"\s*->\s*<\|"XX"\s*->\s*(.*?)\s*,\s*"NN"\s*->\s*(.*?)\|>'
    rows = re.findall(pattern, transfer, re.S)
    born = {number: {'process': process, 'XX': XX.strip(), 'NN': NN.strip()} for number, process, XX, NN in rows}
    assert len(rows) == len(born) == 7 and set(born) == {'1', '3', '6', '10', '11', '13', '17'}
    return {'UU_D_dimensional': UU, 'transfer_Born': born}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='.runtime/star-born-references')
    parser.add_argument('--external-root', type=Path, help='offline pinned external checkout parent; no upstream imports')
    args = parser.parse_args()
    output = Path(args.output).resolve()
    assert output.is_relative_to(ROOT) and not output.exists(), 'fresh repository-local output required'
    output.mkdir(parents=True)
    contract = json.loads(CONTRACT.read_text(encoding='utf-8'))
    sources = contract['sources']
    legacy = sources['legacy_code']
    feyn = sources['FeynFacet']
    assets = [('jam3d_dev_lib', legacy, legacy)] + [('FeynFacet', feyn, item) for item in feyn['assets']]
    records, texts = [], {}
    for directory, repository, asset in assets:
        path = PurePosixPath(asset['path'])
        assert not path.is_absolute() and '..' not in path.parts
        remote = urlsplit(repository['repository'])
        assert remote.hostname == 'github.com' and re.fullmatch('[0-9a-f]{40}', repository['commit'])
        url = 'https://raw.githubusercontent.com' + remote.path + '/' + repository['commit'] + '/' + str(path)
        if args.external_root:
            raw = (args.external_root / directory / str(path)).read_bytes()
            data = raw.replace(b'\r\n', b'\n')
            observed = {'mode': 'offline_checkout', 'checkout_sha256': sha(raw), 'normalization': 'CRLF to LF; verify exact Git blob'}
        else:
            response = requests.get(url, timeout=45)
            response.raise_for_status()
            data = response.content
            observed = {'mode': 'pinned_public_download', 'http_status': response.status_code}
        assert len(data) < 2_000_000 and sha(data) == asset['git_blob_sha256'], str(path)
        if 'bytes' in asset:
            assert len(data) == asset['bytes']
        destination = output / 'raw' / directory / str(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
        records.append({'repository': directory, 'commit': repository['commit'], 'source_path': str(path),
                        'url': url, 'raw_path': str(destination), 'bytes': len(data),
                        'git_blob_sha256': sha(data), **observed})
        texts[(directory, str(path))] = data.decode('utf-8')
    data = {'scope': 'pinned_source_expressions_not_native_physics_validation',
            'status': 'PREPARED_NOT_VALIDATED', 'production_ready': False, 'physical_extraction': False,
            'native_FeynFacet_calculation_executed': False, 'legacy_imports_or_pickle_deserialization_executed': False,
            'contract_sha256': sha(CONTRACT.read_bytes()), 'preparer_sha256': sha(Path(__file__).read_bytes()),
            'sources': records,
            'legacy': legacy_records(texts[('jam3d_dev_lib', legacy['path'])]),
            'FeynFacet': feyn_records(texts[('FeynFacet', feyn['assets'][0]['path'])], texts[('FeynFacet', feyn['assets'][1]['path'])])}
    target = output / 'reference-data.json'
    target.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'status': data['status'], 'sources': len(records), 'FeynFacet_UU_channels': 29,
                      'FeynFacet_transfer_classes': 7, 'output': str(target)}), flush=True)


if __name__ == '__main__':
    main()
