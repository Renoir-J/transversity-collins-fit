"""Explicitly synthetic primary/secondary text; no experimental data fixtures."""
import copy
import hashlib
import io
import json
import tarfile

import pytest
from openpyxl import Workbook

from collinsflow import compass_proton_primary as cpp
from collinsflow import historical as hist
from collinsflow.provenance import file_hash


def synthetic_text():
    lines = []
    for harmonic in ("Collins", "Sivers"):
        for charge in ("Positive", "Negative"):
            for projection in ("x", "pt", "z"):
                lines.extend([f"{harmonic} Asymmetry: {charge} Pions - full range",
                              next(h for h, p in cpp.HEADERS.items() if p == projection)])
                coordinates = {"x": [".10", ".20"], "pt": [".20", ".90"], "z": [".30", ".60"]}[projection]
                bounds = {"x": [(".05", ".15"), (".15", ".25")],
                          "pt": [(".10", ".30"), (".30", "inf.")],
                          "z": [(".20", ".40"), (".40", ".80")]}[projection]
                for i, (coordinate, (low, high)) in enumerate(zip(coordinates, bounds)):
                    means = dict(Q2="3.0", y=".30", pT=".20", x=".12", z=".40", W="4.0")
                    means[cpp.COORDINATES[projection]] = coordinate
                    value = ".0125" if harmonic == "Collins" else ".075"
                    if charge == "Negative":
                        value = "-"+value
                    lines.append(f"{coordinate} [ {low}, {high} ) {value} +/- .025 | " +
                                 " ".join(means.values())+" |")
    return "\n".join(lines)+"\n"


def fixture(tmp_path, text=None, extra=None):
    raw = (synthetic_text() if text is None else text).encode()
    path = tmp_path / 'synthetic-primary.tar.gz'
    member = 'synthetic/pions.txt'
    with tarfile.open(path, 'w:gz') as archive:
        info = tarfile.TarInfo(member)
        info.size = len(raw)
        archive.addfile(info, io.BytesIO(raw))
        if extra is not None:
            archive.addfile(extra, io.BytesIO(b'x'*extra.size) if extra.isfile() else None)
    spec = dict(id='synthetic-primary', publication='synthetic-v1', sha256=file_hash(path),
                member=member, member_sha256=hashlib.sha256(raw).hexdigest())
    policy = dict(id='SYNTHETIC-PARTIAL-PRIMARY', publication='synthetic-v1',
        primary_source=dict(archive_sha256=spec['sha256'], member=member, member_sha256=spec['member_sha256']),
        selection=dict(selected_harmonic='Collins', harmonics_in_member=['Collins', 'Sivers'],
                       hadrons=['pi+', 'pi-'], projections=dict(x=2, pt=2, z=2), source_groups=12))
    return path, spec, policy


def secondary(tmp_path, primary, projection='x', mutation=None):
    path = tmp_path/'synthetic-secondary.xlsx'
    headers = ['Q2', 'y', 'pT', 'x', 'z', 'W', 'obs', 'value', 'error_u', 'target', 'hadron', 'col', 'dependence']
    rows = []
    for reference in primary['views'][('pi+', projection)]['points']:
        row = {k: float(v) for k, v in reference['means'].items()}
        row.update(Q2=3.0001, **{cpp.COORDINATES[projection]: float(reference['coordinate'])},
                   obs='AUTcollins', value=float(reference['value']), error_u=float(reference['Sigma']),
                   target='proton', hadron='pi+', col='compass', dependence=projection)
        rows.append(row)
    if mutation:
        mutation(rows)
    book = Workbook()
    sheet = book.active
    sheet.append(headers)
    for row in rows:
        sheet.append([row[h] for h in headers])
    book.save(path)
    book.close()
    spec = dict(id='synthetic-secondary', publication='synthetic-v1', sha256=file_hash(path),
        version='explicitly synthetic fixture', url='https://example.invalid/synthetic-secondary',
        adapter=dict(kind='sidis', experiment='compass', target='proton', hadron='pi+', projection=projection,
            source_obs='AUTcollins', input_unit='fraction', sheet='Sheet', expected_raw_rows=2,
            expected_points=1 if projection == 'pt' else 2, expected_headers=headers,
            uncertainty_layout='reported_unspecified', depolarization_status='publication_review_pending',
            primary_crosscheck_kind='COMPASS_proton_partial', primary_crosscheck_id=primary['id']))
    return path, spec


def view_fixture(tmp_path, projection='x', mutation=None):
    archive, spec, policy = fixture(tmp_path)
    primary = cpp.read_primary(archive, spec, policy)
    path, source = secondary(tmp_path, primary, projection, mutation)
    contract = {'id': 'SYNTHETIC-HISTORICAL', 'SIDIS_projection_cut': {'upper_GeV': .8}}
    view = hist.parse_historical_xlsx(path, source, contract, {},
                                    proton_tables={primary['id']: primary}, proton_policy=policy)
    return view, path, source, primary, policy, contract


def test_primary_selects_labelled_signed_collins_not_sivers_and_retains_unbounded_bin(tmp_path):
    path, spec, policy = fixture(tmp_path)
    result = cpp.read_primary(path, spec, policy)
    assert len(result['views']) == 6
    assert result['views'][('pi-', 'x')]['points'][0]['value'] == cpp.Decimal('-.0125')
    last = result['views'][('pi+', 'pt')]['points'][-1]['bin']
    assert last['high'] is None and last['upper_unbounded'] is True and last['high_literal'] == 'inf.'
    json.dumps(last, allow_nan=False)


def test_partial_crosscheck_keeps_numbers_unknown_errors_and_unequal_means(tmp_path):
    view, *_ = view_fixture(tmp_path)
    check = view['primary_numeric_crosscheck']
    assert check['numeric_cells'] == 6 and check['raw_points'] == 2 and check['mean_cells_compared'] == 12
    assert check['exact_mean_cells'] == 10 and check['unequal_mean_cells_unverified'] == 2
    for point in view['points']:
        assert point['central_value'] == .0125 and point['reported_error_u'] == .025
        assert point['statistical_error'] is None and point['systematic_components'] == {}
        assert point['averages']['Q2_GeV2'] == 3.0001
        mean = point['primary_reported_means']['Q2']
        assert mean == dict(legacy_decimal='3.0001', primary_literal='3.0', exact_decimal_equal=False,
                            unequal_mean_verified=False)
    assert view['depolarization_removed'] is None
    assert view['production_ready'] is view['historical_fit_ready'] is False
    assert any('composition' in b for b in view['blockers'])
    assert not any(b == 'secondary transcription primary numerical crosscheck pending' for b in view['blockers'])


def test_all_raw_rows_crosschecked_before_cut_and_excluded_open_bin_kept(tmp_path):
    view, *_ = view_fixture(tmp_path, 'pt')
    assert len(view['points']) == len(view['excluded']) == 1
    assert view['primary_numeric_crosscheck']['raw_points'] == 2
    excluded = view['excluded'][0]
    assert excluded['primary_projection_bin']['upper_unbounded'] is True
    assert excluded['primary_reported_means']['pT']['exact_decimal_equal'] is True
    with pytest.raises(ValueError, match='disagrees with primary'):
        view_fixture(tmp_path, 'pt', lambda rows: rows[-1].update(value=.0126))


@pytest.mark.parametrize('column,value', [('value', -.0125), ('error_u', .025001), ('x', .100001)])
def test_exact_fields_reject_sign_changes_and_additional_digits(tmp_path, column, value):
    with pytest.raises(ValueError, match='disagrees with primary|coordinate'):
        view_fixture(tmp_path, mutation=lambda rows: rows[0].update({column: value}))


def test_primary_mean_inequality_is_recorded_without_display_interval_acceptance(tmp_path):
    view, *_ = view_fixture(tmp_path, mutation=lambda rows: rows[0].update(Q2=30.))
    mean = view['points'][0]['primary_reported_means']['Q2']
    assert mean['exact_decimal_equal'] is mean['unequal_mean_verified'] is False
    assert view['points'][0]['averages']['Q2_GeV2'] == 30.
    assert view['historical_fit_ready'] is False


@pytest.mark.parametrize('old,new', [
    ('Collins Asymmetry:', 'Sivers Asymmetry:'), ('Positive Pions', 'Positive Kaons'),
    ('full range', 'selected range'), ('Asym +/- Sigma', 'Asym +/- Stat'),
    ('<Q^2> (GeV/c)^2', '<Q^2> (MeV/c)^2'), ('<y> <p_T>', '<p_T> <y>'),
    ('.0125 +/- .025', '.0125 +/- 0'), ('.10 [ .05, .15 )', '.10 [ .15, .20 )'),
    ('.10 [ .05, .15 )', '.10 [ .05, inf. )'), ('3.0 .30 .20 .10 .40 4.0', '3.0 1.30 .20 .10 .40 4.0'),
    ('3.0 .30 .20 .10 .40 4.0', 'NaN .30 .20 .10 .40 4.0')])
def test_primary_identity_units_and_domains_fail_closed(tmp_path, old, new):
    text = synthetic_text()
    assert old in text
    path, spec, policy = fixture(tmp_path, text.replace(old, new))
    with pytest.raises(ValueError):
        cpp.read_primary(path, spec, policy)


def test_missing_or_duplicate_primary_row_rejected(tmp_path):
    text = synthetic_text()
    line = next(line for line in text.splitlines() if '+/- .025' in line)
    for replacement in ('', line+'\n'+line):
        path, spec, policy = fixture(tmp_path, text.replace(line, replacement, 1))
        with pytest.raises(ValueError):
            cpp.read_primary(path, spec, policy)


@pytest.mark.parametrize('field,value', [('harmonic', 'Sivers'), ('target', 'deuteron'),
                                       ('publication', 'different-version'), ('contract_id', 'unreviewed')])
def test_equal_numbers_do_not_certify_wrong_primary_identity(tmp_path, field, value):
    _, path, source, primary, policy, contract = view_fixture(tmp_path)
    primary[field] = value
    with pytest.raises(ValueError, match='identity'):
        hist.parse_historical_xlsx(path, source, contract, {},
                                  proton_tables={primary['id']: primary}, proton_policy=policy)


@pytest.mark.parametrize('field,value', [('harmonic', 'Sivers'), ('hadron', 'pi-'), ('projection', 'pt')])
def test_equal_numbers_do_not_certify_wrong_primary_view(tmp_path, field, value):
    _, path, source, primary, policy, contract = view_fixture(tmp_path)
    primary['views'][('pi+', 'x')][field] = value
    with pytest.raises(ValueError, match='charge/projection'):
        hist.parse_historical_xlsx(path, source, contract, {},
                                  proton_tables={primary['id']: primary}, proton_policy=policy)


def test_missing_primary_or_removed_crosscheck_flag_cannot_skip_required_checks(tmp_path):
    _, path, source, primary, policy, contract = view_fixture(tmp_path)
    with pytest.raises(ValueError, match='identity'):
        hist.parse_historical_xlsx(path, source, contract, {}, proton_policy=policy)
    source['adapter'].pop('primary_crosscheck_kind')
    with pytest.raises(ValueError, match='identity'):
        hist.parse_historical_xlsx(path, source, contract, {},
                                  proton_tables={primary['id']: primary}, proton_policy=policy)


@pytest.mark.parametrize('unsafe', ['../escape', '/absolute', 'bad\\name', 'symlink', 'hardlink', 'device', 'duplicate'])
def test_archive_members_not_extracted_and_unsafe_members_rejected(tmp_path, unsafe):
    info = tarfile.TarInfo('synthetic/pions.txt' if unsafe == 'duplicate' else unsafe)
    if unsafe == 'symlink':
        info.type = tarfile.SYMTYPE
        info.linkname = '/outside'
    elif unsafe == 'hardlink':
        info.type = tarfile.LNKTYPE
        info.linkname = 'synthetic/pions.txt'
    elif unsafe == 'device':
        info.type = tarfile.CHRTYPE
    path, spec, policy = fixture(tmp_path, extra=info)
    with pytest.raises(ValueError, match='unsafe or duplicate'):
        cpp.read_primary(path, spec, policy)
    assert not (tmp_path/'synthetic').exists()


def test_primary_archive_and_member_hashes_are_bound(tmp_path):
    path, spec, policy = fixture(tmp_path)
    original = path.read_bytes()
    path.write_bytes(original+b'x')
    with pytest.raises(ValueError, match='changed primary.*archive'):
        cpp.read_primary(path, spec, policy)
    path.write_bytes(original)
    spec['member_sha256'] = policy['primary_source']['member_sha256'] = '0'*64
    with pytest.raises(ValueError, match='changed primary.*member'):
        cpp.read_primary(path, spec, policy)


def test_oversized_primary_member_is_rejected_before_parsing(tmp_path):
    path, spec, policy = fixture(tmp_path, 'x'*2_000_001)
    with pytest.raises(ValueError, match='oversized primary.*member'):
        cpp.read_primary(path, spec, policy)


def test_excessive_archive_member_count_is_rejected(tmp_path):
    path = tmp_path/'too-many-synthetic-members.tar.gz'
    with tarfile.open(path, 'w:gz') as archive:
        for i in range(257):
            archive.addfile(tarfile.TarInfo(f'synthetic-{i}.txt'), io.BytesIO())
    spec = dict(sha256=file_hash(path), member='synthetic-0.txt', member_sha256=hashlib.sha256(b'').hexdigest())
    with pytest.raises(ValueError, match='oversized primary.*archive'):
        cpp.primary_text(path, spec)


@pytest.mark.parametrize('mutation', ['statistical', 'means_verified', 'systematic'])
def test_partial_candidate_cannot_claim_certified_means_or_uncertainties(tmp_path, mutation):
    view, *_, contract = view_fixture(tmp_path)
    contract['expected'] = dict(views=1, points=2, SIDIS_views=1, SIDIS_points=2, SIA_views=0, SIA_points=0)
    candidate = dict(schema_version=1, scope=hist.SCOPE, production_ready=False, physical_extraction=False,
                     historical_fit_ready=False, contract_id=contract['id'], views=[view])
    hist.validate_historical_candidate(candidate, contract)
    changed = copy.deepcopy(candidate)
    point = changed['views'][0]['points'][0]
    if mutation == 'statistical':
        point['statistical_error'] = .025
    elif mutation == 'systematic':
        point['systematic_components']['invented'] = .01
    else:
        point['primary_reported_means']['Q2']['unequal_mean_verified'] = True
    with pytest.raises(ValueError, match='cannot certify'):
        hist.validate_historical_candidate(changed, contract)
