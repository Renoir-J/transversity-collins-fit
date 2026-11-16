"""Explicitly synthetic schema/ownership controls, never experimental validation."""
import copy
from pathlib import Path

import pytest

from collinsflow import historical_sia_points as sia
from collinsflow.provenance import atomic_json, read_json

ROOT = Path(__file__).resolve().parents[1]


def synthetic_candidate():
    contract = read_json(ROOT/'config/historical-sia-point-contract.json')
    initial = read_json(ROOT/'config/initial-contract.json')
    primary = contract['primary_definition']
    views = []
    for selected in contract['selection']['source_views']:
        rows = []
        for i in range(9):
            p = .05+.1*i
            rows.append(dict(id=f"{selected['id']}:r{i+2}",
                averages=dict(z1=.3, z2=.4, Q2_GeV2=112.36, pT_GeV=p),
                central_value=(-1)**i*.0125, statistical_error=.05,
                systematic_components=dict(reported_systematic=.002),
                uncertainty_source_columns=dict(statistical='stat_u', reported_systematic='sys_u'),
                locator=dict(sheet='explicitly-synthetic-unit-fixture', excel_row=i+2),
                bins=dict(pT_GeV=[max(0., p-.04), p+.04]), primary_bins=dict(pT=[max(0., p-.04), p+.04]),
                primary_locator=dict(publication='1309.5278v2', member='paper.tex', table='tab:totPt',
                    table_row=i+1, source_line=100+i),
                angular_source=dict(column='S2/1+C2', value=.6, average_type='mean_of_angular_ratio', frame='hadron2',
                    averaging_measure='explicitly synthetic mathematical measure; no experimental acceptance',
                    measure_status='joint_acceptance_unresolved')))
        crosscheck = dict(publication='1309.5278v2', frame='hadron2', sha256=primary['archive_sha256'],
            member=primary['member'], member_sha256=primary['member_sha256'], table='tab:totPt', raw_points=9,
            numeric_cells=81, compared_columns=primary['checked_columns'])
        views.append(dict(id=selected['id'], source_id=selected['id'], process='SIA', projection='pT',
            ratio=selected['ratio'], frame='hadron2', publication='1309.5278', source_sha256=selected['sha256'],
            source_observable=f"A{selected['ratio']}-0-PT", input_units='percent', output_units='fraction',
            overlap_group='babar_1309.5278', production_ready=False, historical_fit_ready=False,
            blockers=['explicitly synthetic schema fixture; no raw/primary replay or data admission'],
            primary_numeric_crosscheck=crosscheck, points=rows, excluded=[]))
    return dict(scope='historical_benchmark_candidate', production_ready=False, physical_extraction=False,
                historical_fit_ready=False, fixture_scope='explicitly_synthetic_unit_fixture', views=views), contract, initial


def test_preserve_signed_fractions_error_components_angular_measure_and_point_ownership():
    candidate, contract, initial = synthetic_candidate()
    before = copy.deepcopy(candidate)
    points = sia.request_points(candidate, contract, initial)
    assert len(points) == 18 and [p['ratio'] for p in points] == ['UL']*9+['UC']*9
    assert points[1]['measured'] == dict(central_value=-.0125, statistical_error=.05,
                                       systematic_components=dict(reported_systematic=.002))
    assert all(p['source_angular']['average_type'] == 'mean_of_angular_ratio' for p in points)
    assert all(p['legacy_Q2_status'] == 'source_literal_not_checked_by_selected_primary_table' for p in points)
    assert all(p['kinematics']['Q2_GeV2'] == 112.36 for p in points)
    points[0]['measured']['systematic_components']['reported_systematic'] = 100.
    points[0]['source_angular']['averaging_measure'] = 'changed'
    assert candidate == before


@pytest.mark.parametrize('scope', ['production', 'synthetic_closure', 'experimental_candidate'])
def test_other_scopes_do_not_enter_real_source_preparation(scope):
    candidate, contract, initial = synthetic_candidate()
    candidate['scope'] = scope
    with pytest.raises(ValueError, match='scope'):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('flag', sia.FALSE_FLAGS)
def test_contract_cannot_self_authorize_admission_or_coverage(flag):
    candidate, contract, initial = synthetic_candidate()
    contract[flag] = True
    with pytest.raises(ValueError, match='approve'):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('field,value', [
    ('frame', 'thrust'), ('projection', 'z1_z2'), ('publication', 'unreviewed-version'), ('ratio', 'UC'),
    ('source_observable', 'AUTcollins'), ('source_sha256', '0'*64), ('output_units', 'percent'),
    ('input_units', 'fraction'), ('overlap_group', ''), ('production_ready', True)])
def test_wrong_frame_ratio_version_units_or_admission_rejects(field, value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0][field] = value
    with pytest.raises(ValueError, match='selection'):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('mutation', ['omit_view', 'omit_point', 'duplicate', 'reorder', 'primary_columns', 'mean_geometry'])
def test_source_selection_and_exact_geometry_are_not_inferred_or_rounded(mutation):
    candidate, contract, initial = synthetic_candidate()
    view = candidate['views'][0]
    if mutation == 'omit_view':
        candidate['views'].pop()
    elif mutation == 'omit_point':
        view['points'].pop()
    elif mutation == 'duplicate':
        view['points'][1] = copy.deepcopy(view['points'][0])
    elif mutation == 'reorder':
        view['points'].reverse()
    elif mutation == 'primary_columns':
        view['primary_numeric_crosscheck']['compared_columns'] = ['pT', 'value']
    else:
        view['points'][0]['averages']['z1'] += 1e-10
    with pytest.raises(ValueError):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('field,value', [('average_type', 'ratio_of_angular_means'), ('value', True),
    ('value', float('nan')), ('frame', 'thrust'), ('measure_status', 'acceptance_validated')])
def test_angular_definition_and_acceptance_limit_are_source_fields(field, value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0]['points'][0]['angular_source'][field] = value
    with pytest.raises(ValueError):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('field,value', [('Q2_GeV2', 10.58**2), ('Q2_GeV2', True),
    ('z1', 0.), ('pT_GeV', float('inf')), ('pT_GeV', 5.)])
def test_no_legacy_energy_replacement_or_invalid_means(field, value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0]['points'][0]['averages'][field] = value
    with pytest.raises(ValueError):
        sia.request_points(candidate, contract, initial)


def test_existing_request_bytes_survive_rejected_overwrite(tmp_path):
    output = tmp_path/'request.json'
    output.write_text('preserved original evidence', encoding='utf-8')
    with pytest.raises(ValueError, match='fresh'):
        sia.prepare_request('absent-manifest', 'absent-cache', 'absent-candidates', output)
    assert output.read_text(encoding='utf-8') == 'preserved original evidence'


def test_existing_diagnostic_artifacts_survive_before_numerical_dispatch(tmp_path, monkeypatch):
    monkeypatch.setattr(sia, 'verify_request', lambda p: dict(source_preparation=dict(manifest=str(ROOT/'config/historical-sources.json'))))
    output = tmp_path/'diagnostic'
    output.mkdir()
    sentinel = output/'predictions.log'
    sentinel.write_text('preserved interrupted output', encoding='utf-8')
    with pytest.raises(ValueError, match='fresh'):
        sia.run_point_job(tmp_path/'request.json', tmp_path/'native.json', tmp_path/'parton.json', output)
    assert sentinel.read_text(encoding='utf-8') == 'preserved interrupted output'


def test_cached_production_scope_never_dispatches_an_auditor(tmp_path, monkeypatch):
    monkeypatch.setattr(sia, 'verify_request', lambda p: dict(source_preparation=dict(manifest=str(ROOT/'config/historical-sources.json'))))
    monkeypatch.setattr(sia.subprocess, 'run', lambda *a, **k: pytest.fail('numerical dispatch after wrong-scope cache'))
    atomic_json(tmp_path/'invocation.json', dict(status='PASS', scope='production'))
    with pytest.raises(ValueError, match='scope'):
        sia.verify_point_job(tmp_path/'request.json', tmp_path)
