"""Synthetic representation controls; no primary or numerical admission."""
import copy
from pathlib import Path

import pytest

from collinsflow import historical_sia_integrals as sia
from collinsflow.provenance import atomic_json, file_hash, read_json

ROOT = Path(__file__).resolve().parents[1]


def synthetic_candidate():
    contract = read_json(ROOT/'config/historical-sia-integral-contract.json')
    initial = read_json(ROOT/'config/initial-contract.json')
    primary = contract['primary_definition']
    views = []
    for selected in contract['selection']['source_views']:
        systematic = 'syt_u' if selected['ratio'] == 'UL' else 'sys_u'
        rows = []
        for i in range(36):
            z1, z2 = .2+.1*(i//6), .2+.1*(i % 6)
            bins = dict(z1=[z1-.04, z1+.04], z2=[z2-.04, z2+.04])
            rows.append(dict(id=f"{selected['id']}:r{i+2}",
                averages=dict(z1=z1, z2=z2, Q2_GeV2=112.36, pT_GeV=.3+.001*i),
                central_value=(-1)**i*.0125, statistical_error=.005,
                systematic_components=dict(reported_systematic=.002),
                uncertainty_source_columns=dict(statistical='stat_u', reported_systematic=systematic),
                locator=dict(sheet='explicitly-synthetic-unit-fixture', excel_row=i+2),
                bins=copy.deepcopy(bins), primary_bins=bins,
                primary_locator=dict(publication='1309.5278v2', member='paper.tex', table=selected['primary_table'],
                    table_row=i+1, source_line=100+i),
                angular_source=dict(column='S2/1+C2', value=.6, average_type='mean_of_angular_ratio', frame='hadron2',
                    averaging_measure='explicitly synthetic mathematical measure; no experimental acceptance',
                    measure_status='joint_acceptance_unresolved')))
        cross = dict(publication='1309.5278v2', frame='hadron2', sha256=primary['archive_sha256'],
            member='paper.tex', member_sha256=primary['member_sha256'], table=selected['primary_table'],
            raw_points=36, numeric_cells=396,
            compared_columns=['S2/1+C2', 'pT', 'stat_u', systematic, 'value', 'z1', 'z1min', 'z1max', 'z2', 'z2min', 'z2max'])
        views.append(dict(id=selected['id'], source_id=selected['id'], source_sha256=selected['sha256'],
            process='SIA', projection='z1_z2', ratio=selected['ratio'], frame='hadron2', publication='1309.5278',
            source_observable=f"A{selected['ratio']}-0-PT-INT", input_units='percent', output_units='fraction',
            prediction_prescription='P-integrated source view; theory adapter pending', overlap_group='babar_1309.5278',
            production_ready=False, historical_fit_ready=False, blockers=['explicitly synthetic representation fixture'],
            primary_numeric_crosscheck=cross, points=rows, excluded=[]))
    return dict(scope='historical_benchmark_candidate', production_ready=False, physical_extraction=False,
                historical_fit_ready=False, fixture_scope='explicitly_synthetic_unit_fixture', views=views), contract, initial


def test_integrated_source_preserves_signed_values_and_P_metadata_separately():
    candidate, contract, initial = synthetic_candidate()
    before = copy.deepcopy(candidate)
    points = sia.request_points(candidate, contract, initial)
    assert len(points) == 72 and [p['ratio'] for p in points] == ['UC']*36+['UL']*36
    assert len({tuple(p['kinematics'].values()) for p in points}) == 36
    assert all(set(p['kinematics']) == {'z1', 'z2', 'Q2_GeV2'} for p in points)
    assert points[1]['measured']['central_value'] == -.0125
    assert points[36]['error_source_columns']['reported_systematic'] == 'syt_u'
    assert points[0]['mapped_P_bounds_GeV'] == [0., .2*3.5]
    assert points[0]['reported_pT_role'] == 'source_metadata_not_model_momentum_or_bound'
    candidate['views'][0]['points'][0]['averages']['pT_GeV'] += .01
    again = sia.request_points(candidate, contract, initial)
    assert again[0]['kinematics'] == points[0]['kinematics']
    assert again[0]['mapped_P_bounds_GeV'] == points[0]['mapped_P_bounds_GeV']
    assert again[0]['reported_pT_mean_GeV'] != points[0]['reported_pT_mean_GeV']
    points[0]['primary_bins']['z1'][0] = -1
    assert before['views'][0]['points'][0]['primary_bins'] == candidate['views'][0]['points'][0]['primary_bins']


@pytest.mark.parametrize('flag', sia.FALSE_FLAGS)
def test_no_self_admission(flag):
    candidate, contract, initial = synthetic_candidate()
    contract[flag] = True
    with pytest.raises(ValueError):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('field,value', [('projection', 'pT'), ('frame', 'thrust'), ('ratio', 'UL'),
    ('source_observable', 'AUC-0-PT'), ('input_units', 'fraction'), ('output_units', 'percent'),
    ('prediction_prescription', 'reported mean P point'), ('source_sha256', '0'*64),
    ('publication', 'unreviewed'), ('overlap_group', ''), ('production_ready', True)])
def test_point_substitution_and_source_mismatch_reject(field, value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0][field] = value
    with pytest.raises(ValueError):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('mutation', ['omit', 'duplicate', 'reorder', 'geometry', 'primary_table', 'primary_column',
                                    'error_column', 'angular_kind', 'Q2', 'z_bin', 'nonfinite'])
def test_literal_source_and_window_provenance_reject_changes(mutation):
    candidate, contract, initial = synthetic_candidate()
    view = candidate['views'][0]
    row = view['points'][0]
    if mutation == 'omit':
        view['points'].pop()
    elif mutation == 'duplicate':
        view['points'][1] = copy.deepcopy(row)
    elif mutation == 'reorder':
        view['points'].reverse()
    elif mutation == 'geometry':
        row['averages']['z1'] += 1e-10
    elif mutation == 'primary_table':
        view['primary_numeric_crosscheck']['table'] = 'tab:totPt'
    elif mutation == 'primary_column':
        view['primary_numeric_crosscheck']['compared_columns'].remove('pT')
    elif mutation == 'error_column':
        row['uncertainty_source_columns']['reported_systematic'] = 'syt_u'
    elif mutation == 'angular_kind':
        row['angular_source']['average_type'] = 'ratio_of_angular_means'
    elif mutation == 'Q2':
        row['averages']['Q2_GeV2'] = 10.58**2
    elif mutation == 'z_bin':
        row['primary_bins']['z1'][1] = row['averages']['z1']-.01
    else:
        row['statistical_error'] = float('nan')
    with pytest.raises(ValueError):
        sia.request_points(candidate, contract, initial)


@pytest.mark.parametrize('bounds', [[.1, 3.5], [0., 3.], [False, 3.5]])
def test_window_is_explicit_and_zero_is_not_a_boolean(bounds):
    _, contract, _ = synthetic_candidate()
    contract['prediction']['qT_window_GeV'] = bounds
    with pytest.raises(ValueError):
        sia.momentum_window(contract)


def test_rehashed_request_window_boolean_is_not_equal_to_float_zero(tmp_path, monkeypatch):
    candidate, contract, initial = synthetic_candidate()
    config = tmp_path/'config'
    config.mkdir()
    for name, value in [('historical-sia-integral-contract.json', contract), ('initial-contract.json', initial),
                        ('sia-qt-window-contract.json', read_json(ROOT/'config/sia-qt-window-contract.json')),
                        ('historical-sources.json', {})]:
        atomic_json(config/name, value)
    source = tmp_path/'candidates'
    source.mkdir()
    atomic_json(source/'candidates.json', candidate)
    atomic_json(source/'report.json', dict(primary_BaBar_raw_points_crosschecked=90,
        primary_BaBar_numeric_cells_crosschecked=954, source_hashes={str(config/'historical-sources.json'): file_hash(config/'historical-sources.json')}))
    monkeypatch.setattr(sia, 'verify_historical_intake', lambda *a: dict(scope='historical_benchmark_candidate', counts={}))
    monkeypatch.setattr(sia, 'verify_primary_upper_cut', lambda *a: {'scope': 'explicitly_synthetic_mock_not_primary_replay'})
    request = tmp_path/'request.json'
    sia.prepare_request(config/'historical-sources.json', tmp_path/'raw', source, request)
    sia.verify_request(request)
    changed = read_json(request)
    changed['momentum_window']['qT_bounds_GeV'][0] = False
    atomic_json(request, changed)
    with pytest.raises(ValueError, match='changed source-derived'):
        sia.verify_request(request)
    with pytest.raises(ValueError, match='fresh'):
        sia.prepare_request(config/'historical-sources.json', tmp_path/'raw', source, request)
