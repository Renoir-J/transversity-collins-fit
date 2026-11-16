"""Synthetic source-layout tests; these values are not experimental measurements."""
import copy
import math
from pathlib import Path

import pytest

from collinsflow.historical_points import request_points
from collinsflow.provenance import read_json

ROOT = Path(__file__).resolve().parents[1]


def synthetic_candidate():
    contract = read_json(ROOT/'config/historical-point-contract.json')
    initial = read_json(ROOT/'config/initial-contract.json')
    views = []
    for selected in contract['selection']['source_views']:
        points = []
        for i in range(4):
            points.append(dict(id=f"{selected['id']}:r{i+2}",
                averages=dict(x=.1+.02*i, y=.4+.01*i, z=.3+.02*i, Q2_GeV2=2.+i, pT_GeV=.2),
                central_value=(-1)**i*.0125, statistical_error=.05,
                systematic_components=dict(experimental=.002, model=.003),
                locator=dict(sheet='synthetic-source-only', excel_row=i+2),
                source_remap=dict(kind='reviewed_legacy_column_remap',
                    primary_tables=['table:kine', 'table:results'], primary_lines=[10+i, 30+i])))
        views.append(dict(id=selected['id'], source_id=selected['id'], process='SIDIS', projection='x',
            publication=contract['selection']['publication'], target='neutron', target_material='3He',
            target_result=contract['selection']['target_result'], hadron=selected['hadron'],
            source_sha256=selected['sha256'], source_observable='AUTcollins', output_units='fraction',
            depolarization_removed=False, production_ready=False, historical_fit_ready=False,
            blockers=['synthetic layout fixture; no primary-source replay or experimental admission'],
            points=points, excluded=[]))
    candidate = dict(scope='historical_benchmark_candidate', production_ready=False,
                     physical_extraction=False, historical_fit_ready=False, views=views)
    return candidate, contract, initial


def test_representation_preserves_configured_order_means_signs_and_separate_errors():
    candidate, contract, initial = synthetic_candidate()
    before = copy.deepcopy(candidate)
    points = request_points(candidate, contract, initial)
    assert len(points) == 8 and candidate == before
    assert [p['id'] for p in points] == [p['id'] for v in candidate['views'] for p in v['points']]
    assert [p['charge'] for p in points] == [1]*4+[-1]*4
    assert points[1]['measured']['central_value'] == -.0125
    assert points[0]['measured']['systematic_components'] == dict(experimental=.002, model=.003)
    assert all(p['original_depolarization_removed'] is False for p in points)
    assert all('synthetic layout fixture' in p['source_blockers'][0] for p in points)
    points[0]['measured']['systematic_components']['model'] = 100.
    assert candidate == before


@pytest.mark.parametrize('scope', ['production', 'synthetic_closure', 'experimental_candidate'])
def test_other_scope_cannot_enter_the_historical_request(scope):
    candidate, contract, initial = synthetic_candidate()
    candidate['scope'] = scope
    with pytest.raises(ValueError, match='scope'):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('flag', ['production_ready', 'physical_extraction', 'historical_fit_ready',
                                  'experimental_data_admitted', 'fit_performed'])
def test_diagnostic_contract_cannot_self_approve_any_fit_or_production_flag(flag):
    candidate, contract, initial = synthetic_candidate()
    contract[flag] = True
    with pytest.raises(ValueError, match='approve'):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('field,value', [
    ('target', '3He'), ('hadron', 'pi-'), ('publication', 'other-version'),
    ('source_sha256', '0'*64), ('depolarization_removed', None),
    ('depolarization_removed', True), ('source_observable', 'AUTsivers'), ('output_units', 'percent')])
def test_primary_identity_and_preserved_source_label_are_required(field, value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0][field] = value
    with pytest.raises(ValueError, match='selection'):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('changed', [None, True])
def test_contract_cannot_erase_or_promote_the_existing_historical_source_label(changed):
    candidate, contract, initial = synthetic_candidate()
    contract['selection']['source_views'][0]['depolarization_removed'] = changed
    with pytest.raises(ValueError, match='selection'):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('mutation', ['omit_view', 'omit_point', 'duplicate_point', 'reorder_points', 'wrong_source_row'])
def test_requested_points_cannot_be_omitted_reordered_or_relabelled(mutation):
    candidate, contract, initial = synthetic_candidate()
    view = candidate['views'][0]
    if mutation == 'omit_view':
        candidate['views'].pop()
    elif mutation == 'omit_point':
        view['points'].pop()
    elif mutation == 'duplicate_point':
        view['points'][1] = copy.deepcopy(view['points'][0])
    elif mutation == 'reorder_points':
        view['points'].reverse()
    else:
        view['points'][0]['locator']['excel_row'] = 100
    with pytest.raises(ValueError):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('field,value', [('x', True), ('x', 0.), ('y', 1.), ('z', math.nan),
                                       ('Q2_GeV2', -1.), ('pT_GeV', -.1)])
def test_invalid_kinematic_inputs_are_not_silently_cast_or_replaced(field, value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0]['points'][0]['averages'][field] = value
    with pytest.raises(ValueError):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('mutation', ['missing_model', 'negative_model', 'missing_stat', 'zero_stat', 'missing_remap'])
def test_missing_or_changed_uncertainty_and_primary_remap_are_explicit_failures(mutation):
    candidate, contract, initial = synthetic_candidate()
    point = candidate['views'][0]['points'][0]
    if mutation == 'missing_model':
        point['systematic_components'].pop('model')
    elif mutation == 'negative_model':
        point['systematic_components']['model'] = -.001
    elif mutation == 'missing_stat':
        point['statistical_error'] = None
    elif mutation == 'zero_stat':
        point['statistical_error'] = 0.
    else:
        point.pop('source_remap')
    with pytest.raises(ValueError):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('mutation', ['wrong_phase', 'wrong_normalization', 'approve_mapping', 'missing_parameter',
                                     'duplicate_parameter', 'nonfinite_parameter'])
def test_complete_fixed_model_definition_is_required(mutation):
    candidate, contract, initial = synthetic_candidate()
    if mutation == 'wrong_phase':
        contract['prediction']['phase'] = 'sin_phi_sum_minus_pi'
    elif mutation == 'wrong_normalization':
        contract['prediction']['normalization'] = 'DNN_removed'
    elif mutation == 'approve_mapping':
        contract['prediction']['source_mapping_status'] = 'approved'
    elif mutation == 'missing_parameter':
        initial['published_parameters'].pop()
    elif mutation == 'duplicate_parameter':
        initial['parameter_order'][-1] = initial['parameter_order'][0]
    else:
        initial['published_parameters'][0] = math.inf
    with pytest.raises(ValueError):
        request_points(candidate, contract, initial)
