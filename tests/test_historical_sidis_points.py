"""Invented source-layout controls, without raw-source replay or physical predictions."""
import copy
from pathlib import Path

import pytest

from collinsflow.historical_sidis_points import FALSE_FLAGS, geometry_key, request_points
from collinsflow.provenance import read_json

ROOT = Path(__file__).resolve().parents[1]


def synthetic_candidate():
    contract = read_json(ROOT/'config/historical-sidis-point-contract.json')
    initial = read_json(ROOT/'config/initial-contract.json')
    views, index = [], 0
    for selected in contract['selection']['source_views']:
        rows = []
        for i in range(selected['points']):
            index += 1
            row = dict(id=f"{selected['id']}:r{i+2}",
                averages=dict(x=.15, y=.3+.001*index, z=.3, Q2_GeV2=4., pT_GeV=.2+index*1e-6),
                central_value=(-1)**i*.0125, statistical_error=.05, systematic_components={},
                locator=dict(sheet='explicitly synthetic layout only', excel_row=i+2), bins={})
            if selected['source_convention_id'] == 'HERMES2010':
                row['systematic_components'] = dict(experimental=.002)
            if selected['source_convention_id'] == 'COMPASS2014':
                row.update(statistical_error=None, reported_error_u=.06,
                           uncertainty_status='reported_error_u_composition_unresolved')
            rows.append(row)
        n_excluded = (3 if selected['source_convention_id'] == 'HERMES2010' else 2) if selected['projection'] == 'pt' else 0
        views.append(dict(id=selected['id'], source_id=selected['id'], process='SIDIS',
            **{k: selected[k] for k in ('publication', 'target', 'hadron', 'projection', 'source_sha256')},
            source_observable='AUTcollins', output_units='fraction', production_ready=False, historical_fit_ready=False,
            depolarization_removed=selected['original_depolarization_removed'],
            blockers=['explicitly synthetic layout; no source replay or admission'], points=rows,
            excluded=[dict(reason='invented excluded-row layout only', excel_row=20+i) for i in range(n_excluded)]))
    candidate = dict(scope='historical_benchmark_candidate', production_ready=False, physical_extraction=False,
                     historical_fit_ready=False, source_fingerprint='explicitly_synthetic_layout', views=views)
    return candidate, contract, initial


def test_complete_exact_means_and_unresolved_errors_are_preserved_without_mutation():
    candidate, contract, initial = synthetic_candidate()
    before = copy.deepcopy(candidate)
    points = request_points(candidate, contract, initial)
    assert len(points) == len(set(map(geometry_key, points))) == 132
    assert candidate == before
    assert [p['id'] for p in points] == [r['id'] for v in candidate['views'] for r in v['points']]
    assert sum(p['measured']['statistical_error'] is None for p in points) == 46
    proton = next(p for p in points if p['publication'] == '1408.4405')
    assert proton['source_row']['reported_error_u'] == .06
    assert proton['model_convention']['phase'] == 'sin_phi_sum_minus_pi'
    assert proton['original_depolarization_removed'] is None
    points[0]['source_row']['systematic_components']['experimental'] = 99.
    assert candidate == before


@pytest.mark.parametrize('flag', FALSE_FLAGS)
def test_no_contract_admission_flag_can_approve_a_request(flag):
    candidate, contract, initial = synthetic_candidate()
    contract[flag] = True
    with pytest.raises(ValueError, match='approve'):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('change', ['scope', 'row', 'geometry', 'view', 'label', 'cut', 'excluded', 'parameter', 'Sigma', 'systematic'])
def test_incomplete_or_changed_source_semantics_are_rejected(change):
    candidate, contract, initial = synthetic_candidate()
    first = candidate['views'][0]
    if change == 'scope':
        candidate['scope'] = 'synthetic_closure'
    elif change == 'row':
        first['points'][0]['id'] = first['points'][1]['id']
    elif change == 'geometry':
        first['points'][1]['averages'] = copy.deepcopy(first['points'][0]['averages'])
    elif change == 'view':
        candidate['views'].pop()
    elif change == 'label':
        first['depolarization_removed'] = True
    elif change == 'cut':
        candidate['views'][2]['points'][0]['averages']['pT_GeV'] = .8
    elif change == 'excluded':
        candidate['views'][2]['excluded'].pop()
    elif change == 'parameter':
        initial['published_parameters'][0] += .01
    elif change == 'Sigma':
        candidate['views'][6]['points'][0]['statistical_error'] = .06
    elif change == 'systematic':
        first['points'][0]['systematic_components']['unapproved_combination'] = .003
    with pytest.raises(ValueError):
        request_points(candidate, contract, initial)


@pytest.mark.parametrize('value', [True, float('nan'), float('inf'), -1., 0., 1.])
def test_source_fraction_cannot_be_boolean_nonfinite_or_outside_domain(value):
    candidate, contract, initial = synthetic_candidate()
    candidate['views'][0]['points'][0]['averages']['x'] = value
    with pytest.raises(ValueError):
        request_points(candidate, contract, initial)
