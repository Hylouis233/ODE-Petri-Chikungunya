import os
import json
import numpy as np
import pandas as pd
from typing import Tuple
from scipy.optimize import differential_evolution
from model_parameters import CHIKVModelParameters
from petri import CHIKVPetriNetModel
from petri_r_value_calculator import PetriNetRCalculator

def load_foshan_data_with_phases():
    data_path = os.path.join(os.path.dirname(__file__), '..', 'fitted_data_foshan.csv')
    data_path = os.path.normpath(data_path)
    data = pd.read_csv(data_path)
    onset_data = data[data['onset/report'] == 'onset'].copy()
    onset_data['date_parsed'] = pd.to_datetime(onset_data['date'])
    onset_data = onset_data.sort_values(by='date_parsed')
    start_date = pd.to_datetime('2025-06-16')
    onset_data['days_from_start'] = (onset_data['date_parsed'] - start_date).dt.days
    daily_cases = pd.to_numeric(onset_data['data'], errors='coerce').fillna(0).values
    days = onset_data['days_from_start'].values
    intervention_dates = {
        'reff1_start': pd.to_datetime('2025-07-09'),
        'reff2_start': pd.to_datetime('2025-07-14'),
        'reff3_start': pd.to_datetime('2025-07-18')
    }
    phase_days = {k: (v - start_date).days for k, v in intervention_dates.items()}
    return days, daily_cases, phase_days

def petri_model_with_three_phase_intervention(t, initial_conditions, beta_0, phase_days,
                                              beta1, beta2, beta3):
    t1 = phase_days['reff1_start']
    t2 = phase_days['reff2_start']
    t3 = phase_days['reff3_start']
    max_t = len(t) - 1
    solutions = []
    current_conditions = initial_conditions.copy()
    state_names = ['S', 'E', 'I_s', 'I_a', 'C', 'R']
    if t1 > 0:
        t_phase0 = t[t <= t1]
        if len(t_phase0) == 0:
            t_phase0 = np.array([0])
        model_phase0 = CHIKVPetriNetModel()
        model_phase0.set_initial_conditions()
        for i, name in enumerate(state_names):
            model_phase0.places[name].tokens = current_conditions[i]
            model_phase0.places[name].initial_tokens = current_conditions[i]
        sol_phase0 = model_phase0.simulate(t_phase0, dt=0.01, beta=beta_0)
        solutions.append(sol_phase0)
        current_conditions = sol_phase0[-1, :].copy()
    if t2 > t1 and t1 < max_t:
        t_phase1 = t[(t > t1) & (t <= t2)]
        if len(t_phase1) > 0:
            t_phase1_adj = np.concatenate([[0], t_phase1 - t1])
            model_phase1 = CHIKVPetriNetModel()
            model_phase1.set_initial_conditions()
            for i, name in enumerate(state_names):
                model_phase1.places[name].tokens = current_conditions[i]
                model_phase1.places[name].initial_tokens = current_conditions[i]
            sol_phase1 = model_phase1.simulate(t_phase1_adj, dt=0.01, beta=beta1)
            solutions.append(sol_phase1[1:, :])
            current_conditions = sol_phase1[-1, :].copy()
    if t3 > t2 and t2 < max_t:
        t_phase2 = t[(t > t2) & (t <= t3)]
        if len(t_phase2) > 0:
            t_phase2_adj = np.concatenate([[0], t_phase2 - t2])
            model_phase2 = CHIKVPetriNetModel()
            model_phase2.set_initial_conditions()
            for i, name in enumerate(state_names):
                model_phase2.places[name].tokens = current_conditions[i]
                model_phase2.places[name].initial_tokens = current_conditions[i]
            sol_phase2 = model_phase2.simulate(t_phase2_adj, dt=0.01, beta=beta2)
            solutions.append(sol_phase2[1:, :])
            current_conditions = sol_phase2[-1, :].copy()
    if t3 < max_t:
        t_phase3 = t[t > t3]
        if len(t_phase3) > 0:
            t_phase3_adj = np.concatenate([[0], t_phase3 - t3])
            model_phase3 = CHIKVPetriNetModel()
            model_phase3.set_initial_conditions()
            for i, name in enumerate(state_names):
                model_phase3.places[name].tokens = current_conditions[i]
                model_phase3.places[name].initial_tokens = current_conditions[i]
            sol_phase3 = model_phase3.simulate(t_phase3_adj, dt=0.01, beta=beta3)
            solutions.append(sol_phase3[1:, :])
    solution = np.vstack(solutions) if solutions else np.array([initial_conditions])
    if len(solution) < len(t):
        last_state = solution[-1, :] if len(solution) > 0 else initial_conditions
        padding = np.tile(last_state, (len(t) - len(solution), 1))
        solution = np.vstack([solution, padding])
    elif len(solution) > len(t):
        solution = solution[:len(t), :]
    E = solution[:, 1]
    model_params = CHIKVModelParameters()
    sigma = model_params.sigma
    p_val = model_params.p
    daily_incidence = p_val * sigma * E
    return daily_incidence, solution

def objective_function_three_phases_petri(params_to_fit, days_observed, cases_observed, phase_days):
    beta_0, beta1, beta2, beta3 = params_to_fit
    if (beta_0 <= 0 or beta_0 > 100.0 or
        beta1 <= 0 or beta1 > 100.0 or
        beta2 <= 0 or beta2 > 100.0 or
        beta3 <= 0 or beta3 > 100.0 or
        beta1 >= beta_0 or beta2 >= beta1 or beta3 >= beta2):
        return 1e10
    max_day = max(days_observed) + 14
    t_span = np.linspace(0, max_day, max_day + 1)
    model = CHIKVPetriNetModel(CHIKVModelParameters())
    model.set_initial_conditions(initial_infected=1)
    initial_conditions = model.get_state_vector()
    daily_incidence, _ = petri_model_with_three_phase_intervention(
        t_span, initial_conditions, beta_0, phase_days, beta1, beta2, beta3
    )
    predicted_cases = np.interp(days_observed, t_span, daily_incidence)
    rmse = np.sqrt(np.mean((predicted_cases - cases_observed) ** 2))
    return rmse

def fit_three_phase_petri_model():
    days_observed, cases_observed, phase_days = load_foshan_data_with_phases()
    bounds = [(0.01, 15.0)] * 4
    result = differential_evolution(
        objective_function_three_phases_petri,
        bounds,
        args=(days_observed, cases_observed, phase_days),
        seed=42,
        maxiter=400,
        popsize=100,
        tol=1e-12,
        workers=-1,
    )
    beta_0, beta1, beta2, beta3 = result.x
    r_calculator = PetriNetRCalculator(CHIKVPetriNetModel(CHIKVModelParameters()))
    r_results_with_ci = r_calculator.calculate_phase_r_values(beta_0, beta1, beta2, beta3)
    phase_names = {
        'phase_0_initial': 'Initial Phase (No Control)',
        'phase_1_intervention': 'First Phase Control',
        'phase_2_intervention': 'Second Phase Control',
        'phase_3_intervention': 'Third Phase Control'
    }
    for phase_key, desc in phase_names.items():
        pr = r_results_with_ci['petri_net_results_with_ci'][phase_key]
        ci_ngm = pr['confidence_intervals']['next_generation_matrix']
    return result.x, days_observed, cases_observed, phase_days, r_results_with_ci

def validate_three_phase_petri_model(fitted_params, days_observed, cases_observed, phase_days):
    beta_0, beta1, beta2, beta3 = fitted_params
    t_span = np.linspace(0, 70, 71)
    model = CHIKVPetriNetModel()
    model.set_initial_conditions(initial_infected=1)
    initial_conditions = model.get_state_vector()
    daily_incidence, solution = petri_model_with_three_phase_intervention(
        t_span, initial_conditions, beta_0, phase_days, beta1, beta2, beta3
    )
    predicted_obs = np.interp(days_observed, t_span, daily_incidence)
    rmse = np.sqrt(np.mean((predicted_obs - cases_observed) ** 2))
    mae = np.mean(np.abs(predicted_obs - cases_observed))
    ss_res = np.sum((cases_observed - predicted_obs) ** 2)
    ss_tot = np.sum((cases_observed - np.mean(cases_observed)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    return t_span, daily_incidence, solution

def save_petri_net_results(fitted_params, phase_days, validation_metrics, r_results):
    beta_0, beta1, beta2, beta3 = fitted_params
    mp = CHIKVModelParameters()
    p = mp.p
    gamma = mp.gamma
    epsilon = mp.epsilon
    effective_infectious_period = p / gamma + (1 - p) * epsilon / gamma
    R0_initial = beta_0 * effective_infectious_period
    results = {
        'model_type': 'petri_net_three_phase_intervention_CHIKV',
        'fitting_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'fitted_parameters': {
            'beta_0': float(beta_0),
            'initial_infected': 1,
            'beta1': float(beta1),
            'beta2': float(beta2),
            'beta3': float(beta3)
        },
        'epidemiological_parameters': {
            'R0_initial': float(R0_initial),
            'R_eff1': float(beta1 * effective_infectious_period),
            'R_eff2': float(beta2 * effective_infectious_period),
            'R_eff3': float(beta3 * effective_infectious_period),
            'effective_infectious_period': float(effective_infectious_period)
        },
        'intervention_timeline': {
            'phase1_start_day': int(phase_days['reff1_start']),
            'phase2_start_day': int(phase_days['reff2_start']),
            'phase3_start_day': int(phase_days['reff3_start']),
            'phase1_date': '2025-07-09',
            'phase2_date': '2025-07-14',
            'phase3_date': '2025-07-18'
        },
        'control_effectiveness': {
            'phase1_reduction_percent': float((1 - beta1 / beta_0) * 100),
            'phase2_reduction_percent': float((1 - beta2 / beta_0) * 100),
            'phase3_reduction_percent': float((1 - beta3 / beta_0) * 100)
        },
        'validation_metrics': validation_metrics,
        'r_value_results': r_results
    }
    results_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'petri_net_fitting_results.json')
    results_path = os.path.normpath(results_path)
    with open(results_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

def save_plotting_data(t_span, daily_incidence, days_obs, cases_obs, fitted_params, phase_days, solution):
    plotting_data = {
        'model_type': 'petri_net_three_phase_intervention_CHIKV',
        'saved_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'time_series': {
            't_span': t_span.tolist(),
            'daily_incidence': daily_incidence.tolist(),
            'days_observed': days_obs.tolist(),
            'cases_observed': cases_obs.tolist()
        },
        'fitted_parameters': {
            'beta_0': float(fitted_params[0]),
            'beta1': float(fitted_params[1]),
            'beta2': float(fitted_params[2]),
            'beta3': float(fitted_params[3])
        },
        'phase_days': {
            'reff1_start': int(phase_days['reff1_start']),
            'reff2_start': int(phase_days['reff2_start']),
            'reff3_start': int(phase_days['reff3_start'])
        },
        'model_parameters': {
            'p': float(CHIKVModelParameters().p),
            'gamma': float(CHIKVModelParameters().gamma),
            'epsilon': float(CHIKVModelParameters().epsilon)
        },
        'data_summary': {
            'total_predicted_cases': float(np.sum(daily_incidence)),
            'total_observed_cases': float(np.sum(cases_obs)),
            'simulation_days': len(t_span),
            'observation_days': len(days_obs)
        }
    }
    plotting_data_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'petri_net_plotting_data.json')
    plotting_data_path = os.path.normpath(plotting_data_path)
    with open(plotting_data_path, 'w', encoding='utf-8') as f:
        json.dump(plotting_data, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    fitted_params, days_obs, cases_obs, phase_days, r_results = fit_three_phase_petri_model()
    t_span, daily_incidence, solution = validate_three_phase_petri_model(
        fitted_params, days_obs, cases_obs, phase_days)
    predicted_at_obs = np.interp(np.array(days_obs, dtype=float),
                                 np.array(t_span, dtype=float),
                                 np.array(daily_incidence, dtype=float))
    cases_obs_numeric = np.array(cases_obs, dtype=float)
    validation_metrics = {
        'rmse': float(np.sqrt(np.mean((predicted_at_obs - cases_obs_numeric) ** 2))),
        'mae': float(np.mean(np.abs(predicted_at_obs - cases_obs_numeric))),
        'total_predicted': float(np.sum(daily_incidence)),
        'total_observed': float(np.sum(cases_obs_numeric)),
        'relative_error_percent': float((np.sum(daily_incidence) - np.sum(cases_obs_numeric)) / np.sum(cases_obs_numeric) * 100)
    }
    save_petri_net_results(fitted_params, phase_days, validation_metrics, r_results)
    save_plotting_data(t_span, daily_incidence, days_obs, cases_obs, fitted_params, phase_days, solution)