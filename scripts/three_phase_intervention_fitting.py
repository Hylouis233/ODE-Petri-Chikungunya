import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
import pandas as pd
import sys
import os
import multiprocessing as mp
import json

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from ODE import CHIKVTransmissionModel
from model_parameters import CHIKVModelParameters
from r_value_calculator import RValueCalculator

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_foshan_data_with_phases():
    data_path = os.path.join(os.path.dirname(current_dir), 'fitted_data_foshan.csv')
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
    
    phase_days = {}
    for phase, date in intervention_dates.items():
        phase_days[phase] = (date - start_date).days
    
    return days, daily_cases, phase_days

def model_with_three_phase_intervention(t, initial_conditions, beta_0, phase_days, beta1, beta2, beta3):
    model_params = CHIKVModelParameters()
    
    t1 = phase_days['reff1_start']
    t2 = phase_days['reff2_start']
    t3 = phase_days['reff3_start']
    
    max_t = len(t) - 1
    
    solutions = []
    current_conditions = initial_conditions.copy()
    
    if t1 > 0:
        t_phase0 = t[t <= t1]
        if len(t_phase0) == 0:
            t_phase0 = np.array([0])
        
        model_params.beta = beta_0
        model_phase0 = CHIKVTransmissionModel(model_params)
        sol_phase0 = model_phase0.solve_model(t_phase0, current_conditions)
        solutions.append(sol_phase0)
        current_conditions = sol_phase0[-1, :].copy()
    
    if t2 > t1 and t1 < max_t:
        t_phase1 = t[(t > t1) & (t <= t2)]
        if len(t_phase1) > 0:
            t_phase1_adj = np.concatenate([[0], t_phase1 - t1])
            
            model_params.beta = beta1
            model_phase1 = CHIKVTransmissionModel(model_params)
            sol_phase1 = model_phase1.solve_model(t_phase1_adj, current_conditions)
            solutions.append(sol_phase1[1:, :])
            current_conditions = sol_phase1[-1, :].copy()
    
    if t3 > t2 and t2 < max_t:
        t_phase2 = t[(t > t2) & (t <= t3)]
        if len(t_phase2) > 0:
            t_phase2_adj = np.concatenate([[0], t_phase2 - t2])
            
            model_params.beta = beta2
            model_phase2 = CHIKVTransmissionModel(model_params)
            sol_phase2 = model_phase2.solve_model(t_phase2_adj, current_conditions)
            solutions.append(sol_phase2[1:, :])
            current_conditions = sol_phase2[-1, :].copy()
    
    if t3 < max_t:
        t_phase3 = t[t > t3]
        if len(t_phase3) > 0:
            t_phase3_adj = np.concatenate([[0], t_phase3 - t3])
            
            model_params.beta = beta3
            model_phase3 = CHIKVTransmissionModel(model_params)
            sol_phase3 = model_phase3.solve_model(t_phase3_adj, current_conditions)
            solutions.append(sol_phase3[1:, :])
    
    if solutions:
        solution = np.vstack(solutions)
    else:
        solution = np.array([initial_conditions])
    
    if len(solution) < len(t):
        last_state = solution[-1, :] if len(solution) > 0 else initial_conditions
        padding = np.tile(last_state, (len(t) - len(solution), 1))
        solution = np.vstack([solution, padding])
    elif len(solution) > len(t):
        solution = solution[:len(t), :]
    
    E = solution[:, 1]
    sigma = model_params.sigma
    p_val = model_params.p
    daily_incidence = p_val * sigma * E
    
    return daily_incidence, solution

def objective_function_three_phases(params_to_fit, days_observed, cases_observed, phase_days):
    beta_0, beta1, beta2, beta3 = params_to_fit
    
    if (beta_0 <= 0 or beta_0 > 100.0 or 
        beta1 <= 0 or beta1 > 100.0 or
        beta2 <= 0 or beta2 > 100.0 or 
        beta3 <= 0 or beta3 > 100.0 or
        beta1 >= beta_0 or beta2 >= beta1 or beta3 >= beta2):
        return 1e10
    
    max_day = max(days_observed) + 14
    t_span = np.linspace(0, max_day, max_day + 1)
    
    model_params = CHIKVModelParameters()
    model_params.beta = beta_0
    model = CHIKVTransmissionModel(model_params)
    initial_conditions = model.set_initial_conditions(initial_infected=1)
    
    daily_incidence, _ = model_with_three_phase_intervention(
        t_span, initial_conditions, beta_0, 
        phase_days, beta1, beta2, beta3
    )
    
    predicted_cases = np.interp(days_observed, t_span, daily_incidence)
    
    rmse = np.sqrt(np.mean((predicted_cases - cases_observed) ** 2))
    
    return rmse

def fit_three_phase_model():
    days_observed, cases_observed, phase_days = load_foshan_data_with_phases()
    
    bounds = [
        (0.01, 15.0),
        (0.01, 15.0),
        (0.01, 15.0),
        (0.01, 15.0),
    ]
    
    max_workers = mp.cpu_count()
    
    result = differential_evolution(
        objective_function_three_phases,
        bounds,
        args=(days_observed, cases_observed, phase_days),
        seed=42,
        maxiter=400,
        popsize=100,
        tol=1e-12,
        workers=max_workers
    )
    
    beta_0, beta1, beta2, beta3 = result.x
    
    model_params = CHIKVModelParameters()
    r_calculator = RValueCalculator(model_params)
    r_results_with_ci = r_calculator.calculate_phase_r_values(
        beta_0, beta1, beta2, beta3, confidence_level=0.95, n_samples=1000
    )
    
    return result.x, days_observed, cases_observed, phase_days, r_results_with_ci

def validate_three_phase_model(fitted_params, days_observed, cases_observed, phase_days):
    beta_0, beta1, beta2, beta3 = fitted_params
    
    t_span = np.linspace(0, 70, 71)
    
    model_params = CHIKVModelParameters()
    model_params.beta = beta_0
    model = CHIKVTransmissionModel(model_params)
    initial_conditions = model.set_initial_conditions(initial_infected=1)
    
    daily_incidence, solution = model_with_three_phase_intervention(
        t_span, initial_conditions, beta_0,
        phase_days, beta1, beta2, beta3
    )
    
    predicted_obs = np.interp(days_observed, t_span, daily_incidence)
    
    rmse = np.sqrt(np.mean((predicted_obs - cases_observed) ** 2))
    mae = np.mean(np.abs(predicted_obs - cases_observed))
    
    ss_res = np.sum((cases_observed - predicted_obs) ** 2)
    ss_tot = np.sum((cases_observed - np.mean(cases_observed)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    return t_span, daily_incidence, solution


def save_three_phase_results(fitted_params, phase_days, validation_metrics):
    
    beta_0, beta1, beta2, beta3 = fitted_params
    model_params = CHIKVModelParameters()
    p = model_params.p
    gamma = model_params.gamma
    epsilon = model_params.epsilon
    effective_infectious_period = p/gamma + (1-p)*epsilon/gamma
    R0_initial = beta_0 * effective_infectious_period
    
    results = {
        'model_type': 'three_phase_intervention_CHIKV',
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
            'phase1_reduction_percent': float((1-beta1/beta_0)*100),
            'phase2_reduction_percent': float((1-beta2/beta_0)*100),
            'phase3_reduction_percent': float((1-beta3/beta_0)*100)
        },
        'validation_metrics': validation_metrics
    }
    
    results_path = os.path.join(os.path.dirname(current_dir), 'data', 'three_phase_fitting_results.json')
    with open(results_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, ensure_ascii=False, indent=2)


def save_ode_plotting_data(t_span, daily_incidence, days_obs, cases_obs, fitted_params, phase_days, solution):
    
    plotting_data = {
        'model_type': 'three_phase_intervention_CHIKV_ODE',
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
    
    plotting_data_path = os.path.join(os.path.dirname(current_dir), 'data', 'ode_three_phase_plotting_data.json')
    
    with open(plotting_data_path, 'w', encoding='utf-8') as f:
        json.dump(plotting_data, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    fitted_params, days_obs, cases_obs, phase_days, r_results = fit_three_phase_model()
    
    t_span, daily_incidence, solution = validate_three_phase_model(fitted_params, days_obs, cases_obs, phase_days)
    
    predicted_at_obs = np.interp(np.array(days_obs, dtype=float), 
                                 np.array(t_span, dtype=float), 
                                 np.array(daily_incidence, dtype=float))
    cases_obs_numeric = np.array(cases_obs, dtype=float)
    
    validation_metrics = {
        'rmse': float(np.sqrt(np.mean((predicted_at_obs - cases_obs_numeric) ** 2))),
        'total_predicted': float(np.sum(daily_incidence)),
        'total_observed': float(np.sum(cases_obs_numeric))
    }
    save_three_phase_results(fitted_params, phase_days, validation_metrics)
    
    save_ode_plotting_data(t_span, daily_incidence, days_obs, cases_obs, fitted_params, phase_days, solution)