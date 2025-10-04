import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
import json
import os
from model_parameters import CHIKVModelParameters

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class CHIKVTransmissionModel:
    def __init__(self, parameters=None):
        if parameters is None:
            self.params = CHIKVModelParameters()
        else:
            self.params = parameters
        self.p = self.params.get_parameter_dict()
        
    def model_equations(self, y, t):
        S, E, I_s, I_a, C, R = y
        S = max(0, S)
        E = max(0, E)
        I_s = max(0, I_s)
        I_a = max(0, I_a)
        C = max(0, C)
        R = max(0, R)
        N = self.p['N']
        beta = self.p.get('beta', 3.0)
        epsilon = self.p.get('epsilon', 0.375)
        foi = beta * S * (I_s + epsilon * I_a) / N
        dS_dt = self.p['Lambda'] - foi - self.p['mu'] * S
        dE_dt = foi - (self.p['sigma'] + self.p['mu']) * E
        dI_s_dt = self.p['p'] * self.p['sigma'] * E - (self.p['gamma'] + self.p['mu']) * I_s
        dI_a_dt = (1 - self.p['p']) * self.p['sigma'] * E - (self.p['gamma'] + self.p['mu']) * I_a
        dC_dt = self.p['q'] * self.p['gamma'] * I_s - (self.p['omega'] + self.p['mu']) * C
        dR_dt = ((1 - self.p['q']) * self.p['gamma'] * I_s + 
                 self.p['gamma'] * I_a + 
                 self.p['omega'] * C - 
                 self.p['mu'] * R)
        return np.array([dS_dt, dE_dt, dI_s_dt, dI_a_dt, dC_dt, dR_dt])
    
    def set_initial_conditions(self, initial_infected=1):
        N = self.p['N']
        S0 = N - initial_infected
        E0 = 0.0
        I_s0 = initial_infected * self.p['p']
        I_a0 = initial_infected * (1 - self.p['p'])
        C0 = 0.0
        R0 = 0.0
        return np.array([S0, E0, I_s0, I_a0, C0, R0])
    
    def solve_model(self, t_span, initial_conditions=None, **kwargs):
        if initial_conditions is None:
            initial_conditions = self.set_initial_conditions(**kwargs)
        solution = odeint(self.model_equations, initial_conditions, t_span)
        return solution
    
    def extract_compartments(self, solution):
        compartments = {
            'S': solution[:, 0],
            'E': solution[:, 1],
            'I_s': solution[:, 2],
            'I_a': solution[:, 3],
            'C': solution[:, 4],
            'R': solution[:, 5],
        }
        compartments['Total_Infected'] = compartments['I_s'] + compartments['I_a']
        compartments['Total_Human'] = (compartments['S'] + compartments['E'] + 
                                     compartments['I_s'] + compartments['I_a'] + 
                                     compartments['C'] + compartments['R'])
        return compartments
    
    def calculate_daily_incidence(self, solution, t_span):
        E = solution[:, 1]
        daily_incidence = self.p['p'] * self.p['sigma'] * E
        return daily_incidence

def run_foshan_simulation():
    model = CHIKVTransmissionModel()
    foshan_data = pd.read_csv('../fitted_data_foshan.csv')
    onset_data = foshan_data[foshan_data['onset/report'] == 'onset'].copy()
    onset_data['date_parsed'] = pd.to_datetime(onset_data['date'])
    onset_data = onset_data.sort_values('date_parsed')
    start_date = onset_data['date_parsed'].min()
    end_date = onset_data['date_parsed'].max()
    total_days = len(onset_data)
    t_span = np.linspace(0, total_days-1, total_days)
    initial_conditions = model.set_initial_conditions(initial_infected=1)
    solution = model.solve_model(t_span, initial_conditions)
    compartments = model.extract_compartments(solution)
    daily_incidence = model.calculate_daily_incidence(solution, t_span)
    return t_span, solution, compartments, daily_incidence

if __name__ == "__main__":
    t_span, solution, compartments, daily_incidence = run_foshan_simulation()
    results_df = pd.DataFrame({
        'day': t_span,
        'daily_incidence': daily_incidence,
        'cumulative_cases': np.cumsum(daily_incidence),
        'susceptible': compartments['S'],
        'exposed': compartments['E'],
        'infected_symptomatic': compartments['I_s'],
        'infected_asymptomatic': compartments['I_a'],
        'chronic': compartments['C'],
        'recovered': compartments['R']
    })
    save_dir = '../data'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    results_df.to_csv('../data/foshan_model_simulation_results.csv', index=False)
    foshan_data = pd.read_csv('../fitted_data_foshan.csv')
    onset_data = foshan_data[foshan_data['onset/report'] == 'onset'].copy()
    onset_data['date_parsed'] = pd.to_datetime(onset_data['date'])
    onset_data = onset_data.sort_values('date_parsed')
    start_date = onset_data['date_parsed'].min()
    onset_data['days_from_start'] = (onset_data['date_parsed'] - start_date).dt.days
    max_days = min(len(t_span), len(onset_data))
    comparison_df = pd.DataFrame({
        'day': range(max_days),
        'date': onset_data['date'].iloc[:max_days],
        'real_daily_cases': onset_data['data'].iloc[:max_days],
        'model_daily_cases': daily_incidence[:max_days],
        'real_cumulative': np.cumsum(onset_data['data'].iloc[:max_days]),
        'model_cumulative': np.cumsum(daily_incidence[:max_days])
    })
    comparison_df.to_csv('../data/foshan_model_vs_real_comparison.csv', index=False)