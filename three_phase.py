import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import os

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

def load_plotting_data():
    data_dir = "data"
    
    ode_file = os.path.join(data_dir, "ode_three_phase_plotting_data.json")
    with open(ode_file, 'r', encoding='utf-8') as f:
        ode_data = json.load(f)
    
    petri_file = os.path.join(data_dir, "petri_net_plotting_data.json")
    with open(petri_file, 'r', encoding='utf-8') as f:
        petri_data = json.load(f)
    
    return ode_data, petri_data

def plot_combined_three_phase_results(ode_data, petri_data, save_path=None):
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    plot_model_results(ode_data, axes[0], "ODE Model", model_type="ode")
    
    plot_model_results(petri_data, axes[1], "Petri Net Model", model_type="petri")
    
    labels = ['a', 'b', 'c', 'd', 'e', 'f']
    for i, ax in enumerate(axes.flat):
        ax.text(0.05, 0.95, labels[i], transform=ax.transAxes, 
                fontsize=14, fontweight='bold', va='top', ha='left')
    
    plt.tight_layout()
    
    if save_path:
        png_path = save_path
        pdf_path = save_path.replace('.png', '.pdf')
        
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, bbox_inches='tight')
    
    return fig

def plot_model_results(data, axes_row, model_name, model_type="ode"):
    time_series = data['time_series']
    fitted_params = data['fitted_parameters']
    phase_days = data['phase_days']
    
    t_span = time_series['t_span']
    start_date = datetime(2024, 6, 30)
    dates = [start_date + timedelta(days=int(t)) for t in t_span]
    
    intervention_dates = [
        start_date + timedelta(days=phase_days['reff1_start']),
        start_date + timedelta(days=phase_days['reff2_start']),
        start_date + timedelta(days=phase_days['reff3_start'])
    ]
    
    phase_colors = ['lightcoral', 'lightblue', 'lightgreen', 'lightyellow']
    phase_dates = [dates[0]] + intervention_dates + [dates[-1]]
    
    ax1 = axes_row[0]
    observed_days = time_series['days_observed']
    observed_cases = time_series['cases_observed']
    
    observed_dates = [start_date + timedelta(days=int(d)) for d in observed_days]
    
    ax1.plot(dates, time_series['daily_incidence'], 'b-', linewidth=2, label='Fitted curve')
    ax1.scatter(observed_dates, observed_cases, color='red', s=40, alpha=0.8, label='Observed data', zorder=5)
    
    for i, date in enumerate(intervention_dates):
        ax1.axvline(x=date, color='gray', linestyle='--', alpha=0.7)
        ax1.text(date, max(time_series['daily_incidence']) * 0.8, f'Phase {i+1}', 
                rotation=90, ha='right', va='bottom')
    
    for i in range(4):
        ax1.axvspan(phase_dates[i], phase_dates[i+1], alpha=0.2, color=phase_colors[i])
    
    ax1.set_xlabel('Date')
    if model_type == "ode":
        ax1.set_ylabel('Daily incidence')
        ax1.set_title(f'{model_name}: Daily Incidence Fitting')
    else:
        ax1.set_ylabel('Daily incidence rate')
        ax1.set_title(f'{model_name}: Daily Incidence Rate Fitting')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax1.xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45)
    
    ax2 = axes_row[1]
    cumulative_observed = np.cumsum(observed_cases)
    cumulative_fitted = np.cumsum(time_series['daily_incidence'])
    
    ax2.plot(dates, cumulative_fitted, 'b-', linewidth=2, label='Fitted curve')
    ax2.plot(observed_dates, cumulative_observed, 'ro-', markersize=4, linewidth=2, label='Observed data')
    
    for i, date in enumerate(intervention_dates):
        ax2.axvline(x=date, color='gray', linestyle='--', alpha=0.7)
        ax2.text(date, max(cumulative_fitted) * 0.8, f'Phase {i+1}', 
                rotation=90, ha='right', va='bottom')
    
    for i in range(4):
        ax2.axvspan(phase_dates[i], phase_dates[i+1], alpha=0.2, color=phase_colors[i])
    
    ax2.set_xlabel('Date')
    ax2.set_ylabel('Cumulative cases')
    ax2.set_title(f'{model_name}: Cumulative Cases Comparison')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax2.xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)
    
    ax3 = axes_row[2]
    
    model_params = data.get('model_parameters', {'p': 0.75, 'gamma': 0.143, 'epsilon': 0.375})
    p = model_params['p']
    gamma = model_params['gamma']
    epsilon = model_params['epsilon']
    effective_infectious_period = p/gamma + (1-p)*epsilon/gamma
    
    beta_0 = fitted_params['beta_0']
    beta1 = fitted_params['beta1']
    beta2 = fitted_params['beta2']
    beta3 = fitted_params['beta3']
    
    R0_initial = beta_0 * effective_infectious_period
    R_eff1 = beta1 * effective_infectious_period
    R_eff2 = beta2 * effective_infectious_period
    R_eff3 = beta3 * effective_infectious_period
    
    rt_values = []
    for t in t_span:
        if t <= phase_days['reff1_start']:
            rt_values.append(R0_initial)
        elif t <= phase_days['reff2_start']:
            rt_values.append(R_eff1)
        elif t <= phase_days['reff3_start']:
            rt_values.append(R_eff2)
        else:
            rt_values.append(R_eff3)
    
    ax3.plot(dates, rt_values, 'g-', linewidth=2, label='R')
    ax3.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='R = 1')
    
    for i, date in enumerate(intervention_dates):
        ax3.axvline(x=date, color='gray', linestyle='--', alpha=0.7)
        ax3.text(date, max(rt_values) * 0.8, f'Phase {i+1}', 
                rotation=90, ha='right', va='bottom')
    
    for i in range(4):
        ax3.axvspan(phase_dates[i], phase_dates[i+1], alpha=0.2, color=phase_colors[i])
    
    ax3.set_xlabel('Date')
    ax3.set_ylabel('Effective reproduction number R')
    ax3.set_title(f'{model_name}: R Time Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax3.xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45)

def plot_both_models():
    ode_data, petri_data = load_plotting_data()
    
    os.makedirs("figure", exist_ok=True)
    
    save_path = "figure/combined_three_phase_fitting.png"
    fig = plot_combined_three_phase_results(ode_data, petri_data, save_path)
    
    plt.show()
    return fig

if __name__ == "__main__":
    plot_both_models()