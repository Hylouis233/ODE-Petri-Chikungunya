import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class ChikungunyaPlotter:
    
    def __init__(self, data_dir='data', output_dir='figures'):
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.load_all_data()
        
    def load_all_data(self):
        with open(f'{self.data_dir}/ode_three_phase_plotting_data.json', 'r', encoding='utf-8') as f:
            self.ode_data = json.load(f)
        
        with open(f'{self.data_dir}/petri_net_plotting_data.json', 'r', encoding='utf-8') as f:
            self.petri_data = json.load(f)
            
        with open(f'{self.data_dir}/three_phase_comparison_plotting_data.json', 'r', encoding='utf-8') as f:
            self.comparison_data = json.load(f)
            
        with open(f'{self.data_dir}/three_phase_comparison_results.json', 'r', encoding='utf-8') as f:
            self.comparison_results = json.load(f)
            
        with open(f'{self.data_dir}/petri_net_fitting_results.json', 'r', encoding='utf-8') as f:
            self.petri_results = json.load(f)
            
        with open(f'{self.data_dir}/three_phase_fitting_results.json', 'r', encoding='utf-8') as f:
            self.three_phase_results = json.load(f)
            
        with open(f'{self.data_dir}/sensitivity_analysis_plotting_data.json', 'r', encoding='utf-8') as f:
            self.sensitivity_data = json.load(f)
        with open(f'{self.data_dir}/sensitivity_analysis_results.json', 'r', encoding='utf-8') as f:
            self.sensitivity_results = json.load(f)

    def generate_all_plots(self):
        self.plot_model_validation()
    
    def plot_model_validation(self):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Chikungunya Model Validation and Diagnostic Analysis', fontsize=16, fontweight='bold')
        
        ode_time = np.array(self.ode_data['time_series']['t_span'])
        ode_daily = np.array(self.ode_data['time_series']['daily_incidence'])
        observed_cases = np.array(self.ode_data['time_series']['cases_observed'])
        
        petri_time = np.array(self.petri_data['time_series']['t_span'])
        petri_daily = np.array(self.petri_data['time_series']['daily_incidence'])
        
        min_len = min(len(ode_time), len(ode_daily), len(observed_cases), 
                     len(petri_time), len(petri_daily))
        
        ode_time = ode_time[:min_len]
        ode_daily = ode_daily[:min_len]
        petri_time = petri_time[:min_len]
        petri_daily = petri_daily[:min_len]
        observed_cases = observed_cases[:min_len]
        
        ax1 = axes[0, 0]
        ax1.text(0.05, 0.95, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        ode_residuals = observed_cases - ode_daily
        ax1.scatter(ode_time, ode_residuals, alpha=0.6, color='blue', s=30)
        ax1.axhline(y=0, color='red', linestyle='--', alpha=0.7)
        ax1.set_xlabel('Time (Days)')
        ax1.set_ylabel('Residuals (Observed - Predicted)')
        ax1.set_title('ODE Model Residual Analysis')
        ax1.grid(True, alpha=0.3)
        
        mean_residual = np.mean(ode_residuals)
        std_residual = np.std(ode_residuals)
        ax1.text(0.05, 0.85, f'Mean: {mean_residual:.3f}\nStd: {std_residual:.3f}', 
                transform=ax1.transAxes, verticalalignment='top')
        
        ax2 = axes[0, 1]
        ax2.text(0.05, 0.95, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        petri_residuals = observed_cases - petri_daily
        ax2.scatter(petri_time, petri_residuals, alpha=0.6, color='green', s=30)
        ax2.axhline(y=0, color='red', linestyle='--', alpha=0.7)
        ax2.set_xlabel('Time (Days)')
        ax2.set_ylabel('Residuals (Observed - Predicted)')
        ax2.set_title('Petri Net Model Residual Analysis')
        ax2.grid(True, alpha=0.3)
        
        petri_mean_residual = np.mean(petri_residuals)
        petri_std_residual = np.std(petri_residuals)
        ax2.text(0.05, 0.85, f'Mean: {petri_mean_residual:.3f}\nStd: {petri_std_residual:.3f}', 
                transform=ax2.transAxes, verticalalignment='top')
        
        ax3 = axes[1, 0]
        ax3.text(0.05, 0.95, 'c', transform=ax3.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        ode_prediction_interval = 1.96 * std_residual
        petri_prediction_interval = 1.96 * petri_std_residual
        
        ax3.plot(petri_time, observed_cases, 'ko-', label='Observed Values', markersize=4, alpha=0.7)
        ax3.plot(petri_time, ode_daily, 'b-', label='ODE Prediction', linewidth=2)
        ax3.fill_between(petri_time, 
                       ode_daily - ode_prediction_interval, 
                       ode_daily + ode_prediction_interval, 
                       alpha=0.2, color='blue', label='ODE 95% Prediction Interval')
        
        ax3.plot(petri_time, petri_daily, 'g--', label='Petri Net Prediction', linewidth=2)
        ax3.fill_between(petri_time, 
                       petri_daily - petri_prediction_interval, 
                       petri_daily + petri_prediction_interval, 
                       alpha=0.2, color='green', label='Petri Net 95% Prediction Interval')
        
        ax3.set_xlabel('Time (Days)')
        ax3.set_ylabel('Daily Incidence')
        ax3.set_title('Prediction Interval Analysis')
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)
        
        ax4 = axes[1, 1]
        ax4.text(0.05, 0.95, 'd', transform=ax4.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        r_comparisons = self.comparison_data.get('r_value_comparisons', {})
        
        if r_comparisons:
            phases = ['Initial Phase', 'Intervention 1', 'Intervention 2', 'Intervention 3']
            phase_keys = ['phase_0_initial', 'phase_1_intervention', 'phase_2_intervention', 'phase_3_intervention']
            
            ngm_ode_r = []
            ngm_petri_r = []
            ngm_ode_ci_width = []
            ngm_petri_ci_width = []
            
            for phase_key in phase_keys:
                if phase_key in r_comparisons:
                    phase_data = r_comparisons[phase_key]
                    
                    ngm_ode_r.append(phase_data['ngm_method']['ode_r'])
                    ngm_petri_r.append(phase_data['ngm_method']['petri_r'])
                    ngm_ode_ci_width.append(phase_data['ngm_method']['ode_ci_width'])
                    ngm_petri_ci_width.append(phase_data['ngm_method']['petri_ci_width'])
            
            x = np.arange(len(phases))
            width = 0.35
            
            bars1 = ax4.bar(x - width/2, ngm_ode_r, width, label='ODE Model',
                           yerr=[ci/2 for ci in ngm_ode_ci_width], capsize=5, alpha=0.8)
            bars2 = ax4.bar(x + width/2, ngm_petri_r, width, label='Petri Net Model',
                           yerr=[ci/2 for ci in ngm_petri_ci_width], capsize=5, alpha=0.8)
            
            ax4.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='R=1 Threshold')
            ax4.set_xlabel('Intervention Phase')
            ax4.set_ylabel('Reproduction Number R$_0$')
            ax4.set_title('Phase-wise R$_0$ Value Comparison')
            ax4.set_xticks(x)
            ax4.set_xticklabels(phases, rotation=45)
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        plt.savefig(f'{self.output_dir}/model_validation_diagnostics.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{self.output_dir}/model_validation_diagnostics.pdf', bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    plotter = ChikungunyaPlotter()
    plotter.generate_all_plots()