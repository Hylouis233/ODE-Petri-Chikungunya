import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from scipy.stats import pearsonr
import sys
import os
import json
from typing import Dict, List, Tuple, Any
from scripts.petri import CHIKVPetriNetModel
from scripts.petri_r_value_calculator import PetriNetRCalculator
from scripts.r_value_calculator import RValueCalculator
from scripts.model_parameters import CHIKVModelParameters

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['mathtext.fontset'] = 'dejavusans'
sns.set_style("whitegrid")

class SensitivityAnalyzer:
    
    def __init__(self):
        self.base_params = CHIKVModelParameters()
        
        self.petri_model = CHIKVPetriNetModel(self.base_params)
        self.petri_calculator = PetriNetRCalculator(self.petri_model)
        self.ode_calculator = RValueCalculator(self.base_params)
        
        self.load_experimental_r0_values()
        
        self.param_ranges = {
            'p': (0.3, 0.9, 'has_symptoms_probability'),
            'sigma': (0.1, 0.5, 'incubation_rate'),
            'gamma': (0.1, 0.5, 'recovery_rate'),
            'epsilon': (0.1, 0.8, 'asymptomatic_transmissibility'),
            'q': (0.05, 0.3, 'chronicity_probability'),
            'omega': (0.01, 0.1, 'chronic_recovery_rate'),
            'mu': (1e-5, 1e-3, 'death_rate')
        }
    
    def load_experimental_r0_values(self):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'petri_net_fitting_results.json'), 'r', encoding='utf-8') as f:
            petri_data = json.load(f)
        
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'three_phase_comparison_results.json'), 'r', encoding='utf-8') as f:
            ode_data = json.load(f)
        
        self.petri_r0_initial = petri_data['epidemiological_parameters']['R0_initial']
        self.petri_r0_ngm = petri_data['r_value_results']['petri_net_results_with_ci']['phase_0_initial']['point_estimates']['r_next_generation_matrix']
        
        self.ode_r0_initial = ode_data['fitting_results']['ode_r_results']['three_phase_results_with_ci']['phase_0_initial']['point_estimates']['results']['definition_method']
        self.ode_r0_ngm = ode_data['fitting_results']['ode_r_results']['three_phase_results_with_ci']['phase_0_initial']['point_estimates']['results']['next_generation_matrix']
        
        self.base_r0 = self.petri_r0_initial
    
    def single_parameter_sensitivity(self, param_name: str, n_points: int = 50) -> Dict[str, Any]:
        if param_name not in self.param_ranges:
            raise ValueError(f"Unknown parameter: {param_name}")
        
        min_val, max_val, description = self.param_ranges[param_name]
        base_val = getattr(self.base_params, param_name)
        
        param_values = np.linspace(min_val, max_val, n_points)
        
        results = {
            'parameter_name': param_name,
            'parameter_description': description,
            'base_value': base_val,
            'test_range': [min_val, max_val],
            'parameter_values': param_values.tolist(),
            'petri_net_results': {
                'ngm_method': []
            },
            'ode_results': {
                'definition_method': [],
                'ngm_method': []
            }
        }
        
        for param_val in param_values:
            test_params = CHIKVModelParameters()
            setattr(test_params, param_name, param_val)
            
            test_petri_model = CHIKVPetriNetModel(test_params)
            test_petri_calculator = PetriNetRCalculator(test_petri_model)
            
            test_ode_calculator = RValueCalculator(test_params)
            
            experimental_beta = 2.3559611446475555
            
            petri_ngm = test_petri_calculator.calculate_r_next_generation_matrix(experimental_beta)
            
            ode_def_result = test_ode_calculator.calculate_r_definition_method(experimental_beta)
            ode_ngm_result = test_ode_calculator.calculate_r_next_generation_method(experimental_beta)
            
            results['petri_net_results']['ngm_method'].append(petri_ngm)
            results['ode_results']['definition_method'].append(ode_def_result['R0'])
            results['ode_results']['ngm_method'].append(ode_ngm_result['R0_ngm'])
        
        sensitivity_metrics = self.calculate_sensitivity_metrics(results)
        results['sensitivity_metrics'] = sensitivity_metrics
        
        return results
    
    def calculate_sensitivity_metrics(self, results: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
        
        param_values = np.array(results['parameter_values'])
        base_value = results['base_value']
        
        metrics = {}
        
        for model_type in ['petri_net_results', 'ode_results']:
            for method_type in ['definition_method', 'ngm_method']:
                r_values = np.array(results[model_type][method_type])
                
                min_r = np.min(r_values)
                max_r = np.max(r_values)
                mean_r = np.mean(r_values)
                std_r = np.std(r_values)
                
                relative_change = (max_r - min_r) / mean_r * 100
                
                dr_dp = np.gradient(r_values, param_values)
                base_idx = np.argmin(np.abs(param_values - base_value))
                sensitivity_coeff = dr_dp[base_idx] * base_value / r_values[base_idx]
                
                correlation, p_value = pearsonr(param_values, r_values)
                
                method_key = f"{model_type}_{method_type}".replace('_results', '').replace('_method', '')
                
                metrics[method_key] = {
                    'min_r': float(min_r),
                    'max_r': float(max_r),
                    'mean_r': float(mean_r),
                    'std_r': float(std_r),
                    'relative_change_percent': float(relative_change),
                    'sensitivity_coefficient': float(sensitivity_coeff),
                    'correlation': float(correlation),
                    'correlation_p_value': float(p_value)
                }
        
        return metrics
    
    def multi_parameter_sensitivity(self, n_samples: int = 1000) -> Dict[str, Any]:
        param_names = list(self.param_ranges.keys())
        param_samples = {}
        
        for param in param_names:
            min_val, max_val, description = self.param_ranges[param]
            samples = np.random.uniform(min_val, max_val, n_samples)
            param_samples[param] = samples
        
        results = {
            'n_samples': n_samples,
            'parameter_ranges': self.param_ranges,
            'parameter_samples': param_samples,
            'r_value_results': {
                'petri_net_ngm': [],
                'ode_definition': [],
                'ode_ngm': []
            }
        }
        
        for i in range(n_samples):
            test_params = CHIKVModelParameters()
            
            for param in param_names:
                setattr(test_params, param, param_samples[param][i])
            
            test_petri_model = CHIKVPetriNetModel(test_params)
            test_petri_calculator = PetriNetRCalculator(test_petri_model)
            test_ode_calculator = RValueCalculator(test_params)
            
            experimental_beta = 2.3559611446475555
            
            petri_ngm = test_petri_calculator.calculate_r_next_generation_matrix(experimental_beta)
            
            ode_def_result = test_ode_calculator.calculate_r_definition_method(experimental_beta)
            ode_ngm_result = test_ode_calculator.calculate_r_next_generation_method(experimental_beta)
            
            results['r_value_results']['petri_net_ngm'].append(petri_ngm)
            results['r_value_results']['ode_definition'].append(ode_def_result['R0'])
            results['r_value_results']['ode_ngm'].append(ode_ngm_result['R0_ngm'])
        
        valid_indices = []
        for i in range(n_samples):
            if not any(np.isnan(results['r_value_results'][key][i]) for key in results['r_value_results']):
                valid_indices.append(i)
        
        for param in param_names:
            param_samples[param] = param_samples[param][valid_indices]
        
        for key in results['r_value_results']:
            results['r_value_results'][key] = [results['r_value_results'][key][i] for i in valid_indices]
        
        results['valid_samples'] = len(valid_indices)
        results['parameter_samples'] = {param: samples.tolist() for param, samples in param_samples.items()}
        
        global_sensitivity = self.calculate_global_sensitivity(param_samples, results['r_value_results'])
        results['global_sensitivity'] = global_sensitivity
        
        return results
    
    def calculate_global_sensitivity(self, param_samples: Dict[str, np.ndarray], 
                                   r_results: Dict[str, List[float]]) -> Dict[str, Dict[str, Dict[str, float]]]:
        
        sensitivity_results = {}
        
        for r_method, r_values in r_results.items():
            r_array = np.array(r_values)
            sensitivity_results[r_method] = {}
            
            total_variance = np.var(r_array)
            
            for param_name, param_values in param_samples.items():
                median_param = np.median(param_values)
                high_group = r_array[param_values > median_param]
                low_group = r_array[param_values <= median_param]
                
                group_means = [np.mean(low_group), np.mean(high_group)]
                between_group_variance = np.var(group_means) * len(group_means)
                
                sobol_index = between_group_variance / total_variance if total_variance > 0 else 0
                
                correlation, p_value = pearsonr(param_values, r_array)
                
                sensitivity_results[r_method][param_name] = {
                    'sobol_index': float(sobol_index),
                    'correlation': float(correlation),
                    'correlation_p_value': float(p_value)
                }
        
        return sensitivity_results
    
    def create_sensitivity_visualizations(self, single_param_results: Dict[str, Dict], 
                                    multi_param_results: Dict[str, Any]):
        fig_height = 10
        fig, axes = plt.subplots(2, 3, figsize=(18, fig_height))
        
        ax1 = axes[0, 0]
        
        main_params = ['gamma', 'p', 'epsilon']
        colors = plt.cm.tab10(np.linspace(0, 1, len(main_params)))
        
        color_map = {}
        for i, param_name in enumerate(main_params):
            if param_name in single_param_results:
                param_results = single_param_results[param_name]
            param_values = param_results['parameter_values']
            color_map[param_name] = colors[i]
            
            param_values_normalized = np.array(param_values)
            param_min, param_max = np.min(param_values_normalized), np.max(param_values_normalized)
            param_values_normalized = (param_values_normalized - param_min) / (param_max - param_min)
            
    
        ax1.set_xlabel('Normalized Parameter Values')
        ax1.set_ylabel('Basic Reproduction Number R$_0$')
        ax1.set_title('Single Parameter Sensitivity Analysis')
        ax1.set_xlim(-0.05, 1.15)
        ax1.set_xticks(np.arange(0, 1.2, 0.1))
        ax1.set_xticklabels([f'{x:.1f}' for x in np.arange(0, 1.2, 0.1)])
        ax1.legend(loc='lower right')
        ax1.grid(True, alpha=0.3)
        ax1.text(0.02, 0.98, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        ax2 = axes[0, 1]
        
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'petri_net_fitting_results.json'), 'r', encoding='utf-8') as f:
            petri_data = json.load(f)
        
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'three_phase_comparison_results.json'), 'r', encoding='utf-8') as f:
            ode_data = json.load(f)
        
        violin_data = []
        violin_labels = []
        violin_colors = []
        
        petri_ngm_r0 = petri_data['r_value_results']['petri_net_results_with_ci']['phase_0_initial']['point_estimates']['r_next_generation_matrix']
        petri_ngm_ci_data = petri_data['r_value_results']['petri_net_results_with_ci']['phase_0_initial']['confidence_intervals']['next_generation_matrix']
        petri_ngm_ci = [petri_ngm_ci_data['ci_lower'], petri_ngm_ci_data['ci_upper']]
        
        ode_def_r0 = ode_data['fitting_results']['ode_r_results']['three_phase_results_with_ci']['phase_0_initial']['point_estimates']['results']['definition_method']
        ode_def_ci_data = ode_data['fitting_results']['ode_r_results']['three_phase_results_with_ci']['phase_0_initial']['confidence_intervals']['definition_method']
        ode_def_ci = [ode_def_ci_data['ci_lower'], ode_def_ci_data['ci_upper']]
        
        ode_ngm_r0 = ode_data['fitting_results']['ode_r_results']['three_phase_results_with_ci']['phase_0_initial']['point_estimates']['results']['next_generation_matrix']
        ode_ngm_ci_data = ode_data['fitting_results']['ode_r_results']['three_phase_results_with_ci']['phase_0_initial']['confidence_intervals']['next_generation_matrix']
        ode_ngm_ci = [ode_ngm_ci_data['ci_lower'], ode_ngm_ci_data['ci_upper']]
        
        n_samples = 1000
        
        petri_ngm_mean = petri_ngm_r0
        petri_ngm_std = (petri_ngm_ci[1] - petri_ngm_ci[0]) / (2 * 1.96)
        petri_ngm_samples = np.random.normal(petri_ngm_mean, petri_ngm_std, n_samples)
        violin_data.append(petri_ngm_samples)
        violin_labels.append('Petri NGM')
        violin_colors.append('blue')
        
        ode_def_mean = ode_def_r0
        ode_def_std = (ode_def_ci[1] - ode_def_ci[0]) / (2 * 1.96)
        ode_def_samples = np.random.normal(ode_def_mean, ode_def_std, n_samples)
        violin_data.append(ode_def_samples)
        violin_labels.append('ODE Def')
        violin_colors.append('red')
        
        ode_ngm_mean = ode_ngm_r0
        ode_ngm_std = (ode_ngm_ci[1] - ode_ngm_ci[0]) / (2 * 1.96)
        ode_ngm_samples = np.random.normal(ode_ngm_mean, ode_ngm_std, n_samples)
        violin_data.append(ode_ngm_samples)
        violin_labels.append('ODE NGM')
        violin_colors.append('orange')
        
        parts = ax2.violinplot(violin_data, positions=range(len(violin_data)), 
                             showmeans=True, showmedians=True, showextrema=True)
        
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(violin_colors[i])
            pc.set_alpha(0.7)
            pc.set_edgecolor('black')
            pc.set_linewidth(1)
        
        parts['cmeans'].set_color('darkred')
        parts['cmeans'].set_linewidth(2)
        parts['cmedians'].set_color('darkblue')
        parts['cmedians'].set_linewidth(2)
        parts['cmaxes'].set_color('black')
        parts['cmins'].set_color('black')
        parts['cbars'].set_color('black')
        
        ax2.set_xticks(range(len(violin_labels)))
        ax2.set_xticklabels(violin_labels)
        
        ax2.set_ylabel('Basic Reproduction Number R$_0$')
        ax2.set_title('R$_0$ Distribution with Confidence Intervals')
        
        ax2.set_ylim(0, 25)
        
        ax2.grid(True, alpha=0.3, axis='y')
        
        ax2.text(0.02, 0.98, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        fig.delaxes(axes[0, 2])
        ax3 = fig.add_subplot(2, 3, 3, projection='3d')
        
        gamma_range = np.linspace(0.1, 0.5, 25)
        p_range = np.linspace(0.3, 0.9, 25)
        epsilon_range = np.linspace(0.2, 0.8, 25)
        
        experimental_beta = 2.3559611446475555
        
        def calculate_R0(gamma_val, p_val, epsilon_val):
            effective_infectious_period = p_val/gamma_val + (1-p_val)*epsilon_val/gamma_val
            return experimental_beta * effective_infectious_period
        
        Gamma, P = np.meshgrid(gamma_range, p_range)
        
        epsilon_layers = np.linspace(0.2, 0.8, 6)
        
        from matplotlib.colors import Normalize
        norm = Normalize(vmin=epsilon_range.min(), vmax=epsilon_range.max())
        colormap = plt.cm.plasma
        
        surfaces = []
        for i, epsilon_val in enumerate(epsilon_layers):
            R0_surface = np.zeros_like(Gamma)
            for j in range(Gamma.shape[0]):
                for k in range(Gamma.shape[1]):
                    gamma = Gamma[j, k]
                    p = P[j, k]
                    R0_surface[j, k] = calculate_R0(gamma, p, epsilon_val)
            
            color = colormap(norm(epsilon_val))
            
            alpha = 0.4 * (i / (len(epsilon_layers) - 1))
            
            surface = ax3.plot_surface(Gamma, P, R0_surface, 
                                     facecolors=np.full(Gamma.shape + (4,), color), 
                                     alpha=alpha,
                                     linewidth=0.1, 
                                     shade=True, rstride=2, cstride=2)
            surfaces.append(surface)
        
        epsilon_ref = 0.5
        R0_ref_surface = np.zeros_like(Gamma)
        
        Epsilon_gradient = np.zeros_like(Gamma)
        for i in range(Gamma.shape[0]):
            for j in range(Gamma.shape[1]):
                gamma = Gamma[i, j]
                p = P[i, j]
                epsilon_pos = 0.6 * ((gamma - 0.1) / (0.5 - 0.1))
                Epsilon_gradient[i, j] = epsilon_pos
                R0_ref_surface[i, j] = calculate_R0(gamma, p, epsilon_pos)
        
        ref_surface = ax3.plot_surface(Gamma, P, R0_ref_surface, 
                                     facecolors=colormap(norm(Epsilon_gradient)),
                                     alpha=0.9, linewidth=0.3, 
                                     shade=False,
                                     rstride=1, cstride=1)
        
        ax3.set_xlabel('γ', fontsize=10, labelpad=8)
        ax3.set_ylabel('p', fontsize=10, labelpad=8)
        ax3.set_zlabel('R$_0$ value', fontsize=10, labelpad=8)
        ax3.set_title('Sensitivity Analysis of R$_0$', fontsize=11, pad=15)
        
        ax3.set_xlim(0.1, 0.5)
        ax3.set_ylim(0.3, 0.9)
        
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=colormap)
        mappable.set_array(Epsilon_gradient)
        cbar = plt.colorbar(mappable, ax=ax3, shrink=0.6, aspect=25, pad=0.05)
        cbar.set_label('ε', rotation=270, labelpad=18, fontsize=10)
        cbar.ax.tick_params(labelsize=8)
        
        ax3.view_init(elev=25, azim=45)
        ax3.grid(True, alpha=0.3)
        
        ax3.set_xticks(np.linspace(0.1, 0.5, 5))
        ax3.set_yticks(np.linspace(0.3, 0.9, 4))
        ax3.tick_params(axis='both', which='major', labelsize=8)
        
        ax3.tick_params(axis='x', labelsize=8)
        ax3.tick_params(axis='y', labelsize=8)
        ax3.tick_params(axis='z', labelsize=8)
        
        ax3.grid(True, alpha=0.3)
        
        ax3.text2D(0.02, 0.98, 'c', transform=ax3.transAxes, fontsize=14, 
                  fontweight='bold', va='top', ha='left')
        
        if 'gamma' in single_param_results:
            ax4 = axes[1, 0]
            param_results = single_param_results['gamma']
            param_values = param_results['parameter_values']
            param_values_array = np.array(param_values)
            
            petri_ngm = param_results['petri_net_results']['ngm_method']
            param_color = color_map['gamma']
            
            ax4.plot(param_values_array, petri_ngm, color=param_color, linewidth=2, label='Petri Net NGM')
            
            if 'r_value_results' in multi_param_results and 'parameter_samples' in multi_param_results:
                param_samples = multi_param_results['parameter_samples']
                r_values = multi_param_results['r_value_results']['petri_net_ngm']
                gamma_samples = np.array(param_samples['gamma'])
                ax4.scatter(gamma_samples, r_values, alpha=0.3, s=5, color=param_color)
            
            ax4.axvline(x=self.base_params.gamma, color='black', linestyle='--', linewidth=2, alpha=0.8, label='True γ')
            
            ax4.set_xlabel('gamma')
            ax4.set_ylabel('R$_0$')
            ax4.grid(True, alpha=0.3)
            ax4.text(0.02, 0.98, 'd', transform=ax4.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        if 'epsilon' in single_param_results:
            ax5 = axes[1, 1]
            param_results = single_param_results['epsilon']
            param_values = param_results['parameter_values']
            param_values_array = np.array(param_values)
            
            petri_ngm = param_results['petri_net_results']['ngm_method']
            param_color = color_map['epsilon']
            
            ax5.plot(param_values_array, petri_ngm, color=param_color, linewidth=2, label='Petri Net NGM')
            
            if 'r_value_results' in multi_param_results and 'parameter_samples' in multi_param_results:
                param_samples = multi_param_results['parameter_samples']
                r_values = multi_param_results['r_value_results']['petri_net_ngm']
                epsilon_samples = np.array(param_samples['epsilon'])
                ax5.scatter(epsilon_samples, r_values, alpha=0.3, s=5, color=param_color)
            
            ax5.axvline(x=self.base_params.epsilon, color='black', linestyle='--', linewidth=2, alpha=0.8, label='True ε')
            
            ax5.set_xlabel('epsilon')
            ax5.set_ylabel('R$_0$')
            ax5.grid(True, alpha=0.3)
            ax5.text(0.02, 0.98, 'e', transform=ax5.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        if 'p' in single_param_results:
            ax6 = axes[1, 2]
            param_results = single_param_results['p']
            param_values = param_results['parameter_values']
            param_values_array = np.array(param_values)
            
            petri_ngm = param_results['petri_net_results']['ngm_method']
            param_color = color_map['p']
            
            ax6.plot(param_values_array, petri_ngm, color=param_color, linewidth=2, label='Petri Net NGM')
            
            if 'r_value_results' in multi_param_results and 'parameter_samples' in multi_param_results:
                param_samples = multi_param_results['parameter_samples']
                r_values = multi_param_results["r_value_results"]['petri_net_ngm']
                p_samples = np.array(param_samples['p'])
                ax6.scatter(p_samples, r_values, alpha=0.3, s=5, color=param_color)
            
            ax6.axvline(x=self.base_params.p, color='black', linestyle='--', linewidth=2, alpha=0.8, label='True p')
            
            ax6.set_xlabel('p')
            ax6.set_ylabel('R$_0$')
            ax6.grid(True, alpha=0.3)
            ax6.text(0.02, 0.98, 'f', transform=ax6.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
        
        plt.tight_layout()
        
        output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'sensitivity_analysis.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        pdf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'sensitivity_analysis.pdf')
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        
        self.save_sensitivity_plotting_data(single_param_results, multi_param_results)
        
        plt.show()
    
    def save_sensitivity_plotting_data(self, single_param_results: Dict[str, Dict], 
                                     multi_param_results: Dict[str, Any]):
        
        plotting_data = {
            'model_type': 'sensitivity_analysis',
            'saved_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
            'description': 'Sensitivity analysis plotting data: Petri Net vs ODE model R-value calculation',
            
            'single_parameter': {},
            
            'multi_parameter': {
                'r_values': multi_param_results['r_value_results'],
                'global_metrics': multi_param_results['global_sensitivity'],
                'param_samples': multi_param_results['parameter_samples']
            },
            'statistics': {
                'n_single_params': len(single_param_results),
                'n_multi_samples': len(multi_param_results['r_value_results']['petri_net_definition']),
                'valid_sample_ratio': len(multi_param_results['r_value_results']['petri_net_definition']) / 1000
            }
        }
        
        for param_name, param_result in single_param_results.items():
            plotting_data['single_parameter'][param_name] = {
                'parameter_values': param_result['parameter_values'],
                'petri_net_results': {
                    'ngm_method': param_result['petri_net_results']['ngm_method']
                },
                'ode_results': {
                    'definition_method': param_result['ode_results']['definition_method'],
                    'ngm_method': param_result['ode_results']['ngm_method']
                },
                'sensitivity_metrics': param_result['sensitivity_metrics']
            }
        
        plotting_data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'sensitivity_analysis_plotting_data.json')
        
        with open(plotting_data_path, 'w', encoding='utf-8') as f:
            json.dump(plotting_data, f, ensure_ascii=False, indent=2, default=str)
    
    def run_complete_analysis(self):
        
        all_results = {
            'analysis_type': 'sensitivity_analysis',
            'experimental_beta': 2.3559611446475555,
            'experimental_r0_values': {
                'petri_net_ngm': self.petri_r0_ngm,
                'ode_definition': self.ode_r0_initial,
                'ode_ngm': self.ode_r0_ngm,
                'base_r0': self.base_r0
            },
            'parameter_ranges': self.param_ranges
        }
        
        single_param_results = {}
        
        for param_name in self.param_ranges.keys():
            param_result = self.single_parameter_sensitivity(param_name, n_points=50)
            single_param_results[param_name] = param_result
        
        all_results['single_parameter_results'] = single_param_results
        
        multi_param_results = self.multi_parameter_sensitivity(n_samples=1000)
        all_results['multi_parameter_results'] = multi_param_results
        
        self.create_sensitivity_visualizations(single_param_results, multi_param_results)
        
        sensitivity_summary = self.generate_sensitivity_summary(all_results)
        all_results['sensitivity_summary'] = sensitivity_summary
        
        return all_results
    
    def generate_sensitivity_summary(self, results: Dict[str, Any]) -> Dict[str, Any]:
        
        summary = {
            'most_sensitive_parameters': {},
            'method_consistency': {},
            'key_findings': {}
        }
        
        if 'multi_parameter_results' in results and 'global_sensitivity' in results['multi_parameter_results']:
            global_sens = results['multi_parameter_results']['global_sensitivity']
            
            for method, sens_data in global_sens.items():
                sorted_params = sorted(sens_data.items(), key=lambda x: abs(x[1]['sobol_index']), reverse=True)
                summary['most_sensitive_parameters'][method] = {
                    'rank_1': sorted_params[0][0] if sorted_params else None,
                    'rank_2': sorted_params[1][0] if len(sorted_params) > 1 else None,
                    'rank_3': sorted_params[2][0] if len(sorted_params) > 2 else None
                }
        
        if 'single_parameter_results' in results:
            single_results = results['single_parameter_results']
            consistency_scores = []
            
            for param_name, param_result in single_results.items():
                metrics = param_result['sensitivity_metrics']
                ode_coeff = metrics['ode_definition']['sensitivity_coefficient']
                
                pass
            
            summary['method_consistency'] = {
                'average_consistency': np.mean(consistency_scores) if consistency_scores else 0,
                'consistency_scores': consistency_scores
            }
        
        summary['key_findings'] = {
            'highest_impact_parameter': summary['most_sensitive_parameters'].get('petri_net_ngm', {}).get('rank_1'),
            'method_agreement': summary['method_consistency'].get('average_consistency', 0) > 0.9,
            'analysis_reliability': 'high' if len(results.get('multi_parameter_results', {}).get('r_value_results', {}).get('petri_net_ngm', [])) > 800 else 'medium'
        }
        
        return summary
    
    def save_results(self, results: Dict[str, Any], output_path: str):
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2, default=str)


if __name__ == "__main__":
    analyzer = SensitivityAnalyzer()
    
    results = analyzer.run_complete_analysis()
    
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'sensitivity_analysis_results.json')
    analyzer.save_results(results, output_path)