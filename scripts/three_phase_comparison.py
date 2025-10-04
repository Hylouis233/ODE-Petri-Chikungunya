import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import sys
import os
from typing import Dict, List, Tuple, Any, Optional
from petri import CHIKVPetriNetModel
from petri_r_value_calculator import PetriNetRCalculator
from petri_three_phase_intervention_fitting import fit_three_phase_petri_model, load_foshan_data_with_phases
from three_phase_intervention_fitting import fit_three_phase_model
from model_parameters import CHIKVModelParameters
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class ThreePhaseComparator:
    def __init__(self):
        self.params = CHIKVModelParameters()
        self.foshan_data = None
        self.phase_days = None
        self.comparison_results = {}
        
    def load_data(self):
        self.days_observed, self.cases_observed, self.phase_days = load_foshan_data_with_phases()
        return self.days_observed, self.cases_observed, self.phase_days
    
    def fit_both_models(self):
        ode_params, _, _, _, ode_r_results = fit_three_phase_model()
        ode_fitted = True
        petri_params, _, _, _, petri_r_results = fit_three_phase_petri_model()
        petri_fitted = True
        fitting_results = {
            'ode_fitted': ode_fitted,
            'petri_fitted': petri_fitted,
            'ode_parameters': ode_params,
            'petri_parameters': petri_params,
            'ode_r_results': ode_r_results,
            'petri_r_results': petri_r_results
        }
        return fitting_results
    
    def compare_parameter_fitting(self, fitting_results: Dict[str, Any]) -> Dict[str, Any]:
        if not fitting_results['ode_fitted'] or not fitting_results['petri_fitted']:
            return {'comparison_skipped': True}
        ode_params = fitting_results['ode_parameters']
        petri_params = fitting_results['petri_parameters']
        ode_betas = list(ode_params)
        petri_betas = list(petri_params)
        param_differences = {
            'beta_0': {
                'ode': ode_betas[0],
                'petri': petri_betas[0],
                'absolute_diff': abs(ode_betas[0] - petri_betas[0]),
                'relative_diff_percent': abs(ode_betas[0] - petri_betas[0]) / ode_betas[0] * 100
            },
            'beta_1': {
                'ode': ode_betas[1],
                'petri': petri_betas[1],
                'absolute_diff': abs(ode_betas[1] - petri_betas[1]),
                'relative_diff_percent': abs(ode_betas[1] - petri_betas[1]) / ode_betas[1] * 100
            },
            'beta_2': {
                'ode': ode_betas[2],
                'petri': petri_betas[2],
                'absolute_diff': abs(ode_betas[2] - petri_betas[2]),
                'relative_diff_percent': abs(ode_betas[2] - petri_betas[2]) / ode_betas[2] * 100
            },
            'beta_3': {
                'ode': ode_betas[3],
                'petri': petri_betas[3],
                'absolute_diff': abs(ode_betas[3] - petri_betas[3]),
                'relative_diff_percent': abs(ode_betas[3] - petri_betas[3]) / ode_betas[3] * 100
            }
        }
        avg_relative_diff = np.mean([values['relative_diff_percent'] for values in param_differences.values()])
        param_comparison = {
            'parameter_differences': param_differences,
            'average_relative_difference_percent': avg_relative_diff
        }
        return param_comparison
    
    def compare_r_value_calculations(self, fitting_results: Dict[str, Any]) -> Dict[str, Any]:
        if not fitting_results['ode_fitted'] or not fitting_results['petri_fitted']:
            return {'comparison_skipped': True}
        ode_r_results = fitting_results['ode_r_results']
        petri_r_results = fitting_results['petri_r_results']
        phases = ['phase_0_initial', 'phase_1_intervention', 'phase_2_intervention', 'phase_3_intervention']
        r_value_comparisons = {}
        for phase in phases:
            ode_phase = ode_r_results['three_phase_results_with_ci'][phase]
            ode_point = ode_phase['point_estimates']
            if 'r_definition' in ode_point:
                ode_def = ode_point['r_definition']
                ode_ngm = ode_point['r_next_generation_matrix']
            else:
                ode_def = ode_point['results']['definition_method']
                ode_ngm = ode_point['results']['next_generation_matrix']
            ode_beta = ode_point.get('beta', np.nan)
            petri_phase = petri_r_results['petri_net_results_with_ci'][phase]
            petri_point = petri_phase['point_estimates']
            petri_ngm = petri_point['r_next_generation_matrix']
            petri_beta = petri_point.get('beta', np.nan)
            ngm_diff = abs(ode_ngm - petri_ngm) / ode_ngm * 100
            beta_diff = abs(ode_beta - petri_beta) / ode_beta * 100
            ode_def_ci = ode_phase['confidence_intervals']['definition_method']
            ode_ngm_ci = ode_phase['confidence_intervals']['next_generation_matrix']
            petri_ngm_ci = petri_phase['confidence_intervals']['next_generation_matrix']
            ode_def_width = ode_def_ci['ci_upper'] - ode_def_ci['ci_lower']
            ode_ngm_width = ode_ngm_ci['ci_upper'] - ode_ngm_ci['ci_lower']
            petri_ngm_width = petri_ngm_ci['ci_upper'] - petri_ngm_ci['ci_lower']
            r_value_comparisons[phase] = {
                'beta_comparison': {
                    'ode': ode_beta,
                    'petri': petri_beta,
                    'difference_percent': beta_diff
                },
                'definition_method': {
                    'ode_r': ode_def,
                    'ode_ci_width': ode_def_width
                },
                'ngm_method': {
                    'ode_r': ode_ngm,
                    'petri_r': petri_ngm,
                    'difference_percent': ngm_diff,
                    'ode_ci_width': ode_ngm_width,
                    'petri_ci_width': petri_ngm_width
                }
            }
        return {'r_value_comparisons': r_value_comparisons}
    
    def analyze_intervention_effectiveness(self, fitting_results: Dict[str, Any]) -> Dict[str, Any]:
        if not fitting_results['ode_fitted'] or not fitting_results['petri_fitted']:
            return {'analysis_skipped': True}
        ode_params = fitting_results['ode_parameters']
        petri_params = fitting_results['petri_parameters']
        def calculate_reduction_effects(params):
            beta_0, beta_1, beta_2, beta_3 = params
            return {
                'phase_1_reduction': (1 - beta_1/beta_0) * 100,
                'phase_2_reduction': (1 - beta_2/beta_0) * 100,
                'phase_3_reduction': (1 - beta_3/beta_0) * 100,
                'cumulative_reduction': (1 - beta_3/beta_0) * 100
            }
        ode_reductions = calculate_reduction_effects(ode_params)
        petri_reductions = calculate_reduction_effects(petri_params)
        model_params = CHIKVModelParameters()
        p = model_params.p
        gamma = model_params.gamma
        epsilon = model_params.epsilon
        effective_infectious_period = p/gamma + (1-p)*epsilon/gamma
        def calculate_r_progression(params):
            beta_0, beta_1, beta_2, beta_3 = params
            return [
                beta_0 * effective_infectious_period,
                beta_1 * effective_infectious_period,
                beta_2 * effective_infectious_period,
                beta_3 * effective_infectious_period
            ]
        ode_r_progression = calculate_r_progression(ode_params)
        petri_r_progression = calculate_r_progression(petri_params)
        intervention_analysis = {
            'ode_reduction_effects': ode_reductions,
            'petri_reduction_effects': petri_reductions,
            'ode_r_progression': ode_r_progression,
            'petri_r_progression': petri_r_progression,
            'control_effectiveness_comparison': {
                'phase_1_diff': abs(ode_reductions['phase_1_reduction'] - petri_reductions['phase_1_reduction']),
                'phase_2_diff': abs(ode_reductions['phase_2_reduction'] - petri_reductions['phase_2_reduction']),
                'phase_3_diff': abs(ode_reductions['phase_3_reduction'] - petri_reductions['phase_3_reduction'])
            }
        }
        return intervention_analysis
    
    def save_three_phase_comparison_plotting_data(self, fitting_results: Dict[str, Any], 
                                                r_comparison: Dict[str, Any], 
                                                intervention_analysis: Dict[str, Any]):
        plotting_data = {
            'model_type': 'three_phase_comparison',
            'saved_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
            'fitting_results': {
                'ode_fitted': fitting_results.get('ode_fitted', False),
                'petri_fitted': fitting_results.get('petri_fitted', False),
                'ode_parameters': fitting_results.get('ode_parameters', []),
                'petri_parameters': fitting_results.get('petri_parameters', [])
            },
            'r_value_comparisons': r_comparison.get('r_value_comparisons', {}),
            'intervention_analysis': {
                'ode_reduction_effects': intervention_analysis.get('ode_reduction_effects', {}),
                'petri_reduction_effects': intervention_analysis.get('petri_reduction_effects', {}),
                'ode_r_progression': intervention_analysis.get('ode_r_progression', []),
                'petri_r_progression': intervention_analysis.get('petri_r_progression', [])
            }
        }
        plotting_data_path = os.path.join(os.path.dirname(current_dir), 'data', 'three_phase_comparison_plotting_data.json')
        with open(plotting_data_path, 'w', encoding='utf-8') as f:
            json.dump(plotting_data, f, ensure_ascii=False, indent=2, default=str)
    
    def run_complete_three_phase_analysis(self) -> Dict[str, Any]:
        self.load_data()
        fitting_results = self.fit_both_models()
        if not fitting_results['ode_fitted'] or not fitting_results['petri_fitted']:
            return {'analysis_failed': True}
        param_comparison = self.compare_parameter_fitting(fitting_results)
        r_comparison = self.compare_r_value_calculations(fitting_results)
        intervention_analysis = self.analyze_intervention_effectiveness(fitting_results)
        complete_results = {
            'analysis_type': 'three_phase_intervention_comparison',
            'data_source': 'foshan_chikv_outbreak_2025',
            'fitting_results': fitting_results,
            'parameter_comparison': param_comparison,
            'r_value_comparison': r_comparison,
            'intervention_analysis': intervention_analysis
        }
        return complete_results
    
    def save_results(self, results: Dict[str, Any], output_path: str):
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2, default=str)

if __name__ == "__main__":
    comparator = ThreePhaseComparator()
    results = comparator.run_complete_three_phase_analysis()
    if 'analysis_failed' not in results:
        results_path = os.path.join(os.path.dirname(current_dir), 'data', 'three_phase_comparison_results.json')
        comparator.save_results(results, results_path)