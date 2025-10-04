import numpy as np
from scipy.linalg import eig

class RValueCalculator:
    
    def __init__(self, model_params):
        self.params = model_params
        
    def calculate_r_definition_method(self, beta):
        p = self.params.p
        gamma = self.params.gamma
        epsilon = self.params.epsilon
        
        effective_infectious_period = p/gamma + (1-p)*epsilon/gamma
        R0 = beta * effective_infectious_period
        
        return {
            'method': 'definition_method',
            'beta': beta,
            'R0': R0,
            'components': {
                'symptomatic_contribution': beta * (p/gamma),
                'asymptomatic_contribution': beta * ((1-p)*epsilon/gamma),
                'effective_infectious_period': effective_infectious_period
            },
            'parameters_used': {
                'p': p,
                'gamma': gamma,
                'epsilon': epsilon
            }
        }
    
    def calculate_r_next_generation_method(self, beta):
        sigma = self.params.sigma
        p = self.params.p
        gamma = self.params.gamma
        epsilon = self.params.epsilon
        q = self.params.q
        omega = self.params.omega
        mu = self.params.mu
        
        F = np.array([
            [0, beta, epsilon*beta, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ])
        
        V = np.array([
            [sigma + mu, 0, 0, 0],
            [-p*sigma, gamma + mu, 0, 0],
            [-(1-p)*sigma, 0, gamma + mu, 0],
            [0, -q*gamma, 0, omega + mu]
        ])
        
        V_inv = np.linalg.inv(V)
        K = F @ V_inv
        eigenvalues = eig(K)[0]
        R0_ngm = float(np.max(np.abs(eigenvalues)))
        
        return {
            'method': 'next_generation_matrix',
            'beta': beta,
            'R0_ngm': R0_ngm,
            'eigenvalues': eigenvalues.tolist()
        }
    
    def compare_methods(self, beta):
        def_result = self.calculate_r_definition_method(beta)
        ngm_result = self.calculate_r_next_generation_method(beta)
        
        R_def = def_result['R0']
        R_ngm = ngm_result['R0_ngm']
        diff_percent = abs(R_def - R_ngm) / R_ngm * 100
        
        return {
            'beta': beta,
            'results': {
                'definition_method': R_def,
                'next_generation_matrix': R_ngm
            },
            'differences': {
                'definition_vs_ngm_percent': diff_percent
            },
            'detailed_results': {
                'definition_method': def_result,
                'ngm_method': ngm_result
            }
        }
    
    def calculate_phase_r_values(self, beta_0, beta1, beta2, beta3, confidence_level=0.95, n_samples=1000):
        
        phases = {
            'phase_0_initial': beta_0,
            'phase_1_intervention': beta1, 
            'phase_2_intervention': beta2,
            'phase_3_intervention': beta3
        }
        
        results = {}
        
        for phase_name, beta in phases.items():
            results[phase_name] = self.calculate_r(
                beta, confidence_level, n_samples
            )
            
        return {
            'three_phase_results_with_ci': results,
            'intervention_effectiveness': {
                'phase_1_reduction': (1 - beta1/beta_0) * 100,
                'phase_2_reduction': (1 - beta2/beta_0) * 100,
                'phase_3_reduction': (1 - beta3/beta_0) * 100
            },
            'analysis_parameters': {
                'confidence_level': confidence_level,
                'n_samples': n_samples
            }
        }
    
    def calculate_r(self, beta, confidence_level=0.95, n_samples=1000):
        
        point_estimates = self.compare_methods(beta)
        def_samples = []
        ngm_samples = []
        
        for _ in range(n_samples):
            beta_sample = np.random.normal(beta, beta * 0.05)
            beta_sample = max(0.001, beta_sample)
            
            def_result = self.calculate_r_definition_method(beta_sample)
            def_samples.append(def_result['R0'])
            
            ngm_result = self.calculate_r_next_generation_method(beta_sample)
            ngm_samples.append(ngm_result['R0_ngm'])
        
        def_samples = np.array(def_samples)
        ngm_samples = np.array(ngm_samples)
        alpha = 1 - confidence_level
        def_ci_lower = np.percentile(def_samples, 100 * alpha / 2)
        def_ci_upper = np.percentile(def_samples, 100 * (1 - alpha / 2))
        ngm_ci_lower = np.percentile(ngm_samples, 100 * alpha / 2)
        ngm_ci_upper = np.percentile(ngm_samples, 100 * (1 - alpha / 2))
        
        return {
            'point_estimates': point_estimates,
            'confidence_intervals': {
                'confidence_level': confidence_level,
                'n_samples': len(def_samples),
                'definition_method': {
                    'estimate': point_estimates['results']['definition_method'],
                    'ci_lower': def_ci_lower,
                    'ci_upper': def_ci_upper,
                    'ci_range': f"[{def_ci_lower:.3f}, {def_ci_upper:.3f}]"
                },
                'next_generation_matrix': {
                    'estimate': point_estimates['results']['next_generation_matrix'],
                    'ci_lower': ngm_ci_lower,
                    'ci_upper': ngm_ci_upper,
                    'ci_range': f"[{ngm_ci_lower:.3f}, {ngm_ci_upper:.3f}]"
                }
            }
        }