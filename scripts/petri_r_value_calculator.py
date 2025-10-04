import numpy as np
from typing import Dict, Any
from model_parameters import CHIKVModelParameters
from petri import CHIKVPetriNetModel


class PetriNetRCalculator:

    def __init__(self, model: CHIKVPetriNetModel):
        self.model = model
        self.params = model.p

    def calculate_r_next_generation_matrix(self, beta: float) -> float:
        p = self.params['p']
        sigma = self.params['sigma']
        gamma = self.params['gamma']
        mu = self.params['mu']
        epsilon = self.params['epsilon']
        q = self.params['q']
        omega = self.params['omega']

        F = np.array([
            [0, beta, beta * epsilon, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ])

        V = np.array([
            [sigma + mu, 0, 0, 0],
            [-p * sigma, gamma + mu, 0, 0],
            [-(1 - p) * sigma, 0, gamma + mu, 0],
            [0, -q * gamma, 0, omega + mu]
        ])

        V_inv = np.linalg.inv(V)
        K = np.dot(F, V_inv)
        eigenvalues = np.linalg.eigvals(K)
        R0 = float(np.max(np.real(eigenvalues)))
        return R0

    def calculate_phase_r_values(self, beta_0: float, beta1: float, beta2: float,
                                 beta3: float, confidence_level: float = 0.95,
                                 n_samples: int = 1000) -> Dict[str, Any]:
        betas = [beta_0, beta1, beta2, beta3]
        phase_names = ['phase_0_initial', 'phase_1_intervention',
                       'phase_2_intervention', 'phase_3_intervention']

        results: Dict[str, Any] = {}
        for beta, phase_name in zip(betas, phase_names):
            r_ngm = self.calculate_r_next_generation_matrix(beta)

            ngm_samples = []
            for _ in range(n_samples):
                beta_sample = np.random.normal(beta, beta * 0.05)
                beta_sample = max(0.001, beta_sample)
                ngm_samples.append(self.calculate_r_next_generation_matrix(beta_sample))

            ngm_samples = np.array(ngm_samples)

            alpha = 1 - confidence_level
            ngm_ci_lower = np.percentile(ngm_samples, 100 * alpha / 2)
            ngm_ci_upper = np.percentile(ngm_samples, 100 * (1 - alpha / 2))

            results[phase_name] = {
                'point_estimates': {
                    'beta': beta,
                    'r_next_generation_matrix': r_ngm
                },
                'confidence_intervals': {
                    'next_generation_matrix': {
                        'estimate': r_ngm,
                        'ci_lower': ngm_ci_lower,
                        'ci_upper': ngm_ci_upper,
                        'ci_range': f"[{ngm_ci_lower:.3f}, {ngm_ci_upper:.3f}]"
                    }
                }
            }

        return {'petri_net_results_with_ci': results}