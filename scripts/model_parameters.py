import json
import numpy as np

class CHIKVModelParameters:
    
    def __init__(self):
        self.N = 9698900
        self.Lambda = 229
        self.mu = 1.381e-5
        self.sigma = 0.0714
        self.p = 0.75  
        self.gamma = 0.143
        self.q = 0.43
        self.omega = 0.0026
        self.beta = 0.502
        self.epsilon = 0.375
        self.effective_contact_rate = self._estimate_contact_rate()
    
    def _estimate_contact_rate(self):
        attack_rate = 4806 / self.N
        duration = 42
        estimated_R0 = 2.87
        return estimated_R0 * self.gamma
    
    def get_parameter_dict(self):
        return {
            'N': self.N,
            'Lambda': self.Lambda,
            'mu': self.mu,
            'sigma': self.sigma,
            'p': self.p,
            'gamma': self.gamma,
            'q': self.q,
            'omega': self.omega,
            'beta': self.beta,
            'epsilon': self.epsilon
        }
    
    def calculate_basic_reproduction_number(self):
        effective_infectious_period = (
            self.p * (1/self.gamma) + 
            (1-self.p) * self.epsilon * (1/self.gamma)
        )
        R0 = self.beta * effective_infectious_period
        return R0


if __name__ == "__main__":
    params = CHIKVModelParameters()
    param_dict = params.get_parameter_dict()
    with open('../data/simplified_model_parameters.json', 'w', encoding='utf-8') as f:
        json.dump(param_dict, f, indent=4, ensure_ascii=False)