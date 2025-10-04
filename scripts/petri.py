import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import sys
import os
from typing import Dict, Tuple, List, Any, Optional

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from model_parameters import CHIKVModelParameters

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class PetriPlace:
    
    def __init__(self, name: str, initial_tokens: float = 0.0):
        self.name = name
        self.tokens = initial_tokens
        self.initial_tokens = initial_tokens
    
    def add_tokens(self, amount: float):
        self.tokens = max(0, self.tokens + amount)
    
    def remove_tokens(self, amount: float):
        self.tokens = max(0, self.tokens - amount)
    
    def get_tokens(self) -> float:
        return self.tokens
    
    def reset(self):
        self.tokens = self.initial_tokens

class PetriTransition:
    
    def __init__(self, name: str, rate_function, input_places: List[PetriPlace], 
                 output_places: List[PetriPlace], input_weights: List[float], 
                 output_weights: List[float]):
        self.name = name
        self.rate_function = rate_function
        self.input_places = input_places
        self.output_places = output_places
        self.input_weights = input_weights
        self.output_weights = output_weights
    
    def is_enabled(self) -> bool:
        for place, weight in zip(self.input_places, self.input_weights):
            if place.get_tokens() < weight:
                return False
        return True
    
    def fire(self, dt: float, rate_params: Dict[str, float]):
        if not self.is_enabled():
            return 0.0
        
        rate = self.rate_function(self.input_places, self.output_places, rate_params)
        firing_amount = rate * dt
        
        max_firings = float('inf')
        for place, weight in zip(self.input_places, self.input_weights):
            if weight > 0:
                max_firings = min(max_firings, place.get_tokens() / weight)
        
        actual_firing = min(firing_amount, max_firings)
        
        for place, weight in zip(self.input_places, self.input_weights):
            place.remove_tokens(actual_firing * weight)
        
        for place, weight in zip(self.output_places, self.output_weights):
            place.add_tokens(actual_firing * weight)
        
        return actual_firing

class CHIKVPetriNetModel:
    
    def __init__(self, parameters=None):
        if parameters is None:
            self.params = CHIKVModelParameters()
        else:
            self.params = parameters
        
        self.p = self.params.get_parameter_dict()
        
        self.places = {}
        self.transitions = {}
        
        self._initialize_places()
        self._initialize_transitions()
    
    def _initialize_places(self):
        self.places['S'] = PetriPlace('Susceptible', 0.0)
        self.places['E'] = PetriPlace('Exposed', 0.0)
        self.places['I_s'] = PetriPlace('Infected_Symptomatic', 0.0)
        self.places['I_a'] = PetriPlace('Infected_Asymptomatic', 0.0)
        self.places['C'] = PetriPlace('Chronic', 0.0)
        self.places['R'] = PetriPlace('Recovered', 0.0)
    
    def _initialize_transitions(self):
        
        def infection_rate(input_places, output_places, params):
            S = input_places[0].get_tokens()
            I_s = self.places['I_s'].get_tokens()
            I_a = self.places['I_a'].get_tokens()
            N = params.get('N', self.p['N'])
            beta = params.get('beta', 3.0)
            epsilon = params.get('epsilon', self.p['epsilon'])
            
            if N == 0:
                return 0.0
            
            foi = beta * S * (I_s + epsilon * I_a) / N
            return foi
        
        self.transitions['infection'] = PetriTransition(
            name='Infection',
            rate_function=infection_rate,
            input_places=[self.places['S']],
            output_places=[self.places['E']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def birth_rate(input_places, output_places, params):
            return params.get('Lambda', self.p['Lambda'])
        
        self.transitions['birth'] = PetriTransition(
            name='Birth',
            rate_function=birth_rate,
            input_places=[],
            output_places=[self.places['S']],
            input_weights=[],
            output_weights=[1.0]
        )
        
        def s_death_rate(input_places, output_places, params):
            S = input_places[0].get_tokens()
            mu = params.get('mu', self.p['mu'])
            return mu * S
        
        self.transitions['s_death'] = PetriTransition(
            name='S_Death',
            rate_function=s_death_rate,
            input_places=[self.places['S']],
            output_places=[],
            input_weights=[1.0],
            output_weights=[]
        )
        
        def e_to_is_rate(input_places, output_places, params):
            E = input_places[0].get_tokens()
            sigma = params.get('sigma', self.p['sigma'])
            p_val = params.get('p', self.p['p'])
            return sigma * p_val * E
        
        self.transitions['e_to_is'] = PetriTransition(
            name='E_to_Is',
            rate_function=e_to_is_rate,
            input_places=[self.places['E']],
            output_places=[self.places['I_s']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def e_to_ia_rate(input_places, output_places, params):
            E = input_places[0].get_tokens()
            sigma = params.get('sigma', self.p['sigma'])
            p_val = params.get('p', self.p['p'])
            return sigma * (1 - p_val) * E
        
        self.transitions['e_to_ia'] = PetriTransition(
            name='E_to_Ia',
            rate_function=e_to_ia_rate,
            input_places=[self.places['E']],
            output_places=[self.places['I_a']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def e_death_rate(input_places, output_places, params):
            E = input_places[0].get_tokens()
            mu = params.get('mu', self.p['mu'])
            return mu * E
        
        self.transitions['e_death'] = PetriTransition(
            name='E_Death',
            rate_function=e_death_rate,
            input_places=[self.places['E']],
            output_places=[],
            input_weights=[1.0],
            output_weights=[]
        )
        
        def is_to_c_rate(input_places, output_places, params):
            I_s = input_places[0].get_tokens()
            gamma = params.get('gamma', self.p['gamma'])
            q = params.get('q', self.p['q'])
            return gamma * q * I_s
        
        self.transitions['is_to_c'] = PetriTransition(
            name='Is_to_C',
            rate_function=is_to_c_rate,
            input_places=[self.places['I_s']],
            output_places=[self.places['C']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def is_to_r_rate(input_places, output_places, params):
            I_s = input_places[0].get_tokens()
            gamma = params.get('gamma', self.p['gamma'])
            q = params.get('q', self.p['q'])
            return gamma * (1 - q) * I_s
        
        self.transitions['is_to_r'] = PetriTransition(
            name='Is_to_R',
            rate_function=is_to_r_rate,
            input_places=[self.places['I_s']],
            output_places=[self.places['R']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def is_death_rate(input_places, output_places, params):
            I_s = input_places[0].get_tokens()
            mu = params.get('mu', self.p['mu'])
            return mu * I_s
        
        self.transitions['is_death'] = PetriTransition(
            name='Is_Death',
            rate_function=is_death_rate,
            input_places=[self.places['I_s']],
            output_places=[],
            input_weights=[1.0],
            output_weights=[]
        )
        
        def ia_to_r_rate(input_places, output_places, params):
            I_a = input_places[0].get_tokens()
            gamma = params.get('gamma', self.p['gamma'])
            return gamma * I_a
        
        self.transitions['ia_to_r'] = PetriTransition(
            name='Ia_to_R',
            rate_function=ia_to_r_rate,
            input_places=[self.places['I_a']],
            output_places=[self.places['R']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def ia_death_rate(input_places, output_places, params):
            I_a = input_places[0].get_tokens()
            mu = params.get('mu', self.p['mu'])
            return mu * I_a
        
        self.transitions['ia_death'] = PetriTransition(
            name='Ia_Death',
            rate_function=ia_death_rate,
            input_places=[self.places['I_a']],
            output_places=[],
            input_weights=[1.0],
            output_weights=[]
        )
        
        def c_to_r_rate(input_places, output_places, params):
            C = input_places[0].get_tokens()
            omega = params.get('omega', self.p['omega'])
            return omega * C
        
        self.transitions['c_to_r'] = PetriTransition(
            name='C_to_R',
            rate_function=c_to_r_rate,
            input_places=[self.places['C']],
            output_places=[self.places['R']],
            input_weights=[1.0],
            output_weights=[1.0]
        )
        
        def c_death_rate(input_places, output_places, params):
            C = input_places[0].get_tokens()
            mu = params.get('mu', self.p['mu'])
            return mu * C
        
        self.transitions['c_death'] = PetriTransition(
            name='C_Death',
            rate_function=c_death_rate,
            input_places=[self.places['C']],
            output_places=[],
            input_weights=[1.0],
            output_weights=[]
        )
        
        def r_death_rate(input_places, output_places, params):
            R = input_places[0].get_tokens()
            mu = params.get('mu', self.p['mu'])
            return mu * R
        
        self.transitions['r_death'] = PetriTransition(
            name='R_Death',
            rate_function=r_death_rate,
            input_places=[self.places['R']],
            output_places=[],
            input_weights=[1.0],
            output_weights=[]
        )
    
    def set_initial_conditions(self, initial_infected=1):
        N = self.p['N']
        
        for place in self.places.values():
            place.reset()
        
        self.places['S'].tokens = N - initial_infected
        self.places['E'].tokens = 0.0
        self.places['I_s'].tokens = initial_infected * self.p['p']
        self.places['I_a'].tokens = initial_infected * (1 - self.p['p'])
        self.places['C'].tokens = 0.0
        self.places['R'].tokens = 0.0
        
        for place in self.places.values():
            place.initial_tokens = place.tokens
    
    def get_state_vector(self) -> np.ndarray:
        return np.array([
            self.places['S'].get_tokens(),
            self.places['E'].get_tokens(),
            self.places['I_s'].get_tokens(),
            self.places['I_a'].get_tokens(),
            self.places['C'].get_tokens(),
            self.places['R'].get_tokens()
        ])
    
    def simulate(self, t_span: np.ndarray, dt: float = 0.01, **kwargs) -> np.ndarray:
        sim_params = self.p.copy()
        sim_params.update(kwargs)
        
        results = np.zeros((len(t_span), 6))
        results[0] = self.get_state_vector()
        
        current_time = t_span[0]
        result_idx = 1
        
        while result_idx < len(t_span) and current_time < t_span[-1]:
            for transition in self.transitions.values():
                if transition.is_enabled():
                    transition.fire(dt, sim_params)
            
            current_time += dt
            
            if result_idx < len(t_span) and current_time >= t_span[result_idx]:
                results[result_idx] = self.get_state_vector()
                result_idx += 1
        
        if result_idx < len(t_span):
            final_state = self.get_state_vector()
            for i in range(result_idx, len(t_span)):
                results[i] = final_state
        
        return results
    
    def calculate_daily_incidence(self, solution: np.ndarray) -> np.ndarray:
        E = solution[:, 1]
        daily_incidence = self.p['p'] * self.p['sigma'] * E
        return daily_incidence
    
    def extract_compartments(self, solution: np.ndarray) -> Dict[str, np.ndarray]:
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