from bloodsim import BloodSim
import pickle
import types
import numpy as np
import pandas as pd


class Blood:
    def __init__(self, n_day, lab_memory, sim):
        self.n_day = n_day
        self.lab_memory = {'crp': lab_memory, 'crp_mean': np.array(blood.pop_blood['crp'].iloc[:, 1:].mean())}
        self.sim = types.SimpleNamespace()
        self.sim.results = types.SimpleNamespace()
        self.sim.results.new_infections = np.array(sim.results.new_infections)
        pass


def make_input_file(n_days = 300, n_dayperson = 250):
    n_dayperson_dict = {
        'day': list(range(n_days)),
        'n_person': [n_dayperson] * n_days,
    }
    
    return pd.DataFrame(n_dayperson_dict)


if __name__=='__main__':
    pop_size = int(1e6)
    n_dayperson = 1250
    seeds = 1
    
    for i in range(seeds):
        blood = BloodSim(n_person_per_day = make_input_file(n_days = 300, n_dayperson = n_dayperson), end_day = 300, pop_size = pop_size, guest_strategy = 'random', rand_seed = i)
        blood.sim_run()
        lossy_blood = Blood(n_day = len(blood.lab_memory['crp'].keys()), lab_memory = blood.lab_memory['crp'], sim = blood.sim)

        with open(f'synthetic_data/synthetic_blood_{n_dayperson}_{i}.pkl', 'wb') as f:
            pickle.dump(lossy_blood, f)
