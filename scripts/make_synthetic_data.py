from bloodsim import BloodSim
import os
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor


def make_input_file(n_days = 300, n_dayperson = 250):
    n_dayperson_dict = {
        'day': list(range(n_days)),
        'n_person': [n_dayperson] * n_days,
    }
    
    return pd.DataFrame(n_dayperson_dict)


def run_seed(args):
    i, n_dayperson = args
    pop_size = int(1e6)
    n_days = 300

    blood = BloodSim(n_person_per_day = make_input_file(n_days = n_days, n_dayperson = n_dayperson), end_day = n_days, pop_size = pop_size, guest_strategy = 'random', rand_seed = i)
    blood.sim_run()

    crp = pd.DataFrame(blood.lab_memory['crp']).melt(var_name='day', value_name='crp')
    crp['seed'] = i

    new_infections = pd.DataFrame({
        'day': np.arange(len(blood.sim.results.new_infections)),
        'new_infections': np.array(blood.sim.results.new_infections),
        'seed': i,
    })

    return crp, new_infections


if __name__=='__main__':
    seeds = 1000
    workers = 25
    n_dayperson = 1250
    
    os.makedirs('synthetic_data', exist_ok=True)

    with ProcessPoolExecutor(max_workers = workers) as executor:
        results = list(executor.map(run_seed, [(i, n_dayperson) for i in range(seeds)]))

    crp = pd.concat([x[0] for x in results], ignore_index=True)
    new_infections = pd.concat([x[1] for x in results], ignore_index=True)

    crp.to_csv(f'synthetic_data/crp_{n_dayperson}.csv', index=False)
    new_infections.to_csv(f'synthetic_data/new_infections_{n_dayperson}.csv', index=False)
