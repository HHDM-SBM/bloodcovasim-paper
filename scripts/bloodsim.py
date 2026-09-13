import pandas as pd
import numpy as np
import glob
import os
import random
import warnings
import covasim as cv
from scipy.ndimage.filters import gaussian_filter1d



class Person:
    ''' The Person class handles the individual dynamics of CRP concentration in the blood plasma of a specific person based on the progression of their disease as recorded in Covasim. '''

    def __init__(self, n_day, person_id, blood_parameters, person_states):
        self.n_day = n_day
        self.person_id = str(person_id)
        self.blood_parameters = blood_parameters
        self.state = 'sus'
        self.new_state_day = 0
        self.normal_values = {}
        c_1 = person_states
        self.change_state_dates = [int(x) if not np.isnan(x) else np.nan for x in c_1]
        self.init_person_blood()


    def init_person_blood(self):
        self.personal_dynamic = {}
        for blood_parameter in self.blood_parameters.keys():
            self.normal_values[blood_parameter] = np.random.lognormal(np.log(4), 0.8)
            self.personal_dynamic[blood_parameter] =  self.generate_dynamics(normal_value = self.normal_values[blood_parameter])


    def _update_state(self, time):
        index = max(i for i, date in enumerate(self.change_state_dates) if date == time) + 1
        self.new_state_time = time
        state = ['person_id', 'exposed', 'infectious', 'symptomatic', 'severe', 'critical', 'recovered', 'dead'][index]  
        return state


    def check_state(self, time):
        if time in self.change_state_dates:
            self.state = self._update_state(time)


    def generate_max_state(self):
        names = ['exposed',	'infectious', 'symptomatic', 'severe', 'critical', 'dead']
        
        dates = self.change_state_dates[:5] + self.change_state_dates[6:]
        max_value = np.nanmax(dates)
        max_index = np.where(dates == max_value)[0]
        i = max_index[-1] if max_index.size > 0 else None
        return names[i]


    def generate_max_value(self, max_state):
        match max_state:
            case "symptomatic":
                mean, std = np.log(12), 0.5
            case "severe":
                mean, std = np.log(50), 0.3
            case "critical":
                mean, std = np.log(150), 0.15
            case "dead":
                mean, std = np.log(250), 0.1

        return np.random.lognormal(mean, std)
    
    
    def generate_dynamics(self, normal_value):
        names = ['exposed',	'infectious', 'symptomatic', 'severe', 'critical', 'recovered', 'dead']

        dyn_str = []
        
        if np.all(np.isnan(self.change_state_dates)):
            return np.array([normal_value] * self.n_day)
        
        max_state = self.generate_max_state()

        if max_state == 'exposed' or max_state == 'infectious':
            return np.array([normal_value] * self.n_day)

        max_value = self.generate_max_value(max_state)
        peak_value = normal_value + max_value

        dyn_str += [normal_value] * self.change_state_dates[2]

        growth_days = self.change_state_dates[names.index(max_state)] + 1 - len(dyn_str)
        growth_days = max(1, growth_days)
        prob = np.linspace(normal_value, peak_value, growth_days).tolist()
        dyn_str += prob

        if max_state == 'dead':
            dyn_str = gaussian_filter1d(dyn_str, sigma=2).tolist()
            dyn_str += [np.nan] * (self.n_day - len(dyn_str))
        else:
            if np.isnan(self.change_state_dates[5]):
                dyn_str += [peak_value] * (self.n_day - len(dyn_str))
            else:
                decrease_days = self.change_state_dates[5] + 2 - len(dyn_str)
                decrease_days = max(1, decrease_days)
                dyn_str += np.linspace(peak_value, normal_value, decrease_days).tolist()[1:]
                dyn_str += [normal_value] * (self.n_day - len(dyn_str))
            dyn_str = gaussian_filter1d(dyn_str, sigma=2)

        if len(dyn_str) > self.n_day:
            dyn_str = dyn_str[:self.n_day]

        return np.array(dyn_str)
    

    def live_day(self, time):
        self.check_state(time)
        
        
        


class BloodSim():
    ''' The BloodSim class handles the process of conducting blood tests in the population following a Covasim simulation. '''

    def __init__(self, n_person_per_day, start_day=0, rand_seed=0, end_day=300, pop_size=10000, guest_strategy='random', variant=None, use_waning=False):
        if start_day != 0:
            raise ValueError('BloodSim currently supports only start_day = 0')

        self.start_day = start_day
        self.end_day = end_day
        self.time = start_day
        self.pop_size = pop_size
        self.guest_strategy = guest_strategy
        self.population = np.empty(pop_size, dtype=object)
        self.preprocessing_blood_files()
        self.n_person_per_day = n_person_per_day
        self.n_day = end_day-start_day
        self.labor_normal_states_dict = {}
        self.labor_crit_sev_dict = {}
        self.labor_symp_dict = {}
        self.labor_for_random_dict = {}
        self.death_persons_dict = {}
        if variant is None:
            variant = cv.variant('alpha', days=100, n_imports=30)
        self.variant = variant
        self.random_seed = rand_seed
        self.use_waning = use_waning
        self.sim = self.do_covasim()


    def do_covasim(self):
        sim = cv.Sim(variants = self.variant, n_days = self.end_day - self.start_day, pop_size = self.pop_size, rand_seed = self.random_seed, pop_infected=0, rescale=False, use_waning=self.use_waning)
        sim.run()
        print(sim)
        return sim


    def preprocessing_blood_files(self):
        self.blood_parameters = {}
        for file in glob.glob(r'blood/*.xlsx'):
            filename = os.path.split(file)[1]
            wo_ext = os.path.splitext(filename)[0]
            self.blood_parameters[wo_ext] = pd.read_excel(file)

        if len(self.blood_parameters) == 0:
            self.blood_parameters['crp'] = pd.DataFrame()


    def init_population(self):
        recovered = self.sim.people.date_recovered
        exposed = self.sim.people.date_exposed
        death = self.sim.people.date_dead
        infectious = self.sim.people.date_infectious
        severe = self.sim.people.date_severe
        critical = self.sim.people.date_critical
        symptomatic = self.sim.people.date_symptomatic
        
        for person_id in range(self.pop_size):
            
            person_states = [exposed[person_id], infectious[person_id] , symptomatic[person_id], severe[person_id], critical[person_id], recovered[person_id], death[person_id]]
            person = Person(n_day = self.n_day, person_id = person_id, blood_parameters = self.blood_parameters, person_states = person_states)
            
            self.population[person_id] = person


    def update_time(self):
        self.time += 1


    def sim_run(self):
        if self.time == self.start_day:
            self.init_population()

        for day in range(self.start_day, self.end_day):
            self.update_time()
            self.labor_normal = []
            self.labor_symp = []
            self.labor_crit_sev = []
            self.labor_for_random = []
            self.death_persons = []

            for person in self.population:
                person.live_day(time=day)
                if person.state == 'symptomatic':
                    self.labor_symp.append(int(person.person_id))
                elif person.state == 'severe' or person.state == "critical":
                    self.labor_crit_sev.append(int(person.person_id))
                elif person.state == 'dead':
                    self.death_persons.append(int(person.person_id))
                else:
                    self.labor_normal.append(int(person.person_id))
                
                if person.state != 'dead': 
                    self.labor_for_random.append(int(person.person_id))

            self.labor_normal_states_dict[day] = self.labor_normal
            self.labor_crit_sev_dict[day] = self.labor_crit_sev
            self.labor_symp_dict[day] = self.labor_symp
            self.labor_for_random_dict[day] = self.labor_for_random
            self.death_persons_dict[day] = self.death_persons

        self.do_df_pop_parameters()
        self.do_guests()


    def check_parameters(self):
        return list(self.blood_parameters.keys())


    def do_df_pop_parameters(self):
        self.pop_blood = {}
    
        col_names = ['person_id'] + [f'day_{day}' for day in range(self.start_day, self.end_day)]

        for blood_parameter in self.check_parameters():
            data = []

            for person in self.population:
                new_row = [str(person.person_id)] + list(person.personal_dynamic[blood_parameter])
                data.append(new_row)

            df = pd.DataFrame(data, columns=col_names)
        
            self.pop_blood[blood_parameter] = df


    def do_guests(self):
        if self.guest_strategy == 'random':
            self.add_random_id_list()
        else:
            raise ValueError("BloodSim currently supports only guest_strategy = 'random'")

        self.lab_memory = self.get_lab_results()
        return


    def generate_random_numbers(self, row):
        size = int(row['n_person'])
        day = int(row['day'])
        if day not in self.labor_for_random_dict:
            raise ValueError(f'Day {day} is outside simulated range')
        if size > len(self.labor_for_random_dict[day]):
            warnings.warn(f'Not enough people for random sample on day {day}: requested {size}, available {len(self.labor_for_random_dict[day])}')
            size = len(self.labor_for_random_dict[day])
        self.to_see = random.sample(self.labor_for_random_dict[day], size)
        return self.to_see


    def add_random_id_list(self):
        self.n_person_per_day['people_ids'] = self.n_person_per_day.apply(self.generate_random_numbers, axis=1)


    def get_lab_results(self):
        lab_results = {}
        for blood_parameter in self.pop_blood.keys():
            param_days_dict = {}
            for day in range(self.start_day, self.end_day):
                df = self.choose_df_parts(blood_parameter, day)
                day_people_ids = self.n_person_per_day.people_ids[self.n_person_per_day.day == day]
                if len(day_people_ids) != 1:
                    raise ValueError(f'Expected one people_ids row for day {day}, found {len(day_people_ids)}')
                id_list = day_people_ids.iloc[0]
                prom = df[df['person_id'].astype(int).isin(id_list)]
                param_days_dict[day] = list(prom.iloc[:, 1])
            lab_results[blood_parameter] = param_days_dict
        return lab_results


    def choose_df_parts(self, blood_parameter, day):
        df = self.pop_blood[blood_parameter]
        matching_columns = [f'day_{day}'] if f'day_{day}' in df.columns else []
        if len(matching_columns) != 1:
            raise KeyError(f'Expected one blood column for day {day}, found {len(matching_columns)}')
        first_column = df.iloc[:, 0]
        result_df = df[matching_columns]
        result_df.insert(0, 'person_id', first_column)
        return result_df


def replace_zero_to_nan(df):
    first_column = df.iloc[:, 0]
    df = df.iloc[:, 1:]
    mask = (df == 0) | (df.isna())
    rows_to_replace = mask.all(axis=1)
    df.loc[rows_to_replace] = df.loc[rows_to_replace].applymap(lambda x: np.nan if not isinstance(x, str) else x)
    df.insert(0, 'person_id', first_column)
    df.person_id = df.person_id.astype(str)
    return df
