import pandas as pd
import numpy as np
import glob
import os
import random
import covasim as cv
# from profilehooks import profile
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d
import cProfile


def profile(func):
    """Decorator for run function profile"""
    def wrapper(*args, **kwargs):
        profile_filename = func.__name__ + '.prof'
        profiler = cProfile.Profile()
        result = profiler.runcall(func, *args, **kwargs)
        profiler.dump_stats(profile_filename)
        return result
    return wrapper

class Person:
    """
    У каждого агента есть параметры:
    id
    sex - пол (не задано)
    age - возраст (не задано)
    state - состояние агента
    new_state_day - день, в который агент перешел в настоящее состояни
    blood_parameters - словарь с датафреймами изменения показателей крови на каждый день (наверное, надо сделать один общий?)
    personal_dynamic - словрь с ежедневной динамикой каждого параметра крови
    """

    def __init__(self, n_day, person_id, blood_parameters, person_states):
        
        self.n_day = n_day
        self.person_id = str(person_id)
        self.blood_parameters = blood_parameters
        self.state = 'sus'
        
        #        self.sex = sex
        #        self.age = age
        self.new_state_day = 0
        self.normal_values = {}
        #self.normal_pop_values = {'crp': [2, 1, 5],
        #                          'gem': [3, 1, 160]}
        self.normal_pop_values = {'crp': { 'normal': [0, 1.5,5], 
                                           'obesity': [2.02, 0.28, 10], 
                                           'derma': [2.65,0.35,25], 
                                           'heart': [4.25,1.25,200]}}
        self.other_impact = ''
        c_1 = person_states
        self.change_state_dates = [int(x) if not np.isnan(x) else np.nan for x in c_1]
        self.init_person_blood()

    def init_person_blood(self):
        '''
        создание словаря параметров крови конкретной персоны
        '''
        self.personal_dynamic = {}
        for blood_parameter in self.blood_parameters.keys():
            self.normal_values[blood_parameter] = int(self.get_normal_value(blood_parameter))
            self.personal_dynamic[blood_parameter] =  self.generate_dynamic(blood_parameter, normal_value = self.normal_values[blood_parameter], app_time = 0 , max_day_after_symp = 3)
        #     print('1) !!! self.personal_dynamic[blood_parameter]', self.personal_dynamic[blood_parameter])

    def get_normal_value(self, blood_parameter, porog=False):
        other_impact = np.random.choice(['normal', 'obesity', 'derma', 'heart'], 1, p = [0.7, 0.15, 0.1, 0.05])
        self.other_impact = other_impact[0]
        while porog != True:
            normal_values = np.random.normal(self.normal_pop_values[blood_parameter][self.other_impact][0],
                                             self.normal_pop_values[blood_parameter][self.other_impact][1], 1)
            # Преобразование значений в логнормальные
            lognormal_value = np.exp(normal_values)[0]
            if lognormal_value <= self.normal_pop_values[blood_parameter][self.other_impact][2]:
                porog = True
        return lognormal_value

    def _update_state(self, time):
        index = self.change_state_dates.index(time) + 1  # получение индекса необходимого значения
        self.new_state_time = time
        state = ['person_id', 'exposed', 'infectious', 'symptomatic', 'severe', 'critical', 'recovered', 'death'][index]  
        return state  # возврат обновленного состояния


    def check_state(self, time):
        if time in self.change_state_dates:
        #     print('check state')
            self.state = self._update_state(time)

    def generate_max_state(self):
        names = ['exposed',	'infectious',	'symptomatic',	'severe',	'critical',	'recovered']
        max_value = np.nanmax(self.change_state_dates[:5])
        max_index = np.where(self.change_state_dates[:5] == max_value)[0]
        i = max_index[-1] if max_index.size > 0 else None
        return names[i]

    def generate_max_value(self, max_state):
        if max_state == 'symptomatic':
            max_value = random.randint(10, 60)
        elif max_state == 'severe' or max_state == 'critical':
            max_value = random.randint(60, 140)
        else:
            max_value = random.randint(2, 10)
        return max_value

    def generate_dynamic(self, blood_parameter, normal_value, app_time, max_day_after_symp):
        dyn_str = []
        if np.all(np.isnan(self.change_state_dates)):


            for i in range(self.n_day):
                dyn_str.append(normal_value)
            return dyn_str
        
        max_state = self.generate_max_state()
        max_value = normal_value + self.generate_max_value( max_state)
        
        for i in range(self.change_state_dates[1]):            #   EXPOSED  changed
            dyn_str.append(normal_value)        #   EXPOSED
        
        if max_state == 'critical' or max_state == 'severe' or max_state == 'symptomatic' or max_state == 'infectious':
            for i in range(app_time):
                dyn_str.append(normal_value)

            if max_state != 'infectious':
                
                fast_growth_days = self.change_state_dates[2] + max_day_after_symp - len(dyn_str) # до какого дня показатель резко растет
                prob= list(np.linspace(normal_value, max_value, fast_growth_days))
                dyn_str += prob

            if not np.isnan(self.change_state_dates[6]):
                
                dyn_str += list(np.linspace(max_value, max_value+30, self.change_state_dates[6]-len(dyn_str)))
                for i in range(self.n_day-len(dyn_str)):
                    dyn_str.append(np.nan)
                if len(dyn_str) > self.n_day:
                    #print("AAAAAAAAAAAAAAAAAA", self.person_id, len(dyn_str), self.n_day)
                    dyn_str = dyn_str[:self.n_day]
                    
                return dyn_str
                 

            second_change_days = (self.change_state_dates[5] - len(dyn_str))//2
            
            if max_state == 'symptomatic' or max_state == 'infectious':
                
                for i in range(second_change_days):
                    dyn_str.append(max_value)
            else:
                
                dyn_str += list(np.linspace(max_value, max_value+15, second_change_days))
                
                    

            dyn_str += list(np.linspace(dyn_str[-1], normal_value, second_change_days+1))    

        for i in range(self.n_day-len(dyn_str)):
            dyn_str.append(normal_value)

        dyn_str = gaussian_filter1d(dyn_str, sigma=2)
        if len(dyn_str) > self.n_day:
            #print("AAAAAAAAAAAAAAAAAA", self.person_id, len(dyn_str), self.n_day)
            dyn_str = dyn_str[:self.n_day]
        return dyn_str
    
    def live_day(self, time):

        self.check_state(time)
        
        
        


#        self.update_blood(time)


class BloodSim():
    """
    общие функции и функции, объединяющие предыдущие классы

    функции класса:

    sim_run(self) - запуск симуляции. Если первый день - инициализация агентов и их параметров. Затем ежедневное обновление параметров агентов
    init_population(self) - создание массива агентов
    check_parameters(self) - вывести рассматриваемые параметры
    do_df_pop_parameters(self) - объединение показателей крови всей популяции в словарь датафреймов (self.pop_blood)
    preprocessing_blood_files(self) - создание словаря со значениями динамики параметров крови при различных состояниях

    start_day - первый день симуляции
    end_day - последний день симуляции
    time - нынешний день симуляции
    pop_size - размер популяции
    population - массив агентов (np.array)
    """

    def __init__(self, n_person_per_day, start_day=0, rand_seed = 0, end_day=300, pop_size=10000, guest_strategy='random', variant = cv.variant('alpha', days=100, n_imports=30)):
        self.start_day = start_day
        self.end_day = end_day
        self.time = start_day  # день симуляции
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
        self.variant = variant
        self.random_seed = rand_seed
        self.sim = self.do_covasim()

    def do_covasim(self):
        sim = cv.Sim(variants = self.variant, n_days = self.end_day - self.start_day, pop_size = self.pop_size, rand_seed = self.random_seed, pop_infected=0)
        sim.run()
        print(sim)
        return sim


    def preprocessing_blood_files(self):

        '''
        Создание словаря со значениями динамики параметров крови
        '''
        self.blood_parameters = {}  # словарь с дф параметров крови
        for file in glob.glob(r'blood/*.xlsx'):
            filename = os.path.split(file)[1]
            wo_ext = os.path.splitext(filename)[0]
            self.blood_parameters[wo_ext] = pd.read_excel(file)

    def init_population(self):
        '''
        инициализация агентов
        '''
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

    #   @profile(stdout=False, filename='profile_file.prof')
    @profile
    def sim_run(self):

        if self.time == self.start_day:
            self.init_population()

        for day in range(self.end_day):
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
                elif person.state == 'death':
                    self.death_persons.append(int(person.person_id))
                else:
                    self.labor_normal.append(int(person.person_id))
                
                if person.state != 'death': 
                    self.labor_for_random.append(int(person.person_id))

            self.labor_normal_states_dict[day] = self.labor_normal
            self.labor_crit_sev_dict[day] = self.labor_crit_sev
            self.labor_symp_dict[day] = self.labor_symp
            self.labor_for_random_dict[day] = self.labor_for_random
            self.death_persons_dict[day] = self.death_persons

        self.do_df_pop_parameters()
        self.do_guests()

    def check_parameters(self):
        '''
        напечатать имена параметров крови
        '''
        self.preprocessing_blood_files()
        return list(self.blood_parameters.keys())

    def do_df_pop_parameters(self):
        self.pop_blood = {}
    
    # Генерируем список имен столбцов
        col_names = ['person_id'] + [f'day_{day}' for day in range(self.start_day, self.end_day)]

        for blood_parameter in self.check_parameters():
        # Подготавливаем список строк
            data = []

            for person in self.population:
                new_row = [str(person.person_id)] + list(person.personal_dynamic[blood_parameter])
                data.append(new_row)

        # Создаем DataFrame из списка данных с указанием столбцов
            df = pd.DataFrame(data, columns=col_names)
        
            self.pop_blood[blood_parameter] = df

    def do_guests(self):
        '''
        Выбор стратегии тестирования и создание списка посещающих
        '''

        if self.guest_strategy == 'random':
            self.add_random_id_list()
        elif self.guest_strategy == 'symptomatic':
            self.add_symp_id_list()

        self.lab_memory = self.get_lab_results()
        return

    def add_symp_id_list(self):
        '''
        * сделать три
        '''
        
        self.n_person_per_day['people_ids'] = self.n_person_per_day.apply(self.generate_numbers_with_coef, axis=1)
        
        pass

    def generate_numbers_with_coef(self, row):
        '''
        Фунция на замену. Недостаток: выбирает с заданной вероятность СПИСКИ людей (в результате если болеет только один человек, а пришли 5 - он почти однозначно пройдет проверку, даже если население 1000+)
        Предполагаемое решение: проводит пересчет веса в зависимости от длин списков людей с симптомами/в крит состоянии
        Выбрать определенный
        '''
        
        count = int(row['n_person'])
        day = int(row['day'])
        names = ['  crit_sev  ', '  symp  ', '  normal  ']
        lists = [self.labor_crit_sev_dict[day], self.labor_symp_dict[day],
                 self.labor_normal_states_dict[day]]  # убирать с этой строки
        weights = [0.35, 0.34, 0.32]
        result = []
        people_group_count = self.get_people_group_count(lists, weights, count)

        for group in range(3):
            random_values = random.sample(lists[group], people_group_count[group])
        #    print('DAY:', day, names[group], '  count:    ', people_group_count[group], lists[group])
            result += random_values
        #print('RESULT:  ' , result)

        '''
        non_empty_lists = [lst for lst in lists if lst]
        non_empty_weights = [weights[i] for i, lst in enumerate(lists) if lst]


        for i in range(count):

            non_empty_lists = [lst for lst in non_empty_lists if lst]
            non_empty_weights = [weights[i] for i, lst in enumerate(non_empty_lists) if lst]

            if not non_empty_lists:
                raise ValueError("Все списки пусты. Невозможно сделать выбор.")

            chosen_list_index = random.choices(range(len(non_empty_lists)), weights=non_empty_weights, k=1)[0]
            chosen_number = random.choice(non_empty_lists[chosen_list_index])

            result.append(chosen_number)
            non_empty_lists[chosen_list_index].remove(chosen_number)
        '''
        return result

    def get_people_group_count(self, lists, weights, count):

        crit_count = int(count * weights[0])
        if crit_count > len(lists[0]):
            crit_count = len(lists[0])

        symp_count = int(count * weights[1] / (weights[1] + weights[2]))
        if symp_count > len(lists[1]):
            symp_count = len(lists[1])

        normal_count = count - symp_count - crit_count
        result_list = [crit_count, symp_count, normal_count]

        return result_list

    def generate_random_numbers(self, row):
        size = int(row['n_person'])
        day = int(row['day'])
        self.to_see = random.sample(self.labor_for_random_dict[day], size)
        
#        self.to_see = np.random.randint(self.pop_size, size=size).tolist()
        return self.to_see

    def add_random_id_list(self):
        '''
        добавление в таблицу посетителей айдишники посетивших на каждый день
        '''
        self.n_person_per_day['people_ids'] = self.n_person_per_day.apply(self.generate_random_numbers, axis=1)

    def get_lab_results(self):
        '''
        получение результатов лабораторной диагностики посетителей клиники
        '''
        lab_results = {}
        for blood_parameter in self.pop_blood.keys():
            param_days_dict = {}
            for day in range(self.end_day - self.start_day):
                df = self.choose_df_parts(blood_parameter, day)
                id_list = list(self.n_person_per_day.people_ids[self.n_person_per_day.day == day])[0]
                prom = df[df['person_id'].astype(int).isin(id_list)]
                param_days_dict[day] = list(prom.iloc[:, 1])
            lab_results[blood_parameter] = param_days_dict
        return lab_results

    def choose_df_parts(self, blood_parameter, day):
        df = self.pop_blood[blood_parameter]
        matching_columns = [col for col in df.columns if col.split('_')[-1] == str(day)]
        first_column = df.iloc[:, 0]
        result_df = df[matching_columns]
        result_df.insert(0, 'person_id', first_column)
        return result_df


def replace_zero_to_nan(df):
    '''
    Строки, где все ячейки (кроме id) принимают значение 0 или nan - заменить на nan
    '''
    first_column = df.iloc[:, 0]
    df = df.iloc[:, 1:]
    mask = (df == 0) | (df.isna())
    rows_to_replace = mask.all(axis=1)
    df.loc[rows_to_replace] = df.loc[rows_to_replace].applymap(lambda x: np.nan if not isinstance(x, str) else x)
    df.insert(0, 'person_id', first_column)
    df.person_id = df.person_id.astype(str)
    return df




def do_plot(df_dict):
    '''
    Построение графиков параметров крови по словарю датафреймов.
    '''
    for key_word in df_dict.keys():
        df = df_dict[key_word]
        mean_values = df.iloc[:, 1:].mean()
        plt.plot(mean_values, marker='o', label=key_word, alpha=0.5)

    nbis_value = df.shape[1]
    if nbis_value > 10:
        nbis_value = 10

    plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=nbis_value))
    plt.legend()
    plt.show()
