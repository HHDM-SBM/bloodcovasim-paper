import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def calculate_aic(n, rss, k):
    rss = max(float(rss), np.finfo(float).eps)
    return n * np.log(rss / n) + 2 * k

def random_noise_model(x, a):
    return a

def linear_growth_model(x, a, b):
    return a * x + b

def exponential_growth_model(x, a, b, c):
    return a * np.exp(b * x) + c


if __name__ == "__main__":
    np.random.seed(0)
    window_size = 28

    x = np.arange(window_size)
    y = np.exp(0.1 * (x+25)) + np.random.normal(0, 10, size=x.shape)

    fig, ax = plt.subplots()
    ax.plot(x, y, 'o', ms=5, color='black')
    ax.set_xlabel('Day', fontsize=16)
    ax.set_ylabel('CRP concentration [mg/L]', fontsize=16)

    
    popt_exponential, pcov_exponential = curve_fit(exponential_growth_model, x, y, p0=(1, 0.1, 0))
    
    popt_linear, pcov_linear = curve_fit(linear_growth_model, x, y)
    
    popt_random, pcov_random = curve_fit(random_noise_model, x, y)
    
    

    rss_random = np.sum((y - random_noise_model(x, *popt_random))**2)
    rss_linear = np.sum((y - linear_growth_model(x, *popt_linear))**2)
    rss_exponential = np.sum((y - exponential_growth_model(x, *popt_exponential))**2)

    aic_random = calculate_aic(window_size, rss_random, 1)
    aic_linear = calculate_aic(window_size, rss_linear, 2)
    aic_exponential = calculate_aic(window_size, rss_exponential, 2)

    print(aic_random, aic_linear, aic_exponential)

    ax.plot(x, exponential_growth_model(x, *popt_exponential), '-', label=f'AIC = {aic_exponential:.0f}', color='#D55E00')
    ax.plot(x, linear_growth_model(x, *popt_linear), '-', label=f'AIC = {aic_linear:.0f}', color='#F0E442')
    ax.plot(x, np.ones_like(x)*random_noise_model(x, *popt_random), '-', label=f'AIC = {aic_random:.0f}', color='#009E73')

    ax.legend(loc='upper center', fontsize=16, frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig('led_graph.png')