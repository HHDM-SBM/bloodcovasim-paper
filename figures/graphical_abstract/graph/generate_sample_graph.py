import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["svg.fonttype"] = "none"


if __name__ == "__main__":

    np.random.seed(0)


    x = np.arange(115)
    y = np.ones_like(x) * 0.5 + np.random.normal(0, 1, size=x.shape) + np.exp(0.025 * (x-64))
    y_mean = [np.mean(y[i-14:i]) for i in range(14, len(y))]


    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(x, y, marker="o", markersize=5, linestyle="", alpha=0.5, color="black") 
    ax.plot(x[14:], y_mean, color="black", linewidth=2)
    ax.axvspan(xmin=50, xmax=78, alpha=0.3, color="black")
    ax.set_xlim(14, 114)
    ax.set_xlabel("Day", fontsize=16)
    ax.set_ylabel("CRP concentration [mg/L]", fontsize=16)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig("sample_graph.svg", format="svg")