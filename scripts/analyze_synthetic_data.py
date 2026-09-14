import types
from pathlib import Path

import numpy as np
import pandas as pd

from run_bootstrap import (
    METRIC_COLORS,
    PLOT_COLORS,
    TEST_LABELS,
    TEST_MARKERS,
    TESTS,
    configure_plot_font,
    do_tests,
    get_tests,
    prepare_value_scale,
    summarize_detection_rate_metrics,
    value_axis_label,
)


class SyntheticBlood:
    def __init__(self, crp_by_day, new_infections, random_seed=None, n_dayperson=None):
        self.lab_memory = {"crp": crp_by_day}
        self.sim = types.SimpleNamespace()
        self.sim.results = types.SimpleNamespace()
        self.sim.results.new_infections = np.asarray(new_infections)
        self.random_seed = random_seed
        self.n_dayperson = n_dayperson


DEFAULT_DAY_ZERO = pd.Timestamp("2000-01-01")
CRP_PLOT_COLOR = PLOT_COLORS["value_smoothing"]
INFECTIONS_PLOT_COLOR = PLOT_COLORS["secondary_axis"]


def load_synthetic_csv_bloods(
    n_dayperson,
    random_seeds=range(30),
    results_dir=Path("../synthetic_data"),
):
    results_dir = Path(results_dir)
    crp_path = results_dir / f"crp_{n_dayperson}.csv"
    infections_path = results_dir / f"new_infections_{n_dayperson}.csv"

    crp = pd.read_csv(crp_path, dtype={"day": int, "crp": float, "seed": int})
    infections = pd.read_csv(
        infections_path,
        dtype={"day": int, "new_infections": float, "seed": int},
    )
    random_seeds = list(random_seeds)
    crp = crp.loc[crp["seed"].isin(random_seeds)]
    infections = infections.loc[infections["seed"].isin(random_seeds)]

    if crp.empty:
        raise ValueError(f"No CRP rows found for seeds {random_seeds} in {crp_path}")

    max_day = int(max(crp["day"].max(), infections["day"].max()))
    bloods = []

    for seed in random_seeds:
        seed_crp = crp.loc[crp["seed"] == seed]
        seed_infections = infections.loc[infections["seed"] == seed]

        if seed_crp.empty or seed_infections.empty:
            continue

        crp_by_day = {
            int(day): values["crp"].to_numpy()
            for day, values in seed_crp.groupby("day", sort=True)
        }
        new_infections = (
            seed_infections
            .set_index("day")["new_infections"]
            .reindex(range(max_day + 1), fill_value=0)
            .to_numpy()
        )
        bloods.append(
            SyntheticBlood(
                crp_by_day=crp_by_day,
                new_infections=new_infections,
                random_seed=seed,
                n_dayperson=n_dayperson,
            )
        )

    if not bloods:
        raise ValueError(f"No synthetic runs loaded from {results_dir} for n_dayperson={n_dayperson}")

    return bloods


def date_to_sim_day(date, day_zero=DEFAULT_DAY_ZERO):
    return (pd.to_datetime(date) - day_zero).dt.days


def blood_crp_to_dataframe(blood, day_zero=DEFAULT_DAY_ZERO):
    rows = []

    for day, values in sorted(blood.lab_memory["crp"].items()):
        date = day_zero + pd.Timedelta(int(day), unit="D")
        for value in values:
            rows.append({"date": date, "value": value})

    return pd.DataFrame(rows)


def get_synthetic_tp_window(
    blood,
    day_zero=DEFAULT_DAY_ZERO,
    infection_threshold=100,
    label_delay_days=7,
):
    new_infections = np.asarray(blood.sim.results.new_infections)
    above_threshold = np.flatnonzero(new_infections > infection_threshold)

    if len(above_threshold) == 0:
        raise ValueError("No day with blood.sim.results.new_infections > infection_threshold")

    tp_start_day = int(above_threshold[0]) + label_delay_days
    tp_end_day = int(np.argmax(new_infections)) + label_delay_days

    if tp_end_day < tp_start_day:
        tp_start_day, tp_end_day = tp_end_day, tp_start_day

    return (
        day_zero + pd.Timedelta(tp_start_day, unit="D"),
        day_zero + pd.Timedelta(tp_end_day, unit="D"),
    )



def load_synthetic_bloods(
    n_dayperson,
    random_seeds=range(30),
    results_dir=Path("../synthetic_data"),
):
    return load_synthetic_csv_bloods(
        n_dayperson=n_dayperson,
        random_seeds=random_seeds,
        results_dir=results_dir,
    )



def analyze_synthetic_n_dayperson(
    n_dayperson,
    random_seeds=range(30),
    results_dir=Path("../synthetic_data"),
    metrics_output_path=None,
    rates_output_path=None,
    plot_output_path=None,
    window_size=14,
    alpha=0.05,
    chow_alpha=1e-5,
    led_min_aic_delta=0.0,
    align_plot_to_tp_start=False,
    log_values=True,
    log_offset=1.0,
    plot_end_day=None,
    chow_input="mean",
    led_input="mean",
    include_daily_detectors=True,
    tests_to_plot=None,
):
    if metrics_output_path is None:
        metrics_output_path = f"synthetic_detection_rate_metrics_{n_dayperson}.csv"
    if rates_output_path is None:
        rates_output_path = f"synthetic_detection_rates_{n_dayperson}.csv"
    if plot_output_path is None:
        plot_output_path = f"synthetic_detection_plot_{n_dayperson}.png"

    bloods = load_synthetic_bloods(
        n_dayperson=n_dayperson,
        random_seeds=random_seeds,
        results_dir=results_dir,
    )

    rates_df = synthetic_detection_rates(
        bloods,
        window_size=window_size,
        alpha=alpha,
        chow_alpha=chow_alpha,
        led_min_aic_delta=led_min_aic_delta,
        align_to_tp_start=align_plot_to_tp_start,
        log_values=log_values,
        log_offset=log_offset,
        chow_input=chow_input,
        led_input=led_input,
        include_daily_detectors=include_daily_detectors,
    )
    metrics_df = summarize_synthetic_detection_rates(
        rates_df,
        bloods,
        align_to_tp_start=align_plot_to_tp_start,
    )

    metrics_df.to_csv(metrics_output_path, index=False)
    rates_df.to_csv(rates_output_path, index=False)
    figure = plot_synthetic_detections(
        bloods,
        save_path=plot_output_path,
        window_size=window_size,
        alpha=alpha,
        chow_alpha=chow_alpha,
        led_min_aic_delta=led_min_aic_delta,
        align_to_tp_start=align_plot_to_tp_start,
        log_values=log_values,
        log_offset=log_offset,
        chow_input=chow_input,
        led_input=led_input,
        include_daily_detectors=include_daily_detectors,
        tests_to_plot=tests_to_plot,
        detection_rates=rates_df,
        end_day=plot_end_day,
    )

    return rates_df, metrics_df, figure



def summarize_synthetic_detection_rates(
    rates,
    bloods,
    day_zero=DEFAULT_DAY_ZERO,
    align_to_tp_start=False,
    smoothing_span=7,
):
    tp_windows = [get_synthetic_tp_window(blood, day_zero=day_zero) for blood in bloods]
    if align_to_tp_start:
        tp_start = day_zero
        tp_end = day_zero + max(tp_end - tp_start for tp_start, tp_end in tp_windows)
    else:
        tp_start = min(window[0] for window in tp_windows)
        tp_end = max(window[1] for window in tp_windows)

    return summarize_detection_rate_metrics(
        rates,
        tp_start=tp_start,
        tp_end=tp_end,
        reference_date=tp_start,
        smoothing_span=smoothing_span,
    )


def synthetic_detection_rates(
    bloods,
    day_zero=DEFAULT_DAY_ZERO,
    window_size=14,
    alpha=0.05,
    chow_alpha=1e-5,
    led_min_aic_delta=0.0,
    align_to_tp_start=False,
    log_values=True,
    log_offset=1.0,
    chow_input="mean",
    led_input="mean",
    include_daily_detectors=True,
):
    rates = []

    for blood_index, blood in enumerate(bloods):
        random_seed = getattr(blood, "random_seed", blood_index)
        data = prepare_value_scale(
            blood_crp_to_dataframe(blood, day_zero=day_zero),
            log_values=log_values,
            log_offset=log_offset,
        )
        start = data.date.min() + pd.Timedelta(3 * window_size, unit="D")
        end = data.date.max()
        results = do_tests(
            data,
            start=start,
            end=end,
            window_size=window_size,
            alpha=alpha,
            chow_alpha=chow_alpha,
            led_min_aic_delta=led_min_aic_delta,
            chow_input=chow_input,
            led_input=led_input,
            include_daily_detectors=include_daily_detectors,
        )
        if align_to_tp_start:
            tp_start, _ = get_synthetic_tp_window(blood, day_zero=day_zero)
            aligned_detection_date = (
                day_zero + (results["detection_date"] - tp_start)
            )
        else:
            aligned_detection_date = results["detection_date"]

        for test in get_tests(include_daily_detectors=include_daily_detectors):
            test_results = results.dropna(subset=[test])
            test_rates = pd.DataFrame(
                {
                    "date": aligned_detection_date.reindex(test_results.index).to_numpy(),
                    "detected": test_results[test].astype(bool).to_numpy(),
                }
            )
            test_rates["test"] = test
            test_rates["random_seed"] = random_seed
            rates.append(test_rates)

    rates_df = pd.concat(rates, ignore_index=True)
    return (
        rates_df
        .groupby(["date", "test"], as_index=False)
        .agg(detection_rate=("detected", "mean"), n_iterations=("detected", "size"))
    )


def plot_synthetic_detections(
    bloods,
    save_path=None,
    day_zero=DEFAULT_DAY_ZERO,
    window_size=14,
    alpha=0.05,
    chow_alpha=1e-5,
    led_min_aic_delta=0.0,
    align_to_tp_start=False,
    crp_seed_index=0,
    log_values=True,
    log_offset=1.0,
    smoothing_span=7,
    font_size=22,
    tests_to_plot=None,
    detection_rates=None,
    end_day=None,
    chow_input="mean",
    led_input="mean",
    include_daily_detectors=True,
):
    import matplotlib.pyplot as plt

    configure_plot_font("Inter")

    infection_rows = []
    tp_windows = []

    for blood_index, blood in enumerate(bloods):
        random_seed = getattr(blood, "random_seed", blood_index)
        tp_start, tp_end = get_synthetic_tp_window(blood, day_zero=day_zero)
        infections = np.asarray(blood.sim.results.new_infections)
        infections = infections.copy()
        if len(infections) > 100:
            infections[100] = 0
        infection_dates = day_zero + pd.to_timedelta(np.arange(len(infections)), unit="D")
        if align_to_tp_start:
            infection_dates = day_zero + (infection_dates - tp_start)
            tp_windows.append((day_zero, day_zero + (tp_end - tp_start)))
        else:
            tp_windows.append((tp_start, tp_end))

        infection_rows.append(pd.DataFrame({
            "date": infection_dates,
            "new_infections": infections,
            "random_seed": random_seed,
        }))

    tp_start = min(window[0] for window in tp_windows)
    tp_end = max(window[1] for window in tp_windows)
    infection_df = pd.concat(infection_rows, ignore_index=True)
    infection_summary = (
        infection_df
        .groupby("date", as_index=False)
        .agg(
            new_infections_mean=("new_infections", "mean"),
            new_infections_std=("new_infections", "std"),
        )
    )
    infection_summary["new_infections_std"] = infection_summary["new_infections_std"].fillna(0)

    crp_blood = bloods[crp_seed_index]
    crp_data = prepare_value_scale(
        blood_crp_to_dataframe(crp_blood, day_zero=day_zero),
        log_values=log_values,
        log_offset=log_offset,
    )
    if align_to_tp_start:
        crp_tp_start, _ = get_synthetic_tp_window(crp_blood, day_zero=day_zero)
        crp_data["date"] = day_zero + (crp_data["date"] - crp_tp_start)
    crp_daily_mean = (
        crp_data
        .groupby("date", as_index=False)["value"]
        .mean()
        .rename(columns={"value": "crp_mean"})
    )
    crp_daily_mean["crp_rolling_mean"] = (
        crp_daily_mean["crp_mean"]
        .rolling(window=window_size, min_periods=1)
        .mean()
    )

    if detection_rates is None:
        detection_rates = synthetic_detection_rates(
            bloods,
            day_zero=day_zero,
            window_size=window_size,
            alpha=alpha,
            chow_alpha=chow_alpha,
            led_min_aic_delta=led_min_aic_delta,
            chow_input=chow_input,
            led_input=led_input,
            include_daily_detectors=include_daily_detectors,
            align_to_tp_start=align_to_tp_start,
            log_values=log_values,
            log_offset=log_offset,
        )
    elif isinstance(detection_rates, (str, Path)):
        detection_rates = pd.read_csv(detection_rates, parse_dates=["date"])
    else:
        detection_rates = detection_rates.copy()
    infection_summary["sim_day"] = date_to_sim_day(infection_summary["date"], day_zero=day_zero)
    crp_daily_mean["sim_day"] = date_to_sim_day(crp_daily_mean["date"], day_zero=day_zero)
    detection_rates["sim_day"] = date_to_sim_day(detection_rates["date"], day_zero=day_zero)
    tests_to_plot = TESTS if tests_to_plot is None else list(tests_to_plot)
    detection_rates = detection_rates.loc[detection_rates["test"].isin(tests_to_plot)]

    start = detection_rates.sim_day.min()
    end = detection_rates.sim_day.max()
    if end_day is not None:
        end = min(end, end_day)
    tp_start_day = int((tp_start - day_zero).days)
    tp_end_day = int((tp_end - day_zero).days)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(12, 12),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1]},
    )
    ax_values, ax_detections = axes
    ax_infections = ax_values.twinx()

    for ax in (ax_values, ax_detections):
        ax.axvspan(start, tp_start_day, color=PLOT_COLORS["true_negative_fill"], alpha=0.35)
        ax.axvspan(tp_start_day, tp_end_day + 1, color=PLOT_COLORS["true_positive_fill"], alpha=0.28)
        ax.axvline(tp_start_day, color=PLOT_COLORS["outbreak"], linestyle="--", linewidth=1.2)

    ax_values.text(
        start + (tp_start_day - start) / 2,
        0.97,
        "True negative",
        transform=ax_values.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=font_size,
    )
    ax_values.text(
        tp_start_day + (tp_end_day + 1 - tp_start_day) / 2,
        0.97,
        "True positive",
        transform=ax_values.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=font_size,
    )

    ax_values.annotate(
        "Outbreak",
        xy=(tp_start_day, 0.85),
        xycoords=ax_values.get_xaxis_transform(),
        xytext=(-120, -50),
        textcoords="offset points",
        arrowprops={"arrowstyle": "->", "color": PLOT_COLORS["outbreak"], "linewidth": 1.2},
        fontsize=font_size,
    )

    infection_x = infection_summary["sim_day"]
    infection_mean = infection_summary["new_infections_mean"]
    infection_std = infection_summary["new_infections_std"]
    ax_infections.fill_between(
        infection_x,
        infection_mean - infection_std,
        infection_mean + infection_std,
        color=INFECTIONS_PLOT_COLOR,
        alpha=0.12,
        linewidth=0,
        # label="new infections std",
    )
    ax_infections.plot(
        infection_x,
        infection_mean,
        color=INFECTIONS_PLOT_COLOR,
        linewidth=2,
        label="New infections",
    )
    ax_values.scatter(
        crp_daily_mean["sim_day"],
        crp_daily_mean["crp_mean"],
        color=PLOT_COLORS["value_smoothing"],
        s=18,
        alpha=0.10,
    )
    ax_values.plot(
        crp_daily_mean["sim_day"],
        crp_daily_mean["crp_rolling_mean"],
        color=PLOT_COLORS["value_smoothing"],
        linewidth=2.2,
    )

    ax_values.set_ylabel(value_axis_label(log_values=log_values, log_offset=log_offset), fontsize=font_size)
    ax_infections.set_ylabel("New infections", color=INFECTIONS_PLOT_COLOR, fontsize=font_size)
    ax_values.tick_params(axis="both", labelsize=font_size)
    ax_infections.tick_params(axis="y", colors=INFECTIONS_PLOT_COLOR, labelsize=font_size)
    ax_infections.spines["right"].set_color(INFECTIONS_PLOT_COLOR)
    # ax_values.grid(alpha=0.25)

    value_lines, value_labels = ax_values.get_legend_handles_labels()
    infection_lines, infection_labels = ax_infections.get_legend_handles_labels()
    # ax_values.legend(value_lines + crp_lines, value_labels + crp_labels, loc="best", frameon=False)

    for test, test_rates in detection_rates.groupby("test"):
        test_label = TEST_LABELS.get(test, test)
        test_rates = test_rates.sort_values("sim_day").copy()
        smoothed_detection_rate = test_rates.detection_rate.rolling(
            window=smoothing_span,
            center=True,
            min_periods=smoothing_span,
        ).mean()
        ax_detections.scatter(
            test_rates.sim_day,
            test_rates.detection_rate,
            color=METRIC_COLORS.get(test_label),
            marker=TEST_MARKERS.get(test, "o"),
            s=24,
            alpha=0.45,
        )
        ax_detections.plot(
            test_rates.sim_day,
            smoothed_detection_rate,
            linewidth=1.5,
            color=METRIC_COLORS.get(test_label),
        )
        ax_detections.plot(
            [],
            [],
            linewidth=1.5,
            color=METRIC_COLORS.get(test_label),
            marker=TEST_MARKERS.get(test, "o"),
            markersize=7,
            label=test_label,
        )

    ax_detections.set_ylabel("Detection rate", fontsize=font_size)
    ax_detections.set_xlabel("Simulation day", fontsize=font_size)
    ax_detections.set_ylim(-0.02, 1.02)
    ax_detections.legend(loc="best", ncol=1, frameon=False, fontsize=font_size)
    ax_detections.tick_params(axis="both", labelsize=font_size)
    # ax_detections.grid(alpha=0.25)
    ax_detections.set_xlim(start - 2 * window_size, end)

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, (ax_values, ax_infections, ax_detections)


if __name__ == "__main__":
    for n_dayperson in [1250, 250, 50, 10]:
            
        log_values = True
        window_size = 14
        chow_alpha = 0.05

        analyze_synthetic_n_dayperson(
            random_seeds=range(1000),
            n_dayperson=n_dayperson,
            plot_output_path=f"../figures/synthetic_detection_plot_{n_dayperson}.png",
            metrics_output_path=f"../metrics/raw/synthetic_detection_rate_metrics_{n_dayperson}.csv",
            rates_output_path=f"../metrics/raw/synthetic_detection_rates_{n_dayperson}.csv",
            log_values=log_values,
            log_offset=0.1,
            window_size=window_size,
            chow_alpha=chow_alpha,
            plot_end_day=200,
            include_daily_detectors=False,
        )

        bloods = load_synthetic_bloods(n_dayperson, random_seeds=range(1000))

        plot_synthetic_detections(
            bloods,
            detection_rates=f"../metrics/raw/synthetic_detection_rates_{n_dayperson}.csv",
            save_path=f"../figures/synthetic_detection_plot_{n_dayperson}.png",
            end_day=200,
            chow_alpha=chow_alpha,
            window_size=window_size,
            log_values=log_values,
            log_offset=0.1,
            include_daily_detectors=False,
        )
