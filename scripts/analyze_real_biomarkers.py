from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import anderson_ksamp

from run_bootstrap import (
    LED_aic_summary,
    PLOT_COLORS,
    anchored_weekly_ticks,
    configure_plot_font,
    get_bootstrap_sample,
    get_randomized_period_sample,
    load_real_data,
    prepare_value_scale,
    summarize_detection_rate_metrics,
)


BIOMARKER_LABELS = {
    8: "CRP",
    1001: "SIRI",
    1002: "SII",
    1003: "AISI",
    1009: "NLR",
    1010: "PLR",
    1011: "NMR",
    1012: "LMR",
}

BIOMARKER_COLORS = {
    8: "#CC79A7",
    1001: "#E69F00",
    1002: "#000000",
    1003: "#009E73",
    1009: "#F0E442",
    1010: "#0072B2",
    1011: "#D55E00",
    1012: "#CC79A7",
}

BIOMARKER_MARKERS = {
    8: "*",
    1001: "o",
    1002: "s",
    1003: "^",
    1009: "D",
    1010: "v",
    1011: "P",
    1012: "X",
}


def load_biomarkers_data(
    path=Path("../data/biomarkers_real_data.tsv"),
    crp_path=Path("../data/crp_real_data.tsv"),
    include_crp=True,
):
    data = pd.read_csv(path, sep="\t")
    required = {"date", "parameter_id", "value"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    data = data.loc[data["parameter_id"].isin(BIOMARKER_LABELS)].copy()
    data["date"] = pd.to_datetime(data["date"])
    data["parameter_id"] = data["parameter_id"].astype(int)
    data["value"] = pd.to_numeric(data["value"], errors="raise")
    if include_crp:
        crp = load_real_data(
            crp_path,
            start=data["date"].min(),
            end=data["date"].max(),
        )
        crp["parameter_id"] = 8
        data = pd.concat([data, crp], ignore_index=True)
    return data


def detect_led_and_ad(
    data,
    start=pd.Timestamp("2025-07-01"),
    end=pd.Timestamp("2025-10-01"),
    window_size=14,
    alpha=0.05,
    led_min_aic_delta=0.0,
):
    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    delay = pd.Timedelta(window_size, unit="D")
    candidate_dates = pd.date_range(start - delay, end - delay, freq="D")
    rolling_dates = pd.date_range(candidate_dates.min() - delay, candidate_dates.max() + delay, freq="D")

    rolling_mean = pd.Series(index=rolling_dates, dtype=float)
    for day in rolling_dates:
        in_window = data["date"].between(day - delay, day, inclusive="right")
        rolling_mean.loc[day] = data.loc[in_window, "value"].mean()

    rows = []
    for day in candidate_dates:
        before = data.loc[data["date"].between(day - delay, day, inclusive="left"), "value"]
        after = data.loc[data["date"].between(day, day + delay, inclusive="left"), "value"]
        raw_increase = after.mean() > before.mean()
        if before.empty or after.empty:
            ad_detected = False
        else:
            ad_pvalue = anderson_ksamp([before, after]).significance_level
            ad_detected = bool((ad_pvalue < alpha) and raw_increase)

        led_window = rolling_mean.loc[day - delay:day + delay - pd.Timedelta(days=1)]
        led_before = rolling_mean.loc[day - delay:day - pd.Timedelta(days=1)]
        led_after = rolling_mean.loc[day:day + delay - pd.Timedelta(days=1)]
        led_increase = led_after.mean() > led_before.mean()
        led_summary = LED_aic_summary(np.arange(len(led_window)), led_window.to_numpy())
        led_detected = bool(
            led_increase
            and led_summary["model"] in (1, 2)
            and led_summary["delta_aic"] > led_min_aic_delta
        )
        rows.extend(
            (
                {"date": day + delay, "test": "led", "detected": led_detected},
                {"date": day + delay, "test": "ad", "detected": ad_detected},
            )
        )

    return pd.DataFrame(rows)


def biomarkers_detection_rates(
    data,
    size=250,
    n_iterations=1000,
    start=pd.Timestamp("2025-07-01"),
    end=pd.Timestamp("2025-10-01"),
    window_size=14,
    alpha=0.05,
    led_min_aic_delta=0.0,
    seed=42,
    log_values=False,
    log_offset=0.0,
    change_date=pd.Timestamp("2025-09-01"),
    tp_window_days=15,
    include_randomized_tn=True,
):
    data = prepare_value_scale(data, log_values=log_values, log_offset=log_offset)
    rng = np.random.default_rng(seed)
    all_rates = []
    tp_start = pd.to_datetime(change_date)
    tp_end = tp_start + pd.Timedelta(tp_window_days - 1, unit="D")

    for parameter_id, parameter_data in data.groupby("parameter_id", sort=False):
        counts = None
        tn_counts = None
        for _ in range(n_iterations):
            sample = get_bootstrap_sample(parameter_data[["date", "value"]], size=size, rng=rng)
            detections = detect_led_and_ad(
                sample,
                start=start,
                end=end,
                window_size=window_size,
                alpha=alpha,
                led_min_aic_delta=led_min_aic_delta,
            )
            current = detections.pivot(index="date", columns="test", values="detected").astype(int)
            counts = current if counts is None else counts.add(current, fill_value=0)

            if include_randomized_tn:
                randomized_tn_sample = get_randomized_period_sample(
                    sample,
                    source_start=tp_start,
                    source_end=tp_end,
                    target_start=tp_start - pd.Timedelta(3 * window_size, unit="D"),
                    target_end=tp_end,
                    size=size,
                    rng=rng,
                )
                tn_detections = detect_led_and_ad(
                    randomized_tn_sample,
                    start=tp_start,
                    end=tp_end,
                    window_size=window_size,
                    alpha=alpha,
                    led_min_aic_delta=led_min_aic_delta,
                )
                current_tn = tn_detections.pivot(
                    index="date", columns="test", values="detected"
                ).astype(int)
                tn_counts = (
                    current_tn if tn_counts is None else tn_counts.add(current_tn, fill_value=0)
                )

        rates = counts.div(n_iterations).reset_index().melt(
            id_vars="date", var_name="test", value_name="detection_rate"
        )
        rates["parameter_id"] = int(parameter_id)
        rates["parameter"] = BIOMARKER_LABELS[int(parameter_id)]
        rates["n_iterations"] = n_iterations
        rates["n_dayperson"] = size
        rates["scenario"] = "observed"
        all_rates.append(rates)

        if include_randomized_tn:
            tn_rates = tn_counts.div(n_iterations).reset_index().melt(
                id_vars="date", var_name="test", value_name="detection_rate"
            )
            tn_rates["parameter_id"] = int(parameter_id)
            tn_rates["parameter"] = BIOMARKER_LABELS[int(parameter_id)]
            tn_rates["n_iterations"] = n_iterations
            tn_rates["n_dayperson"] = size
            tn_rates["scenario"] = "true_negative"
            all_rates.append(tn_rates)

    return pd.concat(all_rates, ignore_index=True)


def summarize_biomarker_detection_rates(
    rates,
    change_date=pd.Timestamp("2025-09-01"),
    tp_window_days=15,
    tn_start=None,
    tn_end=None,
    smoothing_span=7,
):
    rates = rates.copy()
    rates["date"] = pd.to_datetime(rates["date"])
    if "scenario" not in rates.columns:
        raise ValueError(
            "Biomarker rates have no randomized true-negative scenario. "
            "Recompute them with biomarkers_detection_rates(..., include_randomized_tn=True)."
        )
    observed_rates = rates.loc[rates["scenario"].eq("observed")].copy()
    negative_rates = rates.loc[rates["scenario"].eq("true_negative")].copy()
    if observed_rates.empty or negative_rates.empty:
        raise ValueError("Both 'observed' and 'true_negative' scenarios are required")
    change_date = pd.to_datetime(change_date)
    if tn_start is None:
        tn_start = observed_rates["date"].min()
    if tn_end is None:
        tn_end = change_date - pd.Timedelta(days=1)
    return summarize_detection_rate_metrics(
        observed_rates,
        tp_start=change_date,
        tp_end=change_date + pd.Timedelta(tp_window_days - 1, unit="D"),
        tn_start=tn_start,
        tn_end=tn_end,
        reference_date=change_date,
        smoothing_span=smoothing_span,
        group_columns=("parameter_id", "parameter", "test"),
        negative_detection_rates=negative_rates,
    )


def plot_biomarkers_detection_rates(
    rates,
    save_path=None,
    smoothing_span=7,
    change_date=pd.Timestamp("2025-09-01"),
    tp_window_days=15,
    start_date=None,
    font_size=18,
):
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    configure_plot_font("Inter")
    if isinstance(rates, (str, Path)):
        rates = pd.read_csv(rates, parse_dates=["date"])
    else:
        rates = rates.copy()
        rates["date"] = pd.to_datetime(rates["date"])

    if "scenario" in rates.columns:
        rates = rates.loc[rates["scenario"].eq("observed")].copy()

    rates = rates.loc[rates["parameter_id"].astype(int).ne(1012)].copy()
    if start_date is not None:
        start_date = pd.to_datetime(start_date).normalize()
        if rates.empty or start_date > rates["date"].max():
            raise ValueError("No detection rates remain in the requested plot date range")

    change_date = pd.to_datetime(change_date).normalize()
    min_date = start_date if start_date is not None else rates["date"].min().normalize()
    max_date = rates["date"].max().normalize()
    weekly_ticks = anchored_weekly_ticks(min_date, max_date, anchor=change_date)
    tp_end = change_date + pd.Timedelta(tp_window_days, unit="D")

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharex=True, sharey=True)
    for ax, test in zip(axes, ("led", "ad")):
        for parameter_id, parameter_rates in rates.loc[rates["test"] == test].groupby("parameter_id"):
            parameter_id = int(parameter_id)
            parameter_rates = parameter_rates.sort_values("date")
            smoothed = parameter_rates["detection_rate"].rolling(
                smoothing_span, center=True, min_periods=smoothing_span
            ).mean()
            color = BIOMARKER_COLORS[parameter_id]
            marker = BIOMARKER_MARKERS[parameter_id]
            is_crp = parameter_id == 8
            ax.scatter(
                parameter_rates["date"], parameter_rates["detection_rate"],
                color=color,
                marker=marker,
                s=44 if is_crp else 26,
                alpha=0.55 if is_crp else 0.35,
                zorder=4 if is_crp else 2,
            )
            ax.plot(
                parameter_rates["date"],
                smoothed,
                color=color,
                linewidth=3 if is_crp else 2,
                zorder=4 if is_crp else 2,
            )
            ax.plot([], [], color=color, marker=marker,
                    linewidth=3 if is_crp else 2, markersize=9 if is_crp else 7,
                    label=BIOMARKER_LABELS[parameter_id])

        ax.axvspan(
            change_date,
            tp_end,
            color=PLOT_COLORS["true_positive_fill"],
            alpha=0.28,
            zorder=0,
        )
        ax.axvline(change_date, color="#000000", linestyle="--", linewidth=1.3)
        ax.set_xlabel("Date", fontsize=font_size)
        ax.set_ylim(-0.02, 1.02)
        if start_date is not None:
            ax.set_xlim(left=start_date)
        ax.tick_params(axis="both", labelsize=font_size - 2)
        ax.set_xticks(weekly_ticks)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.legend(frameon=False, fontsize=font_size - 2, ncol=1)

    axes[0].set_ylabel("Detection rate", fontsize=font_size)
    fig.autofmt_xdate(rotation=30, ha="right")
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, axes


def analyze_real_biomarkers(
    data_path=Path("../data/biomarkers_real_data.tsv"),
    crp_data_path=Path("../data/crp_real_data.tsv"),
    include_crp=True,
    rates_output_path=Path("../metrics/biomarkers_detection_rates.csv"),
    target_metrics_output_path=Path("../metrics/biomarkers_detection_rate_metrics.csv"),
    plot_output_path=Path("../figures/biomarkers_detection_rates.png"),
    plot_start_date=None,
    **rate_kwargs,
):
    data = load_biomarkers_data(
        data_path,
        crp_path=crp_data_path,
        include_crp=include_crp,
    )
    rates = biomarkers_detection_rates(data, **rate_kwargs)
    target_metrics = summarize_biomarker_detection_rates(rates)
    rates.attrs["target_metrics"] = target_metrics
    if rates_output_path is not None:
        Path(rates_output_path).parent.mkdir(parents=True, exist_ok=True)
        rates.to_csv(rates_output_path, index=False)
    if target_metrics_output_path is not None:
        Path(target_metrics_output_path).parent.mkdir(parents=True, exist_ok=True)
        target_metrics.to_csv(target_metrics_output_path, index=False)
    if plot_output_path is not None:
        Path(plot_output_path).parent.mkdir(parents=True, exist_ok=True)
        plot_biomarkers_detection_rates(
            rates,
            save_path=plot_output_path,
            start_date=plot_start_date,
        )
    return rates


if __name__ == "__main__":
    for n_dayperson in [250, 50, 10]:

        rates = analyze_real_biomarkers(
            size=n_dayperson,
            n_iterations=1000,
            window_size=14,
            alpha=0.05,
            seed=42,
            rates_output_path=f"../metrics/raw/biomarkers_detection_rates_{n_dayperson}.csv",
            target_metrics_output_path=f"../metrics/raw/biomarkers_detection_rate_metrics_{n_dayperson}.csv",
            plot_output_path=f"../figures/biomarkers_detection_rates_{n_dayperson}.png",
        )

        rates.to_csv(f"../metrics/raw/biomarkers_detection_rates_{n_dayperson}.csv", index=False)
        rates = pd.read_csv(f"../metrics/raw/biomarkers_detection_rates_{n_dayperson}.csv")

        plot_biomarkers_detection_rates(
            rates,
            start_date=pd.to_datetime("2025-08-15"),
            save_path=f"../figures/biomarkers_detection_rates_{n_dayperson}.png",
            smoothing_span=7,
        )
