import numpy as np
from scipy.stats import mannwhitneyu, anderson_ksamp, ks_2samp
import pandas as pd
import chow_test
from scipy.optimize import curve_fit
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message=".*p-value capped.*")
warnings.filterwarnings("ignore", message=".*p-value floored.*")

TESTS = ["chow", "chow_daily_mean", "led", "led_daily_mean", "ks", "mwu", "ad"]
TEST_LABELS = {
    "chow": "Chow-test",
    "chow_daily_mean": "Chow-test daily",
    "led": "AIC",
    "led_daily_mean": "AIC daily",
    "ks": "Kolmogorov-Smirnov",
    "mwu": "Mann-Whitney",
    "ad": "Anderson-Darling",
}
METRIC_COLORS = {
    "Chow-test": "#0072B2",
    "Chow-test daily": "#44AA99",
    "AIC": "#E69F00",
    "AIC daily": "#882255",
    "Kolmogorov-Smirnov": "#009E73",
    "Mann-Whitney": "#CC79A7",
    "KL-divergence": "#56B4E9",
    "Anderson-Darling": "#332288",
}
TEST_MARKERS = {
    "chow": "o",
    "chow_daily_mean": "v",
    "led": "s",
    "led_daily_mean": "P",
    "ks": "^",
    "mwu": "D",
    "ad": "X",
}
DAILY_TESTS = {"chow_daily_mean", "led_daily_mean"}


def get_tests(include_daily_detectors=True):
    if include_daily_detectors:
        return TESTS
    return [test for test in TESTS if test not in DAILY_TESTS]


PLOT_COLORS = {
    "true_negative_fill": "#D8D8D8",
    "true_positive_fill": "#A6CEE3",
    "outbreak": "#222222",
    "daily_mean": "#B8B8B8",
    "value_smoothing": "#111111",
    "secondary_axis": "#D55E00",
}


def anchored_weekly_ticks(start, end, anchor=pd.Timestamp("2025-09-01")):
    """Return seven-day ticks within a range, anchored to a meaningful date."""
    start = pd.to_datetime(start).normalize()
    end = pd.to_datetime(end).normalize()
    anchor = pd.to_datetime(anchor).normalize()
    if start > end:
        raise ValueError("start must not be later than end")

    first_tick = anchor + pd.Timedelta(
        days=7 * int(np.ceil((start - anchor).days / 7))
    )
    last_tick = anchor + pd.Timedelta(
        days=7 * int(np.floor((end - anchor).days / 7))
    )
    return pd.date_range(first_tick, last_tick, freq="7D")


def load_real_data(
    path=Path("../data/crp_real_data.tsv"),
    parameter_id=None,
    start=pd.Timestamp("2025-07-01"),
    end=pd.Timestamp("2025-10-01"),
):
    data = pd.read_csv(path, sep="\t")
    required = {"date", "value"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if "parameter_id" in data.columns:
        data["parameter_id"] = pd.to_numeric(data["parameter_id"], errors="raise").astype(int)
        available_ids = sorted(data["parameter_id"].unique().tolist())
        if parameter_id is None:
            if len(available_ids) != 1:
                raise ValueError(
                    "The file contains multiple parameters. Pass parameter_id explicitly. "
                    f"Available IDs: {available_ids}"
                )
            parameter_id = available_ids[0]
        parameter_id = int(parameter_id)
        if parameter_id not in available_ids:
            raise ValueError(
                f"parameter_id={parameter_id} is absent. Available IDs: {available_ids}"
            )
        data = data.loc[data["parameter_id"].eq(parameter_id)]
    elif parameter_id is not None:
        raise ValueError("parameter_id was provided, but the file has no 'parameter_id' column")

    data = data.copy()
    data["date"] = pd.to_datetime(data["date"], errors="raise")
    data["value"] = pd.to_numeric(data["value"], errors="raise")
    if start is not None:
        data = data.loc[data["date"] >= pd.to_datetime(start)]
    if end is not None:
        data = data.loc[data["date"] <= pd.to_datetime(end)]
    if data.empty:
        raise ValueError("No observations remain after parameter and date filtering")

    return data[["date", "value"]].sort_values("date").reset_index(drop=True)


def configure_plot_font(font_family="Inter"):
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    font_dirs = [
        Path(__file__).resolve().parent,
        Path.cwd(),
        Path.home() / "AppData/Local/Microsoft/Windows/Fonts",
        Path("C:/Windows/Fonts"),
    ]
    for font_dir in font_dirs:
        if not font_dir.exists():
            continue
        for pattern in (f"{font_family}*.ttf", f"{font_family}*.otf"):
            for font_path in font_dir.glob(pattern):
                font_manager.fontManager.addfont(str(font_path))

    plt.rcParams["font.family"] = font_family


def prepare_value_scale(data, log_values=False, log_offset=0.0):
    data = data.copy()

    if log_values:
        if (data["value"] + log_offset <= 0).any():
            raise ValueError("Log transform requires value + log_offset to be strictly positive")
        data["value"] = np.log(data["value"] + log_offset)

    return data


def value_axis_label(log_values=False, log_offset=0.0):
    if log_values:
        return "CRP, log scale"
    return "CRP concentration"

def calculate_aic(n, rss, k):
    rss = max(float(rss), np.finfo(float).eps)
    return n * np.log(rss / n) + 2 * k

def random_noise_model(x, a):
    return a

def linear_growth_model(x, a, b):
    return a * x + b

def exponential_growth_model(x, a, b, c):
    return a * np.exp(b * x) + c


def LED_aic_summary(x, y):
    window_size = len(x)
    y = np.asarray(y, dtype=float)

    if window_size == 0 or np.isnan(y).any():
        return {
            "model": -1,
            "delta_aic": np.nan,
            "aic_random": np.nan,
            "aic_linear": np.nan,
            "aic_exponential": np.nan,
        }

    try:
        popt_random, _ = curve_fit(random_noise_model, x, y, p0=[np.mean(y)], maxfev=2000)
        popt_linear, _ = curve_fit(linear_growth_model, x, y, p0=[0, np.mean(y)], maxfev=2000)
        popt_exponential, _ = curve_fit(exponential_growth_model, x, y, p0=[1, 0.1, np.mean(y)], maxfev=2000)
    except (RuntimeError, ValueError):
        return {
            "model": -1,
            "delta_aic": np.nan,
            "aic_random": np.nan,
            "aic_linear": np.nan,
            "aic_exponential": np.nan,
        }

    rss_random = np.sum((y - random_noise_model(x, *popt_random))**2)
    rss_linear = np.sum((y - linear_growth_model(x, *popt_linear))**2)
    rss_exponential = np.sum((y - exponential_growth_model(x, *popt_exponential))**2)

    aic_random = calculate_aic(window_size, rss_random, 1)
    aic_linear = calculate_aic(window_size, rss_linear, 2)
    aic_exponential = calculate_aic(window_size, rss_exponential, 3)

    growth_aic_values = [aic_linear, aic_exponential]
    best_growth_model = int(np.argmin(growth_aic_values)) + 1
    best_growth_aic = growth_aic_values[best_growth_model - 1]
    delta_aic = aic_random - best_growth_aic

    return {
        "model": best_growth_model,
        "delta_aic": delta_aic,
        "aic_random": aic_random,
        "aic_linear": aic_linear,
        "aic_exponential": aic_exponential,
    }


def get_bootstrap_sample(data, size=250, rng=None):

    if rng is None:
        rng = np.random.default_rng()

    bootstrap_sample = {
        "date": [],
        "value": [],
    }

    for date in data.date.unique():
        bootstrap_sample["date"] += [date]*size
        bootstrap_sample["value"] += rng.choice(a=data[data.date == date].value, size=size, replace=True).tolist()

    return pd.DataFrame(bootstrap_sample)


def get_randomized_period_sample(
    data,
    source_start,
    source_end,
    target_start,
    target_end,
    size=None,
    rng=None,
):
    if rng is None:
        rng = np.random.default_rng()

    source_start = pd.to_datetime(source_start)
    source_end = pd.to_datetime(source_end)
    target_dates = pd.date_range(target_start, target_end, freq="D")
    source_values = data.loc[
        pd.to_datetime(data.date).between(source_start, source_end),
        "value",
    ].to_numpy()

    if len(source_values) == 0:
        raise ValueError("No values found in the source period for randomized sample")

    rows = []
    daily_sizes = data.groupby("date").size()

    for date in target_dates:
        if size is None:
            n_values = int(daily_sizes.get(date, daily_sizes.median()))
        else:
            n_values = int(size)

        sampled_values = rng.choice(source_values, size=n_values, replace=True)
        rows.append(pd.DataFrame({"date": date, "value": sampled_values}))

    return pd.concat(rows, ignore_index=True)


def do_tests(
    data,
    start=pd.to_datetime("2025-08-01"),
    end=pd.to_datetime("2025-10-01"),
    window_size=14,
    alpha=0.05,
    chow_alpha=0.05,
    led_min_aic_delta=0.0,
    chow_input="mean",
    led_input="mean",
    include_daily_detectors=True,
    include_led_daily_mean=None,
):
    if include_led_daily_mean is not None:
        include_daily_detectors = include_led_daily_mean

    input_aliases = {
        "mean": "mean",
        "rolling_mean": "mean",
        "daily_mean": "mean",
        "variance": "variance",
        "rolling_variance": "variance",
        "daily_variance": "variance",
    }
    if chow_input not in input_aliases:
        raise ValueError(f"Unknown chow_input={chow_input!r}. Use 'mean' or 'variance'")
    if led_input not in input_aliases:
        raise ValueError(f"Unknown led_input={led_input!r}. Use 'mean' or 'variance'")
    chow_mode = input_aliases[chow_input]
    led_mode = input_aliases[led_input]
    rolling_columns = {
        "mean": "window_mean",
        "variance": "window_variance",
    }
    daily_columns = {
        "mean": "daily_mean",
        "variance": "daily_variance",
    }

    detection_delay = pd.Timedelta(window_size, unit="D")
    candidate_start = start - detection_delay
    candidate_end = end - detection_delay
    dates = pd.date_range(
        candidate_start - detection_delay,
        candidate_end + detection_delay,
        freq="D",
    )

    detection_results = pd.DataFrame({
        "date" : dates,
    })

    for day in dates:
        window_data = data.loc[(data.date > day - pd.Timedelta(window_size, unit="D")) & (data.date <= day)]
        detection_results.loc[detection_results.date == day, "window_mean"] = window_data.value.mean()
        detection_results.loc[detection_results.date == day, "window_variance"] = window_data.value.var(ddof=1)
        day_data = data.loc[data.date == day, "value"]
        detection_results.loc[detection_results.date == day, "daily_mean"] = day_data.mean()
        detection_results.loc[detection_results.date == day, "daily_variance"] = day_data.var(ddof=1)


    for day in pd.date_range(candidate_start, candidate_end, freq="D")[::-1]:
        
        before = data.loc[(data.date >= day - pd.Timedelta(window_size, unit="D")) & (data.date < day)]
        after = data.loc[(data.date >= day) & (data.date < day + pd.Timedelta(window_size, unit="D"))]
        raw_windows_available = not before.empty and not after.empty
        if raw_windows_available:
            raw_increase = after.value.mean() > before.value.mean()
            mwu_pvalue = mannwhitneyu(before.value, after.value).pvalue
            ks_pvalue = ks_2samp(before.value, after.value).pvalue
            try:
                ad_pvalue = anderson_ksamp(
                    [before.value, after.value]
                ).significance_level
            except ValueError:
                ad_pvalue = np.nan
        else:
            raw_increase = False
            mwu_pvalue = np.nan
            ad_pvalue = np.nan
            ks_pvalue = np.nan

        detection_results.loc[detection_results.date == day, "mean_delta"] = after.value.mean() - before.value.mean()
        detection_results.loc[detection_results.date == day, "mwu_p_val"] = mwu_pvalue
        detection_results.loc[detection_results.date == day, "ad_p_val"] = ad_pvalue
        detection_results.loc[detection_results.date == day, "ks_p_val"] = ks_pvalue
        detection_results.loc[detection_results.date == day, "mwu"] = (mwu_pvalue < alpha) & raw_increase
        detection_results.loc[detection_results.date == day, "ad"] = (ad_pvalue < alpha) & raw_increase
        detection_results.loc[detection_results.date == day, "ks"] = (ks_pvalue < alpha) & raw_increase
        detection_results.loc[detection_results.date == day, "value"] = after.loc[after.date == day, "value"].mean()


        chow_column = rolling_columns[chow_mode]
        led_column = rolling_columns[led_mode]
        chow_daily_column = daily_columns[chow_mode]
        led_daily_column = daily_columns[led_mode]
        before_means, after_means = detection_results.loc[(detection_results.date >= day - pd.Timedelta(window_size, unit="D")) & (detection_results.date < day)], detection_results.loc[(detection_results.date >= day) & (detection_results.date < day + pd.Timedelta(window_size, unit="D"))]
        chow_increase = after_means[chow_column].mean() > before_means[chow_column].mean()
        led_increase = after_means[led_column].mean() > before_means[led_column].mean()

        chow_before = before_means[chow_column].to_numpy(dtype=float)
        chow_after = after_means[chow_column].to_numpy(dtype=float)
        if np.isnan(chow_before).any() or np.isnan(chow_after).any():
            chow_pvalue = np.nan
        else:
            chow_pvalue = chow_test.p_value(
                y1=chow_before,
                x1=np.arange(len(chow_before)),
                y2=chow_after,
                x2=np.arange(len(chow_before), len(chow_before) + len(chow_after)),
            )
        detection_results.loc[detection_results.date == day, "chow_p_val"] = chow_pvalue

        detection_results.loc[detection_results.date == day, "chow_input"] = chow_mode
        detection_results.loc[detection_results.date == day, "chow_input_delta"] = after_means[chow_column].mean() - before_means[chow_column].mean()
        detection_results.loc[detection_results.date == day, "smoothed_mean_delta"] = after_means.window_mean.mean() - before_means.window_mean.mean()
        detection_results.loc[detection_results.date == day, "chow"] = (detection_results.loc[detection_results.date == day, "chow_p_val"] < chow_alpha) & chow_increase

        daily_increase = after_means[chow_daily_column].mean() > before_means[chow_daily_column].mean()
        if include_daily_detectors:
            chow_daily_before = before_means[chow_daily_column].to_numpy(dtype=float)
            chow_daily_after = after_means[chow_daily_column].to_numpy(dtype=float)
            if np.isnan(chow_daily_before).any() or np.isnan(chow_daily_after).any():
                chow_daily_pvalue = np.nan
            else:
                chow_daily_pvalue = chow_test.p_value(
                    y1=chow_daily_before,
                    x1=np.arange(len(chow_daily_before)),
                    y2=chow_daily_after,
                    x2=np.arange(len(chow_daily_before), len(chow_daily_before) + len(chow_daily_after)),
                )
            detection_results.loc[detection_results.date == day, "chow_daily_mean_input"] = chow_mode
            detection_results.loc[detection_results.date == day, "chow_daily_mean_p_val"] = chow_daily_pvalue
            detection_results.loc[detection_results.date == day, "chow_daily_mean_delta"] = after_means[chow_daily_column].mean() - before_means[chow_daily_column].mean()
            detection_results.loc[detection_results.date == day, "chow_daily_mean"] = (chow_daily_pvalue < chow_alpha) & daily_increase
        else:
            detection_results.loc[detection_results.date == day, "chow_daily_mean"] = np.nan

        window_data = detection_results.loc[(detection_results.date >= day - pd.Timedelta(window_size, unit="D")) & (detection_results.date < day + pd.Timedelta(window_size, unit="D"))]

        led_summary = LED_aic_summary(
            np.arange(len(window_data.date)),
            window_data[led_column],
        )
        led_model = (
            led_summary["model"]
            if led_summary["model"] != -1 and led_summary["delta_aic"] > led_min_aic_delta
            else led_summary["model"] if led_summary["model"] == -1 else 0
        )
        detection_results.loc[detection_results.date == day, "led_input"] = led_mode
        detection_results.loc[detection_results.date == day, "led_input_delta"] = after_means[led_column].mean() - before_means[led_column].mean()
        detection_results.loc[detection_results.date == day, "led_delta_aic"] = led_summary["delta_aic"]
        detection_results.loc[detection_results.date == day, "led_model"] = led_summary["model"]
        detection_results.loc[detection_results.date == day, "led"] = (led_model in [1,2]) & led_increase

        if include_daily_detectors:
            led_daily_summary = LED_aic_summary(
                np.arange(len(window_data.date)),
                window_data[led_daily_column],
            )
            led_daily_model = (
                led_daily_summary["model"]
                if led_daily_summary["model"] != -1 and led_daily_summary["delta_aic"] > led_min_aic_delta
                else led_daily_summary["model"] if led_daily_summary["model"] == -1 else 0
            )
            led_daily_increase = after_means[led_daily_column].mean() > before_means[led_daily_column].mean()
            detection_results.loc[detection_results.date == day, "led_daily_mean_input"] = led_mode
            detection_results.loc[detection_results.date == day, "led_daily_mean_delta_aic"] = led_daily_summary["delta_aic"]
            detection_results.loc[detection_results.date == day, "led_daily_mean_model"] = led_daily_summary["model"]
            detection_results.loc[detection_results.date == day, "led_daily_mean"] = (led_daily_model in [1,2]) & led_daily_increase
        else:
            detection_results.loc[detection_results.date == day, "led_daily_mean"] = np.nan

    detection_results["candidate_date"] = detection_results["date"]
    detection_results["detection_date"] = detection_results["candidate_date"] + detection_delay

    return detection_results



def bootstrap_detection_rates(
    data,
    size=5,
    n_iterations=5,
    start=pd.to_datetime("2025-07-01"),
    end=pd.to_datetime("2025-10-01"),
    window_size=14,
    alpha=0.05,
    chow_alpha=1e-5,
    led_min_aic_delta=0.0,
    seed=42,
    log_values=False,
    log_offset=0.0,
    chow_input="mean",
    led_input="mean",
    include_daily_detectors=True,
    change_date=pd.Timestamp("2025-09-01"),
    tp_window_days=15,
    include_randomized_tn=True,
):
    data = prepare_value_scale(data, log_values=log_values, log_offset=log_offset)
    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    tests = get_tests(include_daily_detectors=include_daily_detectors)
    rng = np.random.default_rng(seed)
    detection_dates = pd.date_range(start, end, freq="D")
    detection_counts = pd.DataFrame(0, index=detection_dates, columns=tests)
    tp_start = pd.to_datetime(change_date)
    tp_end = tp_start + pd.Timedelta(tp_window_days - 1, unit="D")
    tn_detection_dates = pd.date_range(tp_start, tp_end, freq="D")
    tn_detection_counts = pd.DataFrame(0, index=tn_detection_dates, columns=tests)

    for i in range(n_iterations):
        bootstrap_sample = get_bootstrap_sample(data, size=size, rng=rng)
        results = do_tests(
            bootstrap_sample,
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
        results = results.loc[results.detection_date.isin(detection_dates)].set_index("detection_date")

        for test in tests:
            detection_counts[test] += results[test].reindex(detection_dates).eq(True).astype(int)

        if include_randomized_tn:
            randomized_tn_sample = get_randomized_period_sample(
                bootstrap_sample,
                source_start=tp_start,
                source_end=tp_end,
                target_start=tp_start - pd.Timedelta(3 * window_size, unit="D"),
                target_end=tp_end,
                size=size,
                rng=rng,
            )
            tn_results = do_tests(
                randomized_tn_sample,
                start=tp_start,
                end=tp_end,
                window_size=window_size,
                alpha=alpha,
                chow_alpha=chow_alpha,
                led_min_aic_delta=led_min_aic_delta,
                chow_input=chow_input,
                led_input=led_input,
                include_daily_detectors=include_daily_detectors,
            ).set_index("detection_date")
            for test in tests:
                tn_detection_counts[test] += (
                    tn_results[test].reindex(tn_detection_dates).eq(True).astype(int)
                )

    detection_rates = (
        detection_counts
        .div(n_iterations)
        .reset_index(names="date")
        .melt(
            id_vars="date",
            var_name="test",
            value_name="detection_rate",
        )
    )
    detection_rates["n_iterations"] = n_iterations
    detection_rates["scenario"] = "observed"

    if include_randomized_tn:
        tn_detection_rates = (
            tn_detection_counts
            .div(n_iterations)
            .reset_index(names="date")
            .melt(id_vars="date", var_name="test", value_name="detection_rate")
        )
        tn_detection_rates["n_iterations"] = n_iterations
        tn_detection_rates["scenario"] = "true_negative"
        detection_rates = pd.concat([detection_rates, tn_detection_rates], ignore_index=True)

    return detection_rates


def summarize_detection_rate_metrics(
    detection_rates,
    tp_start,
    tp_end,
    tn_start=None,
    tn_end=None,
    tn_dates=None,
    reference_date=None,
    smoothing_span=7,
    detection_threshold=0.5,
    group_columns=("test",),
    negative_detection_rates=None,
):
    rates = detection_rates.copy()
    rates["date"] = pd.to_datetime(rates["date"])
    tp_start = pd.to_datetime(tp_start)
    tp_end = pd.to_datetime(tp_end)
    reference_date = tp_start if reference_date is None else pd.to_datetime(reference_date)
    group_columns = list(group_columns)
    if negative_detection_rates is not None:
        negative_rates = negative_detection_rates.copy()
        negative_rates["date"] = pd.to_datetime(negative_rates["date"])
    else:
        negative_rates = None

    missing = set(group_columns + ["date", "detection_rate"]).difference(rates.columns)
    if missing:
        raise ValueError(f"Missing required detection-rate columns: {sorted(missing)}")
    if smoothing_span < 1:
        raise ValueError("smoothing_span must be at least 1")

    if tn_dates is not None:
        tn_dates = pd.DatetimeIndex(pd.to_datetime(tn_dates))
    elif tn_start is not None or tn_end is not None:
        if tn_start is None or tn_end is None:
            raise ValueError("tn_start and tn_end must be provided together")
        tn_start = pd.to_datetime(tn_start)
        tn_end = pd.to_datetime(tn_end)

    rows = []
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, group in rates.groupby(grouper, sort=False, dropna=False):
        keys = (keys,) if len(group_columns) == 1 else tuple(keys)
        curve = (
            group.groupby("date", as_index=False)["detection_rate"]
            .mean()
            .sort_values("date")
            .reset_index(drop=True)
        )
        curve["smoothed_detection_rate"] = curve["detection_rate"].rolling(
            window=smoothing_span,
            center=True,
            min_periods=smoothing_span,
        ).mean()

        tp_mask = curve["date"].between(tp_start, tp_end)
        if negative_rates is not None:
            negative_group = negative_rates
            for column, key in zip(group_columns, keys):
                negative_group = negative_group.loc[negative_group[column].eq(key)]
            peak_tn = negative_group["detection_rate"].max()
            tn_mask = None
        elif tn_dates is not None:
            tn_mask = curve["date"].isin(tn_dates)
        elif tn_start is not None:
            tn_mask = curve["date"].between(tn_start, tn_end)
        else:
            tn_mask = curve["date"] < tp_start

        peak_tp = curve.loc[tp_mask, "detection_rate"].max()
        if negative_rates is None:
            peak_tn = curve.loc[tn_mask, "detection_rate"].max()
        if pd.isna(peak_tp) or pd.isna(peak_tn):
            peak_ratio = np.nan
        elif peak_tn == 0:
            peak_ratio = np.inf if peak_tp > 0 else np.nan
        else:
            peak_ratio = peak_tp / peak_tn

        crossing_date = pd.NaT
        eligible = curve.index[curve["date"].between(reference_date, tp_end)]
        for index in eligible:
            current_rate = curve.at[index, "smoothed_detection_rate"]
            if pd.isna(current_rate) or current_rate < detection_threshold:
                continue

            crossing_date = curve.at[index, "date"]
            if index > 0:
                previous_rate = curve.at[index - 1, "smoothed_detection_rate"]
                previous_date = curve.at[index - 1, "date"]
                if (
                    pd.notna(previous_rate)
                    and previous_rate < detection_threshold
                    and current_rate > previous_rate
                ):
                    fraction = (
                        (detection_threshold - previous_rate)
                        / (current_rate - previous_rate)
                    )
                    crossing_date = previous_date + fraction * (crossing_date - previous_date)
                    crossing_date = max(crossing_date, reference_date)
            break

        days_to_threshold = (
            (crossing_date - reference_date) / pd.Timedelta(days=1)
            if pd.notna(crossing_date)
            else np.nan
        )
        row = dict(zip(group_columns, keys))
        row.update({
            "peak_detection_rate_tp": peak_tp,
            "peak_tp_to_tn_ratio": peak_ratio,
            "days_to_50pct_detection": float(days_to_threshold),
        })
        rows.append(row)

    return pd.DataFrame(rows)


def summarize_real_detection_rates(
    detection_rates,
    change_date=pd.Timestamp("2025-09-01"),
    tp_window_days=15,
    tn_start=None,
    tn_end=None,
    smoothing_span=7,
):
    change_date = pd.to_datetime(change_date)
    rates = detection_rates.copy()
    rates["date"] = pd.to_datetime(rates["date"])
    if "scenario" not in rates.columns:
        raise ValueError(
            "Real-data rates have no randomized true-negative scenario. "
            "Recompute them with bootstrap_detection_rates(..., include_randomized_tn=True)."
        )
    observed_rates = rates.loc[rates["scenario"].eq("observed")].copy()
    negative_rates = rates.loc[rates["scenario"].eq("true_negative")].copy()
    if observed_rates.empty or negative_rates.empty:
        raise ValueError("Both 'observed' and 'true_negative' scenarios are required")
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
        negative_detection_rates=negative_rates,
    )


def plot_bootstrap_detections(
    data,
    detection_rates,
    change_date=pd.to_datetime("2025-09-01"),
    tp_window_days=15,
    tn_start=pd.to_datetime("2025-07-01"),
    tn_end=pd.to_datetime("2025-07-21"),
    tn_windows=None,
    show_tn_windows=False,
    window_size=14,
    start=None,
    end=None,
    save_path=None,
    log_values=False,
    log_offset=0.0,
    smoothing_span=7,
    font_size=18,
    tests_to_plot=None,
):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    configure_plot_font("Inter")

    data = prepare_value_scale(data, log_values=log_values, log_offset=log_offset)
    data["date"] = pd.to_datetime(data["date"])
    detection_rates = detection_rates.copy()
    detection_rates["date"] = pd.to_datetime(detection_rates["date"])
    if "scenario" in detection_rates.columns:
        detection_rates = detection_rates.loc[detection_rates["scenario"].eq("observed")]
    tests_to_plot = TESTS if tests_to_plot is None else list(tests_to_plot)
    detection_rates = detection_rates.loc[detection_rates["test"].isin(tests_to_plot)]

    if start is None:
        start = min(detection_rates.date.min(), tn_start, change_date)
    else:
        start = pd.to_datetime(start)

    if end is None:
        end = max(
            detection_rates.date.max(),
            tn_end,
            change_date + pd.Timedelta(tp_window_days - 1, unit="D"),
        )
    else:
        end = pd.to_datetime(end)

    plot_dates = pd.date_range(
        start - pd.Timedelta(2 * window_size, unit="D"),
        end,
        freq="D",
    )
    daily_mean = (
        data
        .groupby("date")["value"]
        .mean()
        .reindex(plot_dates)
    )
    rolling_mean = daily_mean.rolling(window=window_size, min_periods=1).mean()

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(12, 12),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1]},
    )
    ax_values, ax_detections = axes

    window_start = change_date
    window_end = change_date + pd.Timedelta(tp_window_days, unit="D")
    if tn_windows is None:
        tn_windows = [(tn_start, tn_end)]

    for ax in axes:
        if show_tn_windows:
            for i, (tn_window_start, tn_window_end) in enumerate(tn_windows):
                ax.axvspan(
                    pd.to_datetime(tn_window_start),
                    pd.to_datetime(tn_window_end) + pd.Timedelta(1, unit="D"),
                    color=PLOT_COLORS["true_negative_fill"],
                    alpha=0.35,
                )
        ax.axvspan(window_start, window_end, color=PLOT_COLORS["true_positive_fill"], alpha=0.28)
        ax.axvline(change_date, color=PLOT_COLORS["outbreak"], linestyle="--", linewidth=1.2)

    if show_tn_windows:
        for tn_window_start, tn_window_end in tn_windows:
            tn_label_x = pd.to_datetime(tn_window_start) + (
                pd.to_datetime(tn_window_end) + pd.Timedelta(1, unit="D") - pd.to_datetime(tn_window_start)
            ) / 2
            ax_values.text(
                tn_label_x,
                0.97,
                "True negative",
                transform=ax_values.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=font_size,
            )
    ax_values.text(
        window_start + (window_end - window_start) / 2,
        0.85,
        "True positive",
        transform=ax_values.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=font_size,
    )

    ax_values.annotate(
        "Back-to-school infections",
        xy=(mdates.date2num(change_date), 0.85),
        xycoords=ax_values.get_xaxis_transform(),
        xytext=(-245, -50),
        textcoords="offset points",
        arrowprops={"arrowstyle": "->", "color": PLOT_COLORS["outbreak"], "linewidth": 1.2},
        fontsize=font_size,
    )

    ax_values.scatter(daily_mean.index, daily_mean.values, color=PLOT_COLORS["value_smoothing"], s=18, alpha=0.10)
    ax_values.plot(rolling_mean.index, rolling_mean.values, color=PLOT_COLORS["value_smoothing"], linewidth=2)
    ax_values.set_ylabel("CRP, log scale", fontsize=font_size)
    ax_values.tick_params(axis="both", labelsize=font_size)

    for test, test_rates in detection_rates.groupby("test"):
        test_label = TEST_LABELS.get(test, test)
        test_rates = test_rates.sort_values("date").copy()
        smoothed_detection_rate = test_rates.detection_rate.rolling(
            window=smoothing_span,
            center=True,
            min_periods=smoothing_span,
        ).mean()
        ax_detections.scatter(
            test_rates.date,
            test_rates.detection_rate,
            color=METRIC_COLORS.get(test_label),
            marker=TEST_MARKERS.get(test, "o"),
            s=24,
            alpha=0.45,
        )
        ax_detections.plot(
            test_rates.date,
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
    ax_detections.set_xlabel("Date", fontsize=font_size)
    ax_detections.set_ylim(-0.02, 1.02)
    ax_detections.legend(loc="best", ncol=1, frameon=False, fontsize=font_size)
    ax_detections.tick_params(axis="both", labelsize=font_size)
    ax_detections.set_xlim(start, end)
    ax_detections.set_xticks(anchored_weekly_ticks(start, end, anchor=change_date))
    ax_detections.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))

    fig.autofmt_xdate(rotation=30, ha="right")
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, axes


if __name__ == "__main__":
    data = load_real_data("../data/crp_real_data.tsv")
    
    for n_dayperson in [250, 50, 10]:
 
        rates = bootstrap_detection_rates(
            data,
            size=n_dayperson,
            n_iterations=1000,
            include_daily_detectors=False,
        )

        rates.to_csv(f"../metrics/raw/bootstrap_detection_rates_{n_dayperson}.csv", index=False)

        rate_metrics = summarize_real_detection_rates(rates)
        rate_metrics.to_csv(
            f"../metrics/raw/bootstrap_detection_rate_metrics_{n_dayperson}.csv",
            index=False,
        )

        rates = pd.read_csv(
            f"../metrics/raw/bootstrap_detection_rates_{n_dayperson}.csv",
            parse_dates=["date"],
        )

        fig, axes = plot_bootstrap_detections(
            data,
            rates,
            start=pd.to_datetime("2025-08-15"),
            save_path=f"../figures/real_bootstrap_detection_plot_{n_dayperson}.png",
            log_values=False,
        )
