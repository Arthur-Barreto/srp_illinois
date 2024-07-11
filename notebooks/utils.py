# pylint: disable=missing-docstring
import numpy as np
import numpy.typing as npt


def autophase(period, duration):
    return np.arange(0, period, period / duration)


def auto_max_min_period(t: npt.NDArray):
    # Compute the baseline.
    total_duration = np.ptp(t)

    # Compute the minimum and maximum periods.
    minimum_period = 2.0 * np.min(np.diff(np.sort(t)))
    maximum_period = 0.5 * total_duration

    return minimum_period, maximum_period, total_duration


def autoperiod(
    minimum_period: float,
    maximum_period: float,
    total_duration: float,
):
    # Check that the minimum and maximum periods are positive.
    assert minimum_period > 0, "minimum_period must be positive"
    assert maximum_period > 0, "maximum_period must be positive"

    # Check that the minimum period is less than the maximum period.
    assert (minimum_period < maximum_period), \
        "minimum_period must be less than maximum_period"

    # Check that the maximum period does not exceed the total duration.
    assert maximum_period <= total_duration, \
        "maximum_period must not exceed total_duration"

    # Convert bounds to frequency
    minimum_frequency = 1.0 / maximum_period
    maximum_frequency = 1.0 / minimum_period

    # Estimate the required frequency spacing
    frequency_resolution = 1.0 / total_duration

    # Because of the sparsity of a transit, this must be much finer than
    # the frequency resolution for a sinusoidal fit. For a sinusoidal fit,
    # df would be 1/baseline (see LombScargle), but here this should be
    # scaled proportionally to the duration in units of baseline.
    # frequency_resolution *= frequency_resolution / maximum_frequency

    # Compute the number of frequencies and the frequency grid
    frequency_grid = np.arange(
        maximum_frequency,
        minimum_frequency,
        -frequency_resolution,
    )
    periods = 1.0 / frequency_grid
    return periods


def autoperiod_astropy(
    t: npt.NDArray,
    y: npt.NDArray,
    dy: npt.NDArray,
    duration: float | list[float],
    minimum_period: float | None = None,
    maximum_period: float | None = None,
    minimum_n_transit: int = 3,
    frequency_factor: float = 1.0,
):
    # Check that the input arrays have the same length.
    all_equal_length = len(t) == len(y) == len(dy)
    assert all_equal_length, "t, y, and dy must have the same length"

    # Check that the duration is a scalar.
    assert isinstance(duration, (float, list)), \
        "duration must be a scalar or a list"

    if isinstance(duration, float):
        duration = [duration]

    # Check that all durations are float.
    assert all(isinstance(d, float) for d in duration), \
        "all durations must be float"

    # Check that the duration is positive.
    assert all(d > 0 for d in duration), "duration must be positive"

    # Check that the minimum period is positive.
    if minimum_period is not None:
        assert minimum_period > 0, "minimum_period must be positive"
    else:
        minimum_period = 2.0 * np.max(duration)

    # Compute the baseline.
    baseline = np.ptp(t)
    min_duration = np.min(duration)

    # Check that the maximum period is positive.
    if maximum_period is not None:
        assert maximum_period > 0, "maximum_period must be positive"
    else:
        maximum_period = baseline

    # Check that the minimum period is less than the maximum period.
    assert (
        minimum_period
        < maximum_period), "minimum_period must be less than maximum_period"

    # Check that the minimum number of transits is positive.
    assert minimum_n_transit > 0, "minimum_n_transit must be positive"

    # Check that the frequency factor is positive.
    assert frequency_factor > 0, "frequency_factor must be positive"

    # Estimate the required frequency spacing
    # Because of the sparsity of a transit, this must be much finer than
    # the frequency resolution for a sinusoidal fit. For a sinusoidal fit,
    # df would be 1/baseline (see LombScargle), but here this should be
    # scaled proportionally to the duration in units of baseline.
    df = frequency_factor * min_duration / baseline**2

    # Convert bounds to frequency
    minimum_frequency = 1.0 / maximum_period
    maximum_frequency = 1.0 / minimum_period

    # Compute the number of frequencies and the frequency grid
    nf = 1 + int(np.round((maximum_frequency - minimum_frequency) / df))
    return 1.0 / (maximum_frequency - df * np.arange(nf))
