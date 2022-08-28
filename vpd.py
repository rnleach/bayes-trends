from collections.abc import Iterable
from datetime import date, timedelta

from formulas import calc_dp, vapor_pressure_liquid_water

def vapor_pressure_deficit(t_c, rh):
    """Given the temperature and humidity, calculate the vapor pressure deficit.

    # Arguments
    t_c - the temperature in Celsius
    rh - the relative humidity in percent, 0 to 100

    # Returns
    The vapor pressure deficit in hectopascals.
    """
    tdc = calc_dp(t_c, rh)

    saturation_vapor_pressure = vapor_pressure_liquid_water(t_c)
    vapor_pressure = vapor_pressure_liquid_water(tdc)

    return saturation_vapor_pressure - vapor_pressure


def get_daily_max_min_vpd(archive, site, starting=None, ending=None, break_hour=6):
    """Get the daily maximum and minimum vapor pressure deficit.

    # Arguments
    archive - an archive from which to retrieve the data.
    site - the site id, e.g. 'kmso', 'smtm8'
    starting - when to start the time series
    end - when to end the time series
    break_hour - subtract this many hours from the valid time to assign the day. This means the 
                 early morning observations will go with the previous day.

    # Returns a generator object that yields a tuple of (date, maximum daily vpd, minimum daily vpd).

    """
    data_gen = archive.get_data_sets_for(site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((vt - timedelta(hours=break_hour), t, rh) for _, vt, t, rh in data_gen)
    # Filter out data that isn't for the month we're interested in.
    data_gen = ((vt, t, rh) for vt, t, rh in data_gen if vt.month == month)
    # Filter out  humidity values below 1%
    data_gen = ((vt, t, rh) for vt, t, rh in data_gen if rh > 0.9)
    # map the temperature and rh to vpd
    data_gen = ((vt, vapor_pressure_deficit(t, rh)) for vt, t, rh in data_gen)

    # Initialize the first time step and the accumulator variables
    vt, vpd = next(data_gen)
    curr_year = vt.year
    curr_month = vt.month
    curr_day = vt.day
    max_vpd = vpd
    min_vpd = vpd

    for vt, vpd in data_gen:

        year = vt.year
        month = vt.month
        day = vt.day

        if day != curr_day or month != curr_month or year != curr_year:

            # Yield the values for the current day, it is done.
            yield (date(curr_year, curr_month, curr_day), max_vpd, min_vpd)

            # Update the current date and reset the accumulator variables
            curr_year = year
            curr_month = month
            curr_day = day
            max_vpd = vpd
            min_vpd = vpd

        # Always update the accumulator variables!
        max_vpd = max(max_vpd, vpd)
        min_vpd = min(min_vpd, vpd)

    # Don't forget to yield the last day's values!
    yield (date(curr_year, curr_month, curr_day), max_vpd, min_vpd)

    return


def get_monthly_average_vpd(archive, site, months, starting=None, ending=None, break_hour=6):
    """Get the monthly average vapor pressure deficit.

    Averaging is done on the hourly observations.

    # Arguments
    archive - an archive from which to retrieve the data.
    site - the site id, e.g. 'kmso', 'smtm8'.
    months - 1 to 12 for the number of the month, or a collection of months.
    starting - when to start the time series.
    end - when to end the time series.
    break_hour - subtract this many hours from the valid time to assign the day. This means the 
                 early morning observations will go with the previous day.

    # Returns
    A generator object that yields a tuple of (year, average vpd) where the average is for the
    requested month(s).

    """
    if isinstance(months, Iterable):
        months = tuple(months)
    else:
        months = (months,)

    # Heuristic test to see if too much missing data.
    avg_hours = len(months) * 30 * 24
    min_allowed_hours = int(round(0.9 * avg_hours))

    data_gen = archive.get_hourly_temp_rh_for(site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((vt - timedelta(hours=break_hour), t, rh) for _, vt, t, rh in data_gen)
    # Filter out data that isn't for the month(s) we're interested in.
    data_gen = ((vt, t, rh) for vt, t, rh in data_gen if vt.month in months)
    # Filter out invalid humidity values less than 1%
    data_gen = ((vt, t, rh) for vt, t, rh in data_gen if rh > 0.9)
    # map the temperature and rh to vpd
    data_gen = ((vt, vapor_pressure_deficit(t, rh)) for vt, t, rh in data_gen)

    # Initialize the first time step and the accumulator variables
    vt, vpd = next(data_gen)
    curr_year = vt.year
    sum_vpd = vpd
    count = 1

    for vt, vpd in data_gen:

        if curr_year != vt.year:

            # Yield the values for the current year
            if count >= min_allowed_hours:
                yield (curr_year, sum_vpd / count)

            # Update the current year and reset the accumulator variables
            curr_year = vt.year
            sum_vpd = 0.0
            count = 0

        # Always update the accumulator variables!
        sum_vpd += vpd
        count += 1

    # Don't forget to yield the last set of values!
    if count >= min_allowed_hours:
        yield (curr_year, sum_vpd / count)

    return


