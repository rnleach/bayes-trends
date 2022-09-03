from collections.abc import Iterable
from datetime import date, timedelta

import matplotlib.pyplot as plt

from formulas import vapor_pressure_liquid_water, specific_humidity

def vapor_pressure_deficit(t_c, td_c):
    """Given the temperature and dew point, calculate the vapor pressure deficit.

    # Arguments
    t_c - the temperature in Celsius
    td_c - the dew point in Celsius.

    # Returns
    The vapor pressure deficit in hectopascals.
    """

    saturation_vapor_pressure = vapor_pressure_liquid_water(t_c)
    vapor_pressure = vapor_pressure_liquid_water(td_c)

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
    data_gen = archive.get_data_for(site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((vt - timedelta(hours=break_hour), t, dew, pres) for _, vt, t, dew, pres in data_gen)
    # Filter out data that isn't for the month we're interested in.
    data_gen = ((vt, t, dew, pres) for vt, t, dew, pres in data_gen if vt.month == month)
    # map the temperature and dew point to vpd
    data_gen = ((vt, vapor_pressure_deficit(t, dew)) for vt, t, dew, pres in data_gen)

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

    data_gen = archive.get_hourly_data_for(site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((vt - timedelta(hours=break_hour), t, dew) for _, vt, t, dew, _ in data_gen)
    # Filter out data that isn't for the month(s) we're interested in.
    data_gen = ((vt, t, dew) for vt, t, dew, _ in data_gen if vt.month in months)
    # map the temperature and dew point to vpd
    data_gen = ((vt, vapor_pressure_deficit(t, dew)) for vt, t, dew in data_gen)

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


def get_averaged_data(archive, site, months, starting=None, ending=None, break_hour=6):
    """Get the averaged temperature, dew point, specific humidity, and vapor pressure deficit.

    Averaging is done on the hourly observations. Average is over the whole month or list of months.

    # Arguments
    archive - an archive from which to retrieve the data.
    site - the site id, e.g. 'kmso', 'smtm8'.
    months - 1 to 12 for the number of the month, or a collection of months.
    starting - when to start the time series.
    end - when to end the time series.
    break_hour - subtract this many hours from the valid time to assign the day. This means the 
                 early morning observations will go with the previous day.

    # Returns
    A generator object that yields a tuple of (year, average temperature, average dew point, 
    average pressure, average specific humidity, average vpd).

    """
    if isinstance(months, Iterable):
        months = tuple(months)
    else:
        months = (months,)

    # Heuristic test to see if too much missing data.
    avg_hours = len(months) * 30 * 24
    min_allowed_hours = int(round(0.9 * avg_hours))

    data_gen = archive.get_hourly_data_for(site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((vt - timedelta(hours=break_hour), t, dew, pres) for _, vt, t, dew, pres in data_gen)
    # Filter out data that isn't for the month(s) we're interested in.
    data_gen = ((vt, t, dew, pres) for vt, t, dew, pres in data_gen if vt.month in months)
    # map the temperature and dew point to vpd
    data_gen = ((vt, t, dew, pres, specific_humidity(dew, pres),  vapor_pressure_deficit(t, dew)) for vt, t, dew, pres in data_gen)

    # Initialize the first time step and the accumulator variables
    vt, t, td, pres, sh, vpd = next(data_gen)
    curr_year = vt.year
    sum_t = t
    sum_td = td
    sum_vpd = vpd
    count = 1

    if pres is not None and sh is not None:
        sum_pres = pres
        sum_sh = sh
        count_pres = 1
    else:
        sum_pres = 0.0
        sum_sh = 0.0
        count_pres = 0

    for vt, t, td, pres, sh, vpd in data_gen:

        if curr_year != vt.year:

            # Yield the values for the current year
            if count >= min_allowed_hours:
                if count_pres >= min_allowed_hours:
                    avg_pres = sum_pres / count_pres
                    avg_sh = sum_sh / count_pres
                else:
                    avg_pres = float('nan')
                    avg_sh = float('nan')

                avg_t = sum_t / count
                avg_td = sum_td / count
                avg_vpd = sum_vpd / count
                yield (curr_year, avg_t, avg_td, avg_pres, avg_sh, avg_vpd)

            # Update the current year and reset the accumulator variables
            curr_year = vt.year
            sum_t = 0.0
            sum_td = 0.0
            sum_vpd = 0.0
            sum_pres = 0.0
            sum_sh = 0.0

            count = 0
            count_pres = 0

        # Always update the accumulator variables!
        sum_t += t
        sum_td += td
        sum_vpd += vpd
        count += 1

        if pres is not None and sh is not None:
            sum_pres += pres
            sum_sh += sh
            count_pres += 1

    # Don't forget to yield the last set of values!
    if count >= min_allowed_hours:
        if count_pres >= min_allowed_hours:
            avg_pres = sum_pres / count_pres
            avg_sh = sum_sh / count_pres
        else:
            avg_pres = float('nan')
            avg_sh = float('nan')

        avg_t = sum_t / count
        avg_td = sum_td / count
        avg_vpd = sum_vpd / count
        yield (curr_year, avg_t, avg_td, avg_pres, avg_sh, avg_vpd)

    return


def normalize_var(df, input_col, output_col):
    """Normalize a Pandas DataFrame column and add it back into the DataFrame as a new column.

    Normalize means subtract the mean and then divide by the standard deviation.

    # Arguments
    df - the Pandas DataFrame
    input_col - The name of the column in df to be normalized.
    output_col - The name of the column to add back to the DataFrame df with the normalized data.

    # Returns a tuple of two functions. The first is a function that can normalize new data values
    and the second is a function that can denormalize data values.
    """
    offset = df[input_col].mean()
    scale = df[input_col].std()
    
    print(f"{input_col} offset = {offset} scale = {scale}")
    
    df[output_col] = (df[input_col] - offset) / scale
    
    def denorm(x):
        return x * scale + offset
    def norm(x):
        return (x - offset) / scale
    
    return (norm, denorm)

def plot_corrs(df):

    fig_corr, axes = plt.subplots(5,5, sharex='col', figsize=(8,8))
    ((ax_ht, *_),(ax_td_t, ax_htd, *_),(ax_pres_t, ax_pres_td, ax_hpres, *_),(ax_sh_t, ax_sh_td, ax_sh_pres, ax_hsh, _),(ax_vpd_t, ax_vpd_td, ax_vpd_pres, ax_vpd_sh, ax_hvpd)) = axes

    fig_corr.delaxes(axes[0][1])
    fig_corr.delaxes(axes[0][2])
    fig_corr.delaxes(axes[0][3])
    fig_corr.delaxes(axes[0][4])
    fig_corr.delaxes(axes[1][2])
    fig_corr.delaxes(axes[1][3])
    fig_corr.delaxes(axes[1][4])
    fig_corr.delaxes(axes[2][3])
    fig_corr.delaxes(axes[2][4])
    fig_corr.delaxes(axes[3][4])

    # Row 1
    ax_ht.hist(df['y_t_obs']);
    ax_ht.set_yticklabels([]);
    ax_ht.set_ylabel("Temp");

    # Row 2
    ax_htd.hist(df['y_td_obs']);
    ax_htd.set_yticklabels([]);

    ax_td_t.scatter(df['y_t_obs'], df['y_td_obs']);
    ax_td_t.set_yticklabels([]);
    ax_td_t.set_ylabel("DP");

    # Row 3
    ax_hpres.hist(df['y_pres_obs']);
    ax_hpres.set_yticklabels([]);

    ax_pres_t.scatter(df['y_t_obs'], df['y_pres_obs']);
    ax_pres_t.set_yticklabels([]);
    ax_pres_t.set_ylabel("Pres");

    ax_pres_td.scatter(df['y_td_obs'], df['y_pres_obs']);
    ax_pres_td.set_yticklabels([]);

    # Row 4
    ax_hsh.hist(df['y_sh_obs']);
    ax_hsh.set_yticklabels([]);

    ax_sh_t.scatter(df['y_t_obs'], df['y_sh_obs']);
    ax_sh_t.set_yticklabels([]);
    ax_sh_t.set_ylabel("S");

    ax_sh_td.scatter(df['y_td_obs'], df['y_sh_obs']);
    ax_sh_td.set_yticklabels([]);

    ax_sh_pres.scatter(df['y_pres_obs'], df['y_sh_obs']);
    ax_sh_pres.set_yticklabels([]);

    # Row 5
    ax_hvpd.hist(df['y_vpd_obs']);
    ax_hvpd.set_yticklabels([]);
    ax_hvpd.set_xlabel("VPD");
    ax_hvpd.set_xticklabels([]);

    ax_vpd_t.scatter(df['y_t_obs'], df['y_vpd_obs']);
    ax_vpd_t.set_yticklabels([]);
    ax_vpd_t.set_ylabel("VPD");
    ax_vpd_t.set_xlabel("Temp");
    ax_vpd_t.set_xticklabels([]);

    ax_vpd_td.scatter(df['y_td_obs'], df['y_vpd_obs']);
    ax_vpd_td.set_yticklabels([]);
    ax_vpd_td.set_xlabel("Dew");
    ax_vpd_td.set_xticklabels([]);

    ax_vpd_pres.scatter(df['y_pres_obs'], df['y_vpd_obs']);
    ax_vpd_pres.set_yticklabels([]);
    ax_vpd_pres.set_xlabel("Pres");
    ax_vpd_pres.set_xticklabels([]);

    ax_vpd_sh.scatter(df['y_sh_obs'], df['y_vpd_obs']);
    ax_vpd_sh.set_yticklabels([]);
    ax_vpd_sh.set_xlabel("S");
    ax_vpd_sh.set_xticklabels([]);

    return
