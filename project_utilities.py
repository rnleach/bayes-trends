from collections import deque
from collections.abc import Iterable
from datetime import date, timedelta
from pathlib import Path
import pickle

import arviz as az
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from formulas import vapor_pressure_liquid_water

from watermark import watermark
print(f"Import versions for libraries imported by {__name__}")
print(watermark(iversions=True, globals_=globals()))

def get_monthly_average_vpd_all_sites(archive, months, starting=None, ending=None, break_hour=6):
    """Get the monthly average vapor pressure deficit for all sites.

    Averaging is done on the hourly observations.

    # Arguments
    archive - an archive from which to retrieve the data.
    months - 1 to 12 for the number of the month, or a collection of months.
    starting - when to start the time series.
    end - when to end the time series.
    break_hour - subtract this many hours from the valid time to assign the day. This means the 
                 early morning observations will go with the previous day.

    # Returns
    A generator object that yields a tuple of (site id, average latitude, average longitude, 
    average elevation, year, average vpd) where the average is for the
    requested month(s). Latitude, longitude, and elevation should rarely change; only in years when
    the station moved.
    """
    sites = archive.get_all_sites()

    for site in sites:

        data = get_monthly_average_vpd(archive, site, months, starting, ending, break_hour)

        # Check for at least 20 values.
        if len(data) < 20:
            print(f"Skipping {site} because it only has {len(data)} values.")
            continue

        # Check to make sure the period of record covers at least 25 years.
        year0 = data[0][3]
        yearx = data[-1][3]
        if yearx - year0 < 25:
            print(f"Skipping {site} because it only covers {yearx - year0} years.")
            continue

        for lat, lon, elev, year, vpd in data:
            yield (site, lat, lon, elev, year, vpd)

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
    A list that contains  tuples of (average latitude, average longitude, average elevation, year,
    average vpd) where the average is for the requested month(s). Latitude, longitude, and 
    elevation should rarely change; only in years when the station moved.
    """
    if isinstance(months, Iterable):
        months = tuple(months)
    else:
        months = (months,)

    # Heuristic test to see if too much missing data.
    avg_hours = len(months) * 30 * 24
    min_allowed_hours = int(round(0.9 * avg_hours))

    data_gen = archive.get_hourly_vpd_data_for(site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((lat, lon, elev, vt - timedelta(hours=break_hour), vpd) 
            for _, vt, vpd, lat, lon, elev, _num_obs in data_gen)
    # Filter out data that isn't for the month(s) we're interested in.
    data_gen = ((lat, lon, elev, vt, vpd) 
            for lat, lon, elev, vt, vpd in data_gen if vt.month in months)

    results = []

    # Initialize the first time step and the accumulator variables
    try:
        lat, lon, elev, vt, vpd = next(data_gen)
    except StopIteration:
        return results
    curr_year = vt.year
    sum_vpd = vpd
    sum_lat = lat
    sum_lon = lon
    sum_elev = elev
    count = 1

    for lat, lon, elev, vt, vpd in data_gen:

        if curr_year != vt.year:

            # Save the values for the current year
            if count >= min_allowed_hours:
                results.append(
                        (sum_lat / count,
                         sum_lon / count,
                         sum_elev / count,
                         curr_year,
                         sum_vpd / count,
                        )
                )

            # Update the current year and reset the accumulator variables
            curr_year = vt.year
            sum_vpd = 0.0
            sum_lat = 0.0
            sum_lon = 0.0
            sum_elev = 0.0
            count = 0

        # Always update the accumulator variables!
        sum_vpd += vpd
        sum_lat += lat
        sum_lon += lon
        sum_elev += elev
        count += 1

    # Don't forget to yield the last set of values!
    if count >= min_allowed_hours:
        results.append(
                (sum_lat / count,
                 sum_lon / count,
                 sum_elev / count,
                 curr_year,
                 sum_vpd / count,
                )
        )

    return results

def get_annual_hours_avg_max_vpd_all_sites(archive, hours=1000, starting=None, ending=None, break_hour=6):
    """Get the annual maximum 'hours' average VPD for all sites.

    Averaging is done on the hourly observations.

    # Arguments
    archive - an archive from which to retrieve the data.
    hours - the length of the period to average.
    starting - when to start the time series.
    end - when to end the time series.
    break_hour - subtract this many hours from the valid time to assign the day. This means the 
                 early morning observations will go with the previous day.

    # Returns
    A generator object that yields a tuple of (site id, average latitude, average longitude, 
    average elevation, year, average vpd) where the average is for the
    requested month(s). Latitude, longitude, and elevation should rarely change; only in years when
    the station moved.
    """
    sites = archive.get_all_sites()

    for site in sites:

        data = get_annual_hours_avg_max_vpd(archive, site, hours, starting, ending, break_hour)

        # Check for at least 20 values.
        if len(data) < 20:
            print(f"Skipping {site} because it only has {len(data)} values.")
            continue

        # Check to make sure the period of record covers at least 25 years.
        year0 = data[0][3]
        yearx = data[-1][3]
        if yearx - year0 < 25:
            print(f"Skipping {site} because it only covers {yearx - year0} years.")
            continue

        for lat, lon, elev, year, vpd in data:
            yield (site, lat, lon, elev, year, vpd)

    return


def get_annual_hours_avg_max_vpd(archive, site, hours=1000, starting=None, ending=None, break_hour=6):
    """Get the annual maximum 'hours' average VPD.

    Averaging is done on the hourly observations.

    # Arguments
    archive - an archive from which to retrieve the data.
    site - the site id, e.g. 'kmso', 'smtm8'.
    hours - the length of the period to average.
    starting - when to start the time series.
    end - when to end the time series.
    break_hour - subtract this many hours from the valid time to assign the day. This means the 
                 early morning observations will go with the previous day.

    # Returns
    A list that contains  tuples of (average latitude, average longitude, average elevation, year,
    max average vpd) where the average has 'hours' period. Latitude, longitude, and 
    elevation should rarely change; only in years when the station moved.
    """
    data_gen = archive.get_xhourly_average_vpd_data_for(hours, site, starting, ending)

    # map the valid time with the break_hour
    data_gen = ((lat, lon, elev, vt - timedelta(hours=break_hour), vpd) 
            for _, vt, vpd, lat, lon, elev in data_gen)

    results = []

    # Initialize the first time step and the accumulator variables
    try:
        lat, lon, elev, vt, vpd = next(data_gen)
    except StopIteration:
        return results
    curr_year = vt.year
    max_avg_vpd = vpd
    avg_lat = lat
    avg_lon = lon
    avg_elev = elev
    count = 0

    # Need at least 90% of the hourly data.
    min_count = int(0.9 * 365 * 24)

    for lat, lon, elev, vt, vpd in data_gen:

        if curr_year != vt.year:

            # Save the values for the current year
            if max_avg_vpd is not None and count > min_count:
                results.append(
                        (avg_lat,
                         avg_lon,
                         avg_elev,
                         curr_year,
                         max_avg_vpd,
                        )
                )

            # Update the current year and reset the accumulator variables
            curr_year = vt.year
            max_avg_vpd = None
            avg_lat = None
            avg_lon = None
            avg_elev = None
            count = 0

        # Always update the accumulator variables!
        if max_avg_vpd is None or vpd > max_avg_vpd:
            avg_lat = lat
            avg_lon = lon
            avg_elev = elev
            max_avg_vpd = vpd
        count += 1

    # Don't forget to yield the last set of values!
    if max_avg_vpd is not None and count > min_count:
        results.append(
                (avg_lat,
                 avg_lon,
                 avg_elev,
                 curr_year,
                 max_avg_vpd,
                )
        )

    return results

class Normalization:
    '''A collection of data and functions related to normalizing a variable.'''

    def __init__(self, shift, scale):
        
        self.shift = shift
        self.scale = scale

    def denorm_scale(self, x):
        return x * self.scale

    def denorm(self, x):
        return x * self.scale + self.shift

    def norm(self, x):
        return (x - self.shift) / self.scale

    @staticmethod
    def denorm_slope(x_norm, y_norm, slope):
        return y_norm.denorm_scale(slope) / x_norm.denorm_scale(1.0)


def normalize_var(df, input_col, output_col):
    """Normalize a Pandas DataFrame column and add it back into the DataFrame as a new column.

    Normalize means subtract the mean and then divide by the standard deviation.

    # Arguments
    df - the Pandas DataFrame
    input_col - The name of the column in df to be normalized.
    output_col - The name of the column to add back to the DataFrame df with the normalized data.

    # Returns a Normalization object.
    """
    offset = df[input_col].mean()
    scale = df[input_col].std()
    
    df[output_col] = (df[input_col] - offset) / scale
    
    return Normalization(offset, scale)

def scale_var(df, input_col, output_col):
    """Normalize a Pandas DataFrame column and add it back into the DataFrame as a new column.

    Normalize means divide by the mean in this case.

    # Arguments
    df - the Pandas DataFrame
    input_col - The name of the column in df to be normalized.
    output_col - The name of the column to add back to the DataFrame df with the normalized data.

    # Returns a Normalization object.
    """
    offset = 0.0
    scale = df[input_col].mean()
    
    df[output_col] = df[input_col] / scale
    
    return Normalization(offset, scale)

def save_project_data(file_name, data_dictionary):
    '''Save the data to disk.

    # Arguments
    data_dictionary has all the objects you want to save.
    '''
    with open(Path('intermediate_data') / file_name, 'wb') as f:
        pickle.dump(data_dictionary, f, pickle.HIGHEST_PROTOCOL)
    return


def load_project_data(file_name):
    '''Load the project data.'''
    with open(Path('intermediate_data') / file_name, 'rb') as f:
        return pickle.load(f)


def save_inference_data(file_name, inference_data):
    inference_data.to_netcdf(Path('intermediate_data') / file_name)

    return


def load_inference_data(file_name):
    return az.InferenceData.from_netcdf(Path('intermediate_data') / file_name)


def plot_map(lats, lons, data=None, label_points=False, colormap=None, title=None, 
             colorbar_label=None, color_min=None, color_max=None, extents=None,
             figsize=None):
    
    # Set up the map background
    stamen_terrain = cimgt.Stamen('terrain-background')
    
    if figsize is None:
        figsize=(12, 12)
    
    fig, ax_map = plt.subplots(1, 1, figsize=figsize,
            subplot_kw={'projection':stamen_terrain.crs})
    
    if extents is None:
        extents = [-127.5, -100, 28, 51]
    
    ax_map.set_extent(extents, crs=ccrs.Geodetic())
        
    ax_map.add_image(stamen_terrain, 10)
    
    # Filter out data outside the extents
    if data is not None:
        min_lon, max_lon, min_lat, max_lat = extents
        filtered_data = []
        for lat, lon, d in zip(lats, lons, data):
            if lat > min_lat and lat < max_lat and lon > min_lon and lon < max_lon:
                filtered_data.append((lat, lon, d))
        lats, lons, data = zip(*tuple(filtered_data))

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    ax_map.add_feature(cfeature.COASTLINE)
    ax_map.add_feature(states_provinces, edgecolor='white')
    ax_map.add_feature(cfeature.BORDERS)
    
    if data is not None and not label_points:
        sc = ax_map.scatter(
            lons, 
            lats, 
            vmin=color_min, 
            vmax=color_max, 
            c=data, 
            cmap=colormap,
            edgecolors='black',
            linewidths=1.5,
            transform=ccrs.PlateCarree(),
            zorder=10,
        )
        plt.colorbar(sc, orientation='horizontal', fraction=0.06, label=colorbar_label)
    else:
        ax_map.scatter(
            lons, 
            lats, 
            transform=ccrs.PlateCarree(), 
            edgecolors='black',
            linewidths=1.5,
            color='white',
            zorder=10
        )

    if title is not None:
        fig.suptitle(title)
        
    if label_points:
        for lon, lat, label in zip(lons, lats, data):
            ax_map.annotate(label, (lon, lat), transform=ccrs.Geodetic(), zorder=12)
            
    ax_map.set_aspect(1)
    ax_map.set_adjustable('datalim')
    
    plt.show();
    
    return


def plot_station_correlation(dataframe, station1, station2, x_denorm, y_denorm, title=None):
    '''Make plots to show how two stations correlate with each other.'''

    fig, (ax_ts, ax_corr) = plt.subplots(2, 1, figsize=(9, 9))
    
    if title is not None:
        fig.suptitle(title)
    
    # The time series plot
    ax_ts.scatter(x_denorm(dataframe.loc[dataframe['site']==station1]['x_obs']), 
                   y_denorm(dataframe.loc[dataframe['site']==station1]['y_vpd_obs']), 
                   alpha=0.6,
                   label=station1)
    
    ax_ts.scatter(x_denorm(dataframe.loc[dataframe['site']==station2]['x_obs']), 
                       y_denorm(dataframe.loc[dataframe['site']==station2]['y_vpd_obs']), 
                       alpha=0.6,
                       label=station2)
    
    ax_ts.set_ylabel('Vapor Pressure\nDeficit (hPa)')
    ax_ts.legend()

    ax_ts.set_xlabel('Year')
    
    # The correllation plot
    group_a = dataframe.loc[dataframe['site']==station1]
    group_b = dataframe.loc[dataframe['site']==station2]
    joint_years = set(group_a['year']) & set(group_b['year'])

    xs = []
    ys = []
    
    for row in joint_years:
        xs.append((group_a.loc[group_a['year'] == row])['avg_vpd'].iloc[0])
        ys.append((group_b.loc[group_b['year'] == row])['avg_vpd'].iloc[0])

    max_val = max(max(xs), max(ys))
    min_val = min(min(xs), min(ys))

    ax_corr.plot([min_val-1, max_val+1],[min_val-1, max_val+1], color='black')
    ax_corr.scatter(xs, ys)
    
    ax_corr.set_xlabel(f"{station1} Vapor Pressure Deficit (hPa)")
    ax_corr.set_ylabel(f"{station2} Vapor Pressure Deficit (hPa)")
   
    ax_corr.annotate(f"Corr: {np.corrcoef(xs, ys)[0][1]:.3f}", (min_val, max_val))
    ax_corr.annotate(f"Elevations:\n{station1} {group_a.iloc[0]['elev']:.1f}m\n{station2} {group_b.iloc[0]['elev']:.1f}m",
                     (max_val-2, min_val) )

    plt.show()
    
    return
