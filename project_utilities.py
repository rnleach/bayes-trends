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

def get_coords_for_sites(archive, sites):
    """Given a list of sites, get the lat, lon, and elevation information to go with it."""
    results = []
    for site in sites:
        loc = archive.get_location_info_for_site(site)
        if loc is None:
            raise Exception(f"This site did not come from the archive {site}.")
        lat, lon, elev = loc
        results.append((site, lat, lon, elev))

    return tuple(results)


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


def get_hourly_vpd_fill_missing(archive, site, starting=None, ending=None):
    """Get an hourly time series of vapor pressure deficit with missing hours filled with NaNs.

    # Arguments
    archive - an archive from which to retrieve the data.
    site - the site id, e.g. 'kmso', 'smtm8'.
    starting - when to start the time series.
    end - when to end the time series.

    # Returns
    A generator that yields tuples of (hour since starting, vpd).
    """
    data_gen = archive.get_hourly_vpd_data_for(site, starting, ending)

    # Initialize the first time step and the accumulator variables
    try:
        site_ , vt_prev, vpd_prev, lat_prev_, lon_prev_, elev_prev_, num_matching_ = next(data_gen)
    except StopIteration:
        return 

    if starting is None:
        starting = vt_prev
    
    for site_, vt, vpd, lat_, lon_, elev_, num_matching_ in data_gen:
        total_hours = (vt_prev - starting).total_seconds() // 3600
        yield (total_hours, vpd_prev)

        while (vt - vt_prev).total_seconds() > 3600:
            vt_prev += timedelta(hours=1)
            total_hours += 1
            yield (total_hours, float("nan"))

        vt_prev, vpd_prev = vt, vpd

    yield (total_hours + 1, vpd_prev)

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

def get_annual_hours_avg_max_vpd_all_sites(archive, starting, ending, hours=1000, break_hour=6):
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

    end_year = int(ending.year)

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

        curr_year = int(starting.year)

        for lat_, lon_, elev_, year, vpd in data:

            while curr_year < int(year):
                yield (site, curr_year, float("nan"))
                curr_year += 1

            yield (site, year, vpd)

            curr_year += 1

        while curr_year <= end_year:
            yield(site, curr_year, float("nan"))
            curr_year += 1

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

        return

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

def shift_and_scale_var(df, input_col, output_col, shift, scale):
    """Normalize a Pandas DataFrame column and add it back into the DataFrame as a new column.

    This function creates a normalization that subtracts the provided offset and then divides by
    the provided scale.

    # Arguments
    df - the Pandas DataFrame
    input_col - The name of the column in df to be normalized.
    output_col - The name of the column to add back to the DataFrame df with the normalized data.
    shift - the value to shift the column by.
    scale - the value to scale the column by.

    # Returns a Normalization object.
    """
    df[output_col] = (df[input_col] - shift) / scale
    
    return Normalization(shift, scale)


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


def haversine_distance(lon1, lat1, lon2, lat2):
    """Calculate the great circle distance between two points."""

    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367.4445 * c
    return km

def common_cases(dataframe, station1, station2):
    """Return VPD from the 2 stations only for the years they have in common."""
    group_a = dataframe.loc[dataframe['site']==station1]
    group_b = dataframe.loc[dataframe['site']==station2]
    joint_years = set(group_a['year']) & set(group_b['year'])

    xs = []
    ys = []
    years = []
    
    for row in joint_years:
        years.append(row)
        xs.append((group_a.loc[group_a['year'] == row])['avg_vpd'].iloc[0])
        ys.append((group_b.loc[group_b['year'] == row])['avg_vpd'].iloc[0])

    return (xs, ys, years)


def station_corr_cov(dataframe, station1, station2):
    """Calculate the correlation between the values at 2 stations."""
    xs, ys, _ = common_cases(dataframe, station1, station2)
    cor = np.corrcoef(xs, ys)[0][1]
    cov = np.cov(xs, ys)[0][1]
    
    return (cor, cov)

def distance_elevation_corr_cov_relationships(site_coords, df):
    num_sites = len(site_coords)

    dff = df.dropna()

    distance_matrix = np.zeros((num_sites, num_sites))
    elevation_matrix = np.zeros((num_sites, num_sites))

    distances = []
    elevations = []
    corrs = []
    covs = []

    for i in range(num_sites):

        station1 = site_coords.iloc[i]['site']
        lat1 = site_coords.iloc[i]['lat']
        lon1 = site_coords.iloc[i]['lon']
        elev1 = site_coords.iloc[i]['elev']

        for j in range(i, num_sites):

            station2 = site_coords.iloc[j]['site']
            lat2 = site_coords.iloc[j]['lat']
            lon2 = site_coords.iloc[j]['lon']
            elev2 = site_coords.iloc[j]['elev']

            corr, cov = station_corr_cov(dff, station1, station2)

            corrs.append(corr)
            covs.append(cov)

            if i == j:
                distance = 0.0
                elev_diff = 0.0
            else:
                distance = haversine_distance(lon1, lat1, lon2, lat2)
                elev_diff = abs(elev2 - elev1)

            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

            elevation_matrix[i, j] = elev_diff
            elevation_matrix[j, i] = elev_diff

            distances.append(distance)
            elevations.append(elev_diff)
            
    return (distance_matrix, elevation_matrix, distances, elevations, corrs, covs)


def plot_station_correlation(df, station1, station2, x_denorm, y_denorm, title=None):
    '''Make plots to show how two stations correlate with each other.'''

    dataframe = df.dropna()

    fig, ((ax_ts, ax_corr), (ax_ts_log, ax_corr_log)) = plt.subplots(2, 2, figsize=(12, 9))
    
    if title is not None:
        fig.suptitle(title)
    
    xs1 = dataframe.loc[dataframe['site']==station1]['x_obs']
    ys1 = dataframe.loc[dataframe['site']==station1]['y_vpd_obs']

    xs2 = dataframe.loc[dataframe['site']==station2]['x_obs']
    ys2 = dataframe.loc[dataframe['site']==station2]['y_vpd_obs']

    # The time series plots
    ax_ts.scatter(x_denorm(xs1), y_denorm(ys1), alpha=0.6, label=station1)
    ax_ts.scatter(x_denorm(xs2), y_denorm(ys2), alpha=0.6, label=station2)
    ax_ts.set_ylabel('Vapor Pressure\nDeficit (hPa)')
    ax_ts.legend()
    ax_ts.set_xlabel('Year')
    
    ax_ts_log.scatter(xs1, np.log(ys1), alpha=0.6, label=station1)
    ax_ts_log.scatter(xs2, np.log(ys2), alpha=0.6, label=station2)
    ax_ts_log.set_ylabel('Log Normalized VPD')
    ax_ts_log.legend()
    ax_ts_log.set_xlabel('Year')
    
    # The correllation plot
    xs, ys, _ = common_cases(dataframe, station1, station2)

    max_val = max(max(xs), max(ys))
    min_val = min(min(xs), min(ys))

    ax_corr.plot([min_val-1, max_val+1],[min_val-1, max_val+1], color='black')
    ax_corr.scatter(xs, ys)
    
    ax_corr.set_xlabel(f"{station1} Vapor Pressure Deficit (hPa)")
    ax_corr.set_ylabel(f"{station2} Vapor Pressure Deficit (hPa)")
   
    ax_corr.annotate(f"Corr: {np.corrcoef(xs, ys)[0][1]:.3f}", (min_val, max_val))

    ax_corr_log.plot(
            [np.log(min_val)-1, np.log(max_val)+1],
            [np.log(min_val)-1, np.log(max_val)+1],
            color='black')
    ax_corr_log.scatter(np.log(xs), np.log(ys))
    
    ax_corr_log.set_xlabel(f"{station1} Log Normalized VPD")
    ax_corr_log.set_ylabel(f"{station2} Log Normalized VPD")
   
    ax_corr_log.annotate(f"Corr: {np.corrcoef(np.log(xs), np.log(ys))[0][1]:.3f}",
            (np.log(min_val), np.log(max_val)))


    plt.show()
    
    return
