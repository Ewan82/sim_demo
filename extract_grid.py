import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.animation as animation
import datetime as dt
import seaborn as sns
import pandas as pd
import joypy


def find_nearest(array, value):
    """
    Find nearest value in an array
    :param array: array of values
    :param value: value for which to find nearest element
    :return: nearest value in array, index of nearest value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def draw_map(low_lat=48.246779775111911, high_lat=48.253191859919177, low_lon=11.711720836040126,
             high_lon=11.724213392466636, res=3.42892e-5/2, ax=None):
    """
    Creates a cylindrical Basemap instance.
    :param low_lat: lower left lat
    :param high_lat: upper right lat
    :param low_lon: lower left lon
    :param high_lon: upper right lon
    :param ax: axis to create instance for
    :return: Basemap instance
    """
    if ax is None:
        m = Basemap(projection='merc', resolution='i',
                    llcrnrlat=low_lat-res, urcrnrlat=high_lat+res,
                    llcrnrlon=low_lon-res, urcrnrlon=high_lon+res, lat_0=48.2497, lon_0=11.7179,
                    suppress_ticks = True)
    else:
        m = Basemap(projection='cyl',resolution='i',
                    llcrnrlat=low_lat - res, urcrnrlat=high_lat + res,
                    llcrnrlon=low_lon - res, urcrnrlon=high_lon + res, ax=ax)
    # draw coastlines, state and country boundaries, edge of map.
    #m.drawcoastlines()
    #m.drawstates()
    #m.drawcountries()
    parallels = np.linspace(low_lat, high_lat, 5.0) #np.arange(0., 81, 0.005)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[True, False, False, False], dashes=[6,900], color='w')
    meridians = np.linspace(low_lon, high_lon, 5)
    m.drawmeridians(meridians, labels=[False, False, False, True], dashes=[6,900], color='w')
    return m


def draw_map_sub(low_lat=4., high_lat=12., low_lon=-9., high_lon=2., ax=None):
    """
    Creates a cylindrical Basemap instance.
    :param low_lat: lower left lat
    :param high_lat: upper right lat
    :param low_lon: lower left lon
    :param high_lon: upper right lon
    :param ax: axis to create instance for
    :return: Basemap instance
    """
    if ax is None:
        m = Basemap(projection='cyl',resolution='i',
                    llcrnrlat=low_lat, urcrnrlat=high_lat,
                    llcrnrlon=low_lon, urcrnrlon=high_lon)
    else:
        m = Basemap(projection='cyl',resolution='i',
                    llcrnrlat=low_lat, urcrnrlat=high_lat,
                    llcrnrlon=low_lon, urcrnrlon=high_lon, ax=ax)
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    parallels = np.arange(0., 81, 5.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False], labelstyle='+/-', fontsize=12)
    meridians = np.arange(0., 361., 5.)
    m.drawmeridians(meridians, labels=[True, False, False, True], labelstyle='+/-', fontsize=12)
    return m


def map_plt(arr, colormap='OrRd', ax=None, v_min=None, v_max=None, points=None):
    # draw map
    arr[arr < -900] = np.nan
    m = draw_map(ax=ax)
    cs = m.imshow(arr, cmap=colormap, interpolation='nearest', vmin=v_min, vmax=v_max)
    if points is not None:
        for p in points:
            m.plot(p[0], p[1], 'o', color='red', markersize=4)
    return ax, cs, m


def map_cbar(arr, cbarmap='Blues', axes=None, title=None):
    """
    Plots a map of error between CCI obs and JULES for either prior or posterior
    :param modelled_sm: array of modelled soil moisture values
    :param colormap: choice of colormap as string
    :return: figure
    """
    if axes is not None:
        ax = axes
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1)
    # plot data
    ax, cs, m = map_plt(arr, colormap=cbarmap, ax=ax, v_min=0, v_max=1)
    # add colorbar.
    print np.nanmean(arr)
    if axes is None:
        cbar = m.colorbar(cs, location='bottom', pad="5%")
        #clevs = [-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15,]
        #cbar.set_ticks(clevs,)
        #cbar.ax.set_xticklabels(clevs, rotation=45)
        #cbar.set_label(r'Bias (m$^3$ m$^{-3}$)')
        if title is not None:
            ax.set_title(title)
        ret_val = fig
    else:
        ret_val = ax, cs
    return ret_val


def extract_s2_data(nc_file, hi_lat=48.253205, lo_lat=48.246775, hi_lon=11.724235, lo_lon=11.711744,
                    time=dt.datetime(2017, 8, 3, 0, 0)):
    """
    Extract S2 data from specified NetCDF file
    :param nc_file: filename of netCDF file
    :param lat: latitude for which to extract S2 data
    :param lon: longitude for which to extract S2 data
    :param time: time to extracte S2 data for
    :return:
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    hi_lat_idx = find_nearest(nc_dat.variables['y'][:], hi_lat)[1]
    print hi_lat_idx
    lo_lat_idx = find_nearest(nc_dat.variables['y'][:], lo_lat)[1]
    print lo_lat_idx
    hi_lon_idx = find_nearest(nc_dat.variables['x'][:], hi_lon)[1]
    print hi_lon_idx
    lo_lon_idx = find_nearest(nc_dat.variables['x'][:], lo_lon)[1]
    print lo_lon_idx
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    date_idx = find_nearest(date, time)[1]
    print date_idx
    band_dat = nc_dat.variables['data'][date_idx, :, hi_lat_idx:lo_lat_idx, lo_lon_idx:hi_lon_idx]
    lats = nc_dat.variables['y'][hi_lat_idx:lo_lat_idx]
    lons = nc_dat.variables['x'][lo_lon_idx:hi_lon_idx]
    nc_dat.close()
    return band_dat, lats, lons


def extract_s1_data(nc_file, hi_lat=48.253205, lo_lat=48.246775, hi_lon=11.724235, lo_lon=11.711744,
                    time=dt.datetime(2017, 8, 3, 0, 0)):
    """
    Extract S2 data from specified NetCDF file
    :param nc_file: filename of netCDF file
    :param lat: latitude for which to extract S2 data
    :param lon: longitude for which to extract S2 data
    :param time: time to extracte S2 data for
    :return:
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    hi_lat_idx = find_nearest(nc_dat.variables['lat'][:], hi_lat)[1]
    print hi_lat_idx
    lo_lat_idx = find_nearest(nc_dat.variables['lat'][:], lo_lat)[1]
    print lo_lat_idx
    hi_lon_idx = find_nearest(nc_dat.variables['lon'][:], hi_lon)[1]
    print hi_lon_idx
    lo_lon_idx = find_nearest(nc_dat.variables['lon'][:], lo_lon)[1]
    print lo_lon_idx
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    date_idx = find_nearest(date, time)[1]
    print date_idx
    vh = nc_dat.variables['sigma0_vh_multi'][date_idx, hi_lat_idx:lo_lat_idx, lo_lon_idx:hi_lon_idx]
    vv = nc_dat.variables['sigma0_vv_multi'][date_idx, hi_lat_idx:lo_lat_idx, lo_lon_idx:hi_lon_idx]
    lats = nc_dat.variables['lat'][hi_lat_idx:lo_lat_idx]
    lons = nc_dat.variables['lon'][lo_lon_idx:hi_lon_idx]
    nc_dat.close()
    return vh, vv, lats, lons


def plot_ndvi(nc_file='/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc'):
    refl, lats, lons = extract_s2_data(nc_file=nc_file, hi_lat=48.253457, lo_lat=48.245823, lo_lon=11.711185,
                                       hi_lon=11.724635, time=dt.datetime(2017,6,26,0,0))
    ndvi = (refl[7, :, :]-refl[3, :, :]) / (refl[7, :, :]+refl[3, :, :])
    print ndvi.shape
    #fig, ax = plt.subplots(nrows=1, ncols=1)
    map = draw_map(low_lat=lats[-1], high_lat=lats[0], low_lon=lons[0], high_lon=lons[-1])
    lon, lat = np.meshgrid(lons, lats)
    xx, yy = map(lon, lat)
    palette = sns.color_palette("colorblind", 11)
    map.pcolormesh(xx, yy, ndvi, cmap='viridis', vmax=1, vmin=0)
    x,y = map(11.717295, 48.251387)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k')
    plt.text(x+40, y+20, 'high')
    x,y = map(11.719048, 48.250805)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k')
    plt.text(x+40, y+20, 'mid')
    x,y = map(11.719875, 48.249410)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k', label='Field sampling points')
    plt.text(x+40, y+20, 'low')
    x,y = map(11.720441, 48.247624)
    map.plot(x, y, '*', color=palette[3], markersize=10, markeredgewidth=0.5, markeredgecolor='k', label='Farmyard retrieval site')
    #x,y = map(11.715308, 48.249028)
    #map.plot(x, y, '*', color=palette[1], markersize=10, markeredgewidth=0.5, markeredgecolor='k', label='Forest retrieval site')
    x, y = map(11.7211, 48.24936)
    map.plot(x, y, '*', color=palette[1], markersize=10, markeredgewidth=0.5, markeredgecolor='k',
             label='Gravel path retrieval site')
    plt.legend(fancybox=True, frameon=True)
    plt.title('NDVI Munich field site 26/06/2017')
    plt.colorbar()
    #map.imshow(ndvi, cmap='viridis', vmax=1, vmin=0)
    #locs = np.arange(0, len(lons), len(lons)/6)
    plt.show()
    return 'd'


def plot_ndvi_ita(nc_file='/export/cloud/nceo/users/if910917/sentinel_data/S2/italy/ITA_181001.nc'):
    refl, lats, lons = extract_s2_data(nc_file=nc_file, hi_lat=41.3798, lo_lat=41.3588, lo_lon=15.4751,
                                       hi_lon=15.5045, time=dt.datetime(2017,3,26,0,0))
    ndvi = (refl[7, :, :]-refl[3, :, :]) / (refl[7, :, :]+refl[3, :, :])
    print ndvi.shape
    #fig, ax = plt.subplots(nrows=1, ncols=1)
    map = draw_map(low_lat=lats[-1], high_lat=lats[0], low_lon=lons[0], high_lon=lons[-1])
    lon, lat = np.meshgrid(lons, lats)
    xx, yy = map(lon, lat)
    palette = sns.color_palette("colorblind", 11)
    map.pcolormesh(xx, yy, ndvi, cmap='viridis', vmax=1, vmin=0)
    x,y = map(15.47727, 41.36014)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k')
    plt.text(x+40, y+20, 'sg01')
    x,y = map(15.486868, 41.37548)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k')
    plt.text(x+40, y+20, 'sg04')
    x,y = map(15.500821, 41.371676)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k', label='Field sampling points')
    plt.text(x+40, y+20, 'sg11')
    plt.legend(fancybox=True, frameon=True)
    plt.title('NDVI Segezia field site 26/03/2017')
    plt.colorbar()
    #map.imshow(ndvi, cmap='viridis', vmax=1, vmin=0)
    #locs = np.arange(0, len(lons), len(lons)/6)
    plt.show()
    return 'd'


def plot_ndvi_pol(nc_file='/export/cloud/nceo/users/if910917/sentinel_data/S2/poland/POL_181001.nc'):
    refl, lats, lons = extract_s2_data(nc_file=nc_file, hi_lat=53.6373, lo_lat=53.6317, lo_lon=22.9769,
                                       hi_lon=22.9867, time=dt.datetime(2017,3,26,0,0))
    ndvi = (refl[7, :, :]-refl[3, :, :]) / (refl[7, :, :]+refl[3, :, :])
    print ndvi.shape
    #fig, ax = plt.subplots(nrows=1, ncols=1)
    map = draw_map(low_lat=lats[-1], high_lat=lats[0], low_lon=lons[0], high_lon=lons[-1])
    lon, lat = np.meshgrid(lons, lats)
    xx, yy = map(lon, lat)
    palette = sns.color_palette("colorblind", 11)
    map.pcolormesh(xx, yy, ndvi, cmap='viridis', vmax=1, vmin=0)
    x,y = map(22.98104, 53.63502)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k')
    plt.text(x+40, y+20, '1')
    x,y = map(22.98128, 53.63398)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k')
    plt.text(x+40, y+20, '2')
    x,y = map(22.98151, 53.63294)
    map.plot(x, y, '*', color='r', markersize=10, markeredgewidth=0.5, markeredgecolor='k', label='Field sampling points')
    plt.text(x+40, y+20, '3')
    plt.legend(fancybox=True, frameon=True)
    plt.title('NDVI Biebrza field site 26/03/2017')
    plt.colorbar()
    #map.imshow(ndvi, cmap='viridis', vmax=1, vmin=0)
    #locs = np.arange(0, len(lons), len(lons)/6)
    plt.show()
    return 'd'


def plot_backscat(nc_file='/export/cloud/nceo/users/if910917/sentinel_data/S1/munich/MNI_2017.nc'):
    vh, vv, lats, lons = extract_s1_data(nc_file=nc_file, hi_lat=48.253457, lo_lat=48.245823, lo_lon=11.711185,
                                       hi_lon=11.724635, time=dt.datetime(2017,3,8,0,0))
    map = draw_map(low_lat=lats[-1], high_lat=lats[0], low_lon=lons[0], high_lon=lons[-1], res=(lats[1]-lats[0])/2)
    lon, lat = np.meshgrid(lons, lats)
    xx, yy = map(lon, lat)
    map.pcolormesh(xx, yy, vv, cmap='viridis',)  # vmax=, vmin=0)
    plt.show()
    return 'done'