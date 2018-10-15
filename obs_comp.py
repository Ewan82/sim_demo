import numpy as np
import signaturesimulator as ss
import matplotlib.pyplot as plt
import matplotlib.mlab as malb
import netCDF4 as nc



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def extract_dates(nc_file):
    nc_dat = nc.Dataset(nc_file, 'r')
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    nc_dat.close()
    return date


def extract_all_refl(nc_file, lat, lon):
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_idx = find_nearest(nc_dat.variables['y'][:], lat)[1]
    lon_idx = find_nearest(nc_dat.variables['x'][:], lon)[1]
    refl = nc_dat.variables['data'][:, :, lat_idx, lon_idx]
    nc_dat.close()
    return refl


def extract_band_refl(nc_file, band, lat, lon):
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_idx = find_nearest(nc_dat.variables['y'][:], lat)[1]
    lon_idx = find_nearest(nc_dat.variables['x'][:], lon)[1]
    band_refl = nc_dat.variables['data'][:, band, lat_idx, lon_idx]
    nc_dat.close()
    return band_refl


def extract_s1_geoms(nc_file, lat, lon, fname=None, time_tuple=None):
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_val, lat_idx = find_nearest(nc_dat.variables['lat'][:], lat)
    lon_val, lon_idx = find_nearest(nc_dat.variables['lon'][:], lon)
    vza = nc_dat.variables['localIncidenceAngle'][:, lon_idx, lat_idx]
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    date_str = [a.strftime("%Y/%m/%d %H:%M") for a in date]
    nc_dat.close()
    if fname is None:
        ret_val = vza
    elif fname is not None and time_tuple is None:
        geom = np.array([date_str, vza, np.zeros(len(vza)), np.zeros(len(vza)), np.zeros(len(vza)),
                         np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza))]).T
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, lat, lon')
        ret_val = 'geom_file saved!'
    else:
        t_idx1 = nc.date2index(time_tuple[0], nc_dat.variables['time'], select='before')
        t_idx2 = nc.date2index(time_tuple[1], nc_dat.variables['time'], select='after')
        date_str = date_str[t_idx1:t_idx2]
        vza = vza[t_idx1:t_idx2]
        geom = np.array([date_str, vza, np.zeros(len(vza)), np.zeros(len(vza)), np.zeros(len(vza)),
                         np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza))]).T
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, lat, lon')
        ret_val = 'geom_file saved!'
    return ret_val


def extract_s2_geoms(nc_file, lat, lon, fname=None, time_tuple=None):
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_val, lat_idx = find_nearest(nc_dat.variables['y'][:], lat)
    lon_val, lon_idx = find_nearest(nc_dat.variables['x'][:], lon)
    vza = nc_dat.variables['VZA_MEAN'][:]
    vaa = nc_dat.variables['VAA_MEAN'][:]
    sza = nc_dat.variables['SZA'][:]
    saa = nc_dat.variables['SAA'][:]
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    date_str = [a.strftime("%Y/%m/%d %H:%M") for a in date]
    nc_dat.close()
    if fname is None:
        ret_val = vza
    elif fname is not None and time_tuple is None:
        geom = np.array([date_str, vza, vaa, sza, saa, np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza))]).T
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, lat, lon')
        ret_val = 'geom_file saved!'
    else:
        t_idx1 = nc.date2index(time_tuple[0], nc_dat.variables['time'], select='before')
        t_idx2 = nc.date2index(time_tuple[1], nc_dat.variables['time'], select='after')
        date_str = date_str[t_idx1:t_idx2]
        vza = vza[t_idx1:t_idx2]
        vaa = vaa[t_idx1:t_idx2]
        sza = sza[t_idx1:t_idx2]
        saa = saa[t_idx1:t_idx2]
        geom = np.array([date_str, vza, vaa, sza, saa, np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza))]).T
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, lat, lon')
        ret_val = 'geom_file saved!'
    return ret_val


def extract_backscat(nc_file, pol, lat, lon):
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_idx = find_nearest(nc_dat.variables['lat'][:], lat)[1]
    lon_idx = find_nearest(nc_dat.variables['lon'][:], lon)[1]
    backscat = nc_dat.variables[pol][:, lat_idx, lon_idx]
    nc_dat.close()
    return backscat


def model_band_refl(site_nml, state_csv, geom_csv):
    sim = ss.Simulator(site_nml=site_nml)
    sim.get_land_state = sim.state_csv(fname=state_csv)
    sim.get_geom = sim.geom_csv(fname=geom_csv)
    sim.run_rt = sim.passive_optical
    sim.run()
    return np.array(sim.spectra.refl)


def model_backscat(site_nml, state_csv, geom_csv, pol='vv'):
    sim = ss.Simulator(site_nml=site_nml)
    sim.get_land_state = sim.state_csv(fname=state_csv)
    sim.get_geom = sim.geom_csv(fname=geom_csv)
    sim.run_rt = sim.active_microwave
    sim.run()
    if pol is 'vv':
        ret_val = np.array(sim.backscat.vv)
    elif pol is 'hv':
        ret_val = np.array(sim.backscat.hv)
    elif pol is 'hh':
        ret_val = np.array(sim.backscat.hh)
    else:
        ret_val = None
        print 'Invalid polarisation specified'
    return ret_val


def plot_ndvi(band4_refl, band8a_refl, dates=None, marker='o', label='sentinel ndvi'):
    plt.plot(dates, (band8a_refl-band4_refl)/(band8a_refl+band4_refl), marker, label=label)
    plt.ylabel('NDVI')
    plt.xlabel('Date')
    return 'plot made!'
