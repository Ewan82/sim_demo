import numpy as np
import matplotlib.mlab as mlab
import netCDF4 as nc
import glob
import utm
import datetime as dt



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


def extract_s1_data(nc_file, lat, lon, fname=None, time_tuple=None, remove_pol=None):
    """
    Extract S1 data from specified NetCDF file
    :param nc_file: filename of netCDF file
    :param lat: latitude for which to extract S1 data
    :param lon: longitude for which to extract S1 data
    :param fname: filename to save extracted S1 data under
    :param time_tuple: tupel containing two times (start and end) to extract S1 data between
    :return:
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_val, lat_idx = find_nearest(nc_dat.variables['lat'][:], lat)
    lon_val, lon_idx = find_nearest(nc_dat.variables['lon'][:], lon)
    vza = nc_dat.variables['localIncidenceAngle'][:, lon_idx, lat_idx]
    hv = nc_dat.variables['sigma0_vh_norm_multi_db'][:, lon_idx, lat_idx]
    vv = nc_dat.variables['sigma0_vv_norm_multi_db'][:, lon_idx, lat_idx]
    sat_flag = nc_dat.variables['satellite'][:]
    mission_lst = ['A', 'B']
    sat_flag = np.array(['S1'+mission_lst[int(x)] for x in sat_flag])
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    date_str = [a.strftime("%Y/%m/%d %H:%M") for a in date]
    if fname is None:
        ret_val = vza
    elif fname is not None and time_tuple is None:
        geom = np.array([date_str, vza, np.zeros(len(vza)), np.zeros(len(vza)), np.zeros(len(vza)),
                         sat_flag, np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza)), hv, vv
                         ]).T
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, sat_flag,lat, lon,  hv, vv')
        ret_val = 'geom_file saved!'
    else:
        t_idx1 = nc.date2index(time_tuple[0], nc_dat.variables['time'], select='before')
        t_idx2 = nc.date2index(time_tuple[1], nc_dat.variables['time'], select='after')
        date_str = date_str[t_idx1:t_idx2]
        vza = vza[t_idx1:t_idx2]
        hv = hv[t_idx1:t_idx2]
        vv = vv[t_idx1:t_idx2]
        sat_flag = sat_flag[t_idx1:t_idx2]
        geom = np.array([date_str, vza, np.zeros(len(vza)), np.zeros(len(vza)), np.zeros(len(vza)),
                         sat_flag, np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza)), hv, vv,
                         ]).T
        if remove_pol is not None:
            geom[:, remove_pol+8] = -99999
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, sat_flag, lat, lon, hv, vv')
        ret_val = 'geom_file saved!'
    nc_dat.close()
    return ret_val


def extract_s2_data(nc_file, lat, lon, fname=None, time_tuple=None, remove_bnds=(0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12)):
    """
    Extract S2 data from specified NetCDF file
    :param nc_file: filename of netCDF file
    :param lat: latitude for which to extract S2 data
    :param lon: longitude for which to extract S2 data
    :param fname: filename to save extracted S2 data under
    :param time_tuple: tupe containing two times (start and end) to extract S2 data between
    :return:
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    lat_val, lat_idx = find_nearest(nc_dat.variables['y'][:], lat)
    lon_val, lon_idx = find_nearest(nc_dat.variables['x'][:], lon)
    vza = nc_dat.variables['VZA_MEAN'][:]
    vaa = nc_dat.variables['VAA_MEAN'][:]
    sza = nc_dat.variables['SZA'][:]
    saa = nc_dat.variables['SAA'][:]
    mission = nc_dat.variables['SPACECRAFT_NAME'][:]
    mission = np.array(['S'+miss[-2:] for miss in mission])
    band_dat = nc_dat.variables['data'][:, :, lat_idx, lon_idx]
    date = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    date_str = [a.strftime("%Y/%m/%d %H:%M") for a in date]
    if time_tuple is None:
        geom = np.array([date_str, vza, vaa, sza, saa, mission, np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza)),
                         band_dat[:,0], band_dat[:,1], band_dat[:,2], band_dat[:,3],
                         band_dat[:, 4], band_dat[:,5], band_dat[:,6], band_dat[:,7], band_dat[:,8], band_dat[:,9],
                         band_dat[:, 10], band_dat[:,11], band_dat[:,12]]).T
        for idx in np.array(remove_bnds)+8:
            geom[:, idx] = -99999
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, sat_flag, lat, lon, b1, b2, b3, b4, b5'
                                                                ', b6, b7, b8a, b8b, b9, b10, b11, b12')
        ret_val = 'geom_file saved!'
    else:
        t_idx1 = nc.date2index(time_tuple[0], nc_dat.variables['time'], select='before')
        t_idx2 = nc.date2index(time_tuple[1], nc_dat.variables['time'], select='after')
        date_str = date_str[t_idx1:t_idx2]
        vza = vza[t_idx1:t_idx2]
        vaa = vaa[t_idx1:t_idx2]
        sza = sza[t_idx1:t_idx2]
        saa = saa[t_idx1:t_idx2]
        mission = mission[t_idx1:t_idx2]
        band_dat = band_dat[t_idx1:t_idx2, :]
        geom = np.array([date_str, vza, vaa, sza, saa, mission, np.array([lat_val]*len(vza)), np.array([lon_val]*len(vza)),
                         band_dat[:,0], band_dat[:,1], band_dat[:,2], band_dat[:,3],
                         band_dat[:, 4], band_dat[:,5], band_dat[:,6], band_dat[:,7], band_dat[:,8], band_dat[:,9],
                         band_dat[:, 10], band_dat[:,11], band_dat[:,12]]).T
        for idx in np.array(remove_bnds)+8:
            geom[:, idx] = -99999
        np.savetxt(fname, geom, fmt='%s', delimiter=',', header='date, vza, vaa, sza, saa, sat_flag, lat, lon, b1, b2, b3, b4, b5'
                                                                ', b6, b7, b8a, b8b, b9, b10, b11, b12')
        ret_val = 'geom_file saved!'
    nc_dat.close()
    return ret_val


def extract_state_arr_mean(state_crop_csv, state_sm_dir, point_no, fdir=None, point_level='high'):
    """
    Extract state array of field data from munich test sites.
    :param state_crop_csv: filename of crop field data csv
    :param state_sm_dir: directory location of soil moisture field data
    :param point_no: field point number, corresponding to munich test site numbering
    :param fdir: directory to save the output
    :param point_level: which point to use choice of: 'low', 'high', 'med', 'mean'
    :return:
    """
    state_dat = mlab.csv2rec(state_crop_csv, skiprows=1, delimiter=',')
    dates = state_dat['none']
    date_str = [a.strftime("%Y/%m/%d %H:%M") for a in dates]
    lai_high = state_dat['lai']
    lai_low = state_dat['lai_1']
    lai_med = state_dat['lai_2']
    lai_mean = state_dat['lai_mean_hants']
    lai_std = state_dat['lai_std']
    canht_high = state_dat['height_cm'] / 100.
    canht_low = state_dat['height_cm_1'] / 100.
    canht_med = state_dat['height_cm_2'] / 100.
    canht_mean = state_dat['height_cm_mean'] / 100.
    canht_std = state_dat['height_cm_std'] / 100.

    sm_loc_dat = mlab.csv2rec(state_sm_dir+'/locations_utm_epsg-32632.csv')
    easting = sm_loc_dat['point_x'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'med')]
    northing = sm_loc_dat['point_y'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'med')]
    file_head_1 = sm_loc_dat['esu_sm'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'med')][0]
    file_head_2 = sm_loc_dat['esu_sm'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'high')][0]
    file_head_3 = sm_loc_dat['esu_sm'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'low')][0]
    lat, lon = utm.to_latlon(easting, northing, 32, 'U')
    sm_csv1 = glob.glob(state_sm_dir + '/' + file_head_1 +'*SM.csv')[0]
    sm_dat1 = mlab.csv2rec(sm_csv1)
    sm_csv2 = glob.glob(state_sm_dir + '/' + file_head_2 + '*SM.csv')[0]
    sm_dat2 = mlab.csv2rec(sm_csv2)
    sm_csv3 = glob.glob(state_sm_dir + '/' + file_head_3 + '*SM.csv')[0]
    sm_dat3 = mlab.csv2rec(sm_csv3)

    #sm_idx = [find_nearest(sm_dat['date'], dt.datetime.combine(x,dt.datetime.min.time()))[1] for x in dates]
    sm_dates = sm_dat1['date'][3:]
    date_str_sm = [a.strftime("%Y/%m/%d %H:%M") for a in sm_dates]
    sm1_port1 = sm_dat1['port1_sm'][3:]
    sm1_port2 = sm_dat1['port2_sm'][3:]
    sm1_mean = np.nanmean((sm1_port1, sm1_port2), axis=0)
    sm2_port1 = sm_dat2['port1_sm']
    sm2_port2 = sm_dat2['port2_sm']
    sm2_mean = np.nanmean((sm2_port1, sm2_port2), axis=0)
    sm3_port1 = sm_dat3['port1_sm'][3:]
    sm3_port2 = sm_dat3['port2_sm'][3:]
    sm3_mean = np.nanmean((sm3_port1, sm3_port2), axis=0)
    sm_mean = np.nanmean((sm1_mean, sm2_mean, sm3_mean), axis=0)
    sm_std = np.nanstd((sm1_port1, sm1_port2, sm2_port1, sm2_port2, sm3_port1, sm3_port2), axis=0)

    save_arr_lai_canht = np.array([date_str, lai_mean, canht_mean, np.array([lat] * len(dates)),
                        np.array([lon] * len(dates)), lai_std, canht_std]).T
    fname_lai_canht = fdir + '/mni_lai_canht_field_' + str(point_no) + '_' + str(point_level) + '.csv'
    np.savetxt(fname_lai_canht, save_arr_lai_canht, fmt='%s', delimiter=',', header='date, lai, canht, lat, lon, lai_std, '
                                                                'canht_std')

    save_arr_sm = np.array([date_str_sm, sm_mean, np.array([lat] * len(sm_dates)), np.array([lon] * len(sm_dates)),
                            sm_std]).T
    fname_sm = fdir + '/mni_sm_field_' + str(point_no) + '_' + str(point_level) + '.csv'
    np.savetxt(fname_sm, save_arr_sm, fmt='%s', delimiter=',', header='date, sm, lat, lon, sm_std')

    return sm_dates, lat, lon


def get_dates_lat_lon(state_sm_dir, point_no, point_level='high'):
    sm_loc_dat = mlab.csv2rec(state_sm_dir+'/locations_utm_epsg-32632.csv')
    easting = sm_loc_dat['point_x'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == point_level)]
    northing = sm_loc_dat['point_y'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == point_level)]
    file_head_1 = sm_loc_dat['esu_sm'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'med')][0]
    file_head_2 = sm_loc_dat['esu_sm'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'high')][0]
    file_head_3 = sm_loc_dat['esu_sm'][(sm_loc_dat['id'] == point_no) & (sm_loc_dat['esu'] == 'low')][0]
    lat, lon = utm.to_latlon(easting, northing, 32, 'U')
    sm_csv1 = glob.glob(state_sm_dir + '/' + file_head_1 +'*SM.csv')[0]
    sm_dat1 = mlab.csv2rec(sm_csv1)
    sm_csv2 = glob.glob(state_sm_dir + '/' + file_head_2 + '*SM.csv')[0]
    sm_dat2 = mlab.csv2rec(sm_csv2)
    sm_csv3 = glob.glob(state_sm_dir + '/' + file_head_3 + '*SM.csv')[0]
    sm_dat3 = mlab.csv2rec(sm_csv3)

    #sm_idx = [find_nearest(sm_dat['date'], dt.datetime.combine(x,dt.datetime.min.time()))[1] for x in dates]
    sm_dates = sm_dat1['date'][3:]
    return sm_dates, lat, lon


def extract_state_sat_csv(field_no, point_level, save_dir):
    #dates, lat, lon = extract_state_arr('/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/' + str(field_no) +
    #                                    '.csv', '/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/'
    #                                            'soil_m/', field_no, save_dir, point_level)
    #print 'extracted state variables!'
    dates, lat, lon = get_dates_lat_lon('/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/soil_m/',
                                        field_no, point_level)
    extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/munich/S1_LMU_site_2017_new.nc', lat, lon,
                    save_dir+'mni_s1_' + str(field_no) + '_' + str(point_level) + '_hv.csv', (dates[0], dates[-1]),
                    remove_pol=1)
    extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/munich/S1_LMU_site_2017_new.nc', lat, lon,
                    save_dir+'mni_s1_' + str(field_no) + '_' + str(point_level) + '_hv_vv.csv', (dates[0], dates[-1]))
    print 'extracted S1 data!'
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + str(field_no) + '_' + str(point_level) + '_b4b8.csv', (dates[0], dates[-1]))
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + str(field_no) + '_' + str(point_level) + '_allbands.csv', (dates[0], dates[-1]),
                    remove_bnds=(0,9,10))
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + str(field_no) + '_' + str(point_level) + '_b4b5b6b7b8b8a.csv',
                    (dates[0], dates[-1]), remove_bnds=(0,9,10,11,12))
    print 'extracted S2 data!'
    return 'all done! :)'


def extract_state_sat_508_latlon_csv(lat, lon, save_dir, sname):
    #dates, lat, lon = extract_state_arr('/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/' + str(field_no) +
    #                                    '.csv', '/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/'
    #                                            'soil_m/', field_no, save_dir, point_level)
    #print 'extracted state variables!'
    dates, med_lat, med_lon = get_dates_lat_lon('/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/soil_m/',
                                        508, 'med')

    extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/munich/S1_LMU_site_2017_new.nc', lat, lon,
                    save_dir+'mni_s1_' + sname + '_hv.csv', (dates[0], dates[-1]),
                    remove_pol=1)
    extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/munich/S1_LMU_site_2017_new.nc', lat, lon,
                    save_dir+'mni_s1_' + sname + '_hv_vv.csv', (dates[0], dates[-1]))
    print 'extracted S1 data!'
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + sname + '_b4b8.csv', (dates[0], dates[-1]))
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + sname + '_allbands.csv', (dates[0], dates[-1]),
                    remove_bnds=(0,9,10))
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + sname + '_b4b5b6b7b8b8a.csv',
                    (dates[0], dates[-1]), remove_bnds=(0,9,10,11,12))
    print 'extracted S2 data!'
    return 'all done! :)'


def extract_state_sat_ita_latlon_csv(lat, lon, save_dir, sname):
    #dates, lat, lon = extract_state_arr('/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/' + str(field_no) +
    #                                    '.csv', '/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/'
    #                                            'soil_m/', field_no, save_dir, point_level)
    #print 'extracted state variables!'
    dates = [dt.datetime(2017, 7,1,0,0), dt.datetime(2017,11,3,0,0)]

    extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/italy/foggia_italy_processed.nc', lat, lon,
                    save_dir+'ita_s1_' + sname + '_hv.csv', (dates[0], dates[-1]),
                    remove_pol=1)
    #extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/italy/foggia_italy_processed.nc', lat, lon,
    #                save_dir+'ita_s1_' + sname + '_hv_vv.csv', (dates[0], dates[-1]))
    print 'extracted S1 data!'
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/italy/ITA_181001.nc', lat, lon,
                    save_dir+'ita_s2_' + sname + '_b4b8.csv', (dates[0], dates[-1]))
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/italy/ITA_181001.nc', lat, lon,
                    save_dir+'ita_s2_' + sname + '_allbands.csv', (dates[0], dates[-1]),
                    remove_bnds=(0,9,10))
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/italy/ITA_181001.nc', lat, lon,
                    save_dir+'ita_s2_' + sname + '_b4b5b6b7b8b8a.csv',
                    (dates[0], dates[-1]), remove_bnds=(0,9,10,11,12))
    print 'extracted S2 data!'
    return 'all done! :)'

def extract_state_sat_csv_2017(field_no, point_level, save_dir):
    dates, lat, lon = extract_state_arr('/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/' + str(field_no) +
                                        '.csv', '/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/'
                                                'soil_m/', field_no, save_dir, point_level)
    print 'extracted state variables!'
    extract_s1_data('/export/cloud/nceo/users/if910917/sentinel_data/S1/munich/S1_LMU_site_2017_new.nc', lat, lon,
                    save_dir+'mni_s1_' + str(field_no) + '_' + str(point_level) + '_2017.csv')
    print 'extracted S1 data!'
    extract_s2_data('/export/cloud/nceo/users/if910917/sentinel_data/S2/munich/MNI_181001.nc', lat, lon,
                    save_dir+'mni_s2_' + str(field_no) + '_' + str(point_level) + '_2017.csv')
    print 'extracted S2 data!'
    return 'all done! :)'