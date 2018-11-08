import numpy as np
import signaturesimulator as ss
import scipy.optimize as spop
import f90nml
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def extract_s1_mod_obs_508(point='high', sitenml='site.nml'):
    field_dir = '/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/field_508/'
    state_name = field_dir+'mni_state_field_508_'+point+'.csv'
    s1_name = field_dir+'mni_s1_508_'+point+'.csv'
    s1_arr = mlab.csv2rec(s1_name, comments='%')
    sim = ss.Simulator(site_nml=sitenml)
    sim.get_land_state=sim.state_csv(state_name)
    sim.get_geom=sim.geom_csv(s1_name)
    sim.run_rt = sim.active_microwave
    sim.run()
    return sim.backscat.date_sat_ob, sim.backscat.hv, sim.backscat.vv, s1_arr['hv'], s1_arr['vv']


def extract_s2_mod_obs_508(point='high', sitenml='site.nml', bnd_indx=(3,7)):
    field_dir = '/export/cloud/nceo/users/if910917/sentinel_data/field_data/munich/field_508/'
    state_name = field_dir+'mni_state_field_508_'+point+'.csv'
    s2_name = field_dir+'mni_s2_508_'+point+'.csv'
    s2_arr = mlab.csv2rec(s2_name, comments='%')
    bnds = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8a', 'b8b', 'b9', 'b10', 'b11', 'b12']
    sim = ss.Simulator(site_nml=sitenml)
    sim.get_land_state=sim.state_csv(state_name)
    sim.get_geom=sim.geom_csv(s2_name)
    sim.run_rt = sim.passive_optical
    sim.run()
    refl = np.array(sim.spectra.refl)
    return sim.spectra.date_sat_ob, refl[:,bnd_indx[0]], refl[:,bnd_indx[1]], s2_arr[bnds[bnd_indx[0]]], \
           s2_arr[bnds[bnd_indx[1]]]


def extract_all_points_s1(sitenml='site.nml'):
    date, hi_hv_mod, hi_vv_mod, hi_hv_ob, hi_vv_ob = extract_s1_mod_obs_508('high', sitenml)
    date, lo_hv_mod, lo_vv_mod, lo_hv_ob, lo_vv_ob = extract_s1_mod_obs_508('low', sitenml)
    date, med_hv_mod, med_vv_mod, med_hv_ob, med_vv_ob = extract_s1_mod_obs_508('med', sitenml)
    date, mean_hv_mod, mean_vv_mod, mean_hv_ob, mean_vv_ob = extract_s1_mod_obs_508('mean', sitenml)
    hv_mod = np.array([hi_hv_mod,lo_hv_mod,med_hv_mod,mean_hv_mod])#.flatten()
    vv_mod = np.array([hi_vv_mod,lo_vv_mod,med_vv_mod,mean_vv_mod])#.flatten()
    hv_ob = np.array([hi_hv_ob,lo_hv_ob,med_hv_ob,mean_hv_ob])#.flatten()
    vv_ob = np.array([hi_vv_ob,lo_vv_ob,med_vv_ob,mean_vv_ob])#.flatten()
    return date, hv_mod, vv_mod, hv_ob, vv_ob


def extract_all_points_s2(sitenml='site.nml'):
    date, hi_b4_mod, hi_b8a_mod, hi_b4_ob, hi_b8a_ob = extract_s2_mod_obs_508('high', sitenml)
    date, lo_b4_mod, lo_b8a_mod, lo_b4_ob, lo_b8a_ob = extract_s2_mod_obs_508('low', sitenml)
    date, med_b4_mod, med_b8a_mod, med_b4_ob, med_b8a_ob = extract_s2_mod_obs_508('med', sitenml)
    date, mean_b4_mod, mean_b8a_mod, mean_b4_ob, mean_b8a_ob = extract_s2_mod_obs_508('mean', sitenml)
    b4_mod = np.array([hi_b4_mod,lo_b4_mod,med_b4_mod,mean_b4_mod])
    b8a_mod = np.array([hi_b8a_mod,lo_b8a_mod,med_b8a_mod,mean_b8a_mod])
    b4_ob = np.array([hi_b4_ob,lo_b4_ob,med_b4_ob,mean_b4_ob])
    b8a_ob = np.array([hi_b8a_ob,lo_b8a_ob,med_b8a_ob,mean_b8a_ob])
    return date, b4_mod, b8a_mod, b4_ob, b8a_ob


def s1_rmse(x0, sitenml='site.nml', point='high'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['s'] = x0[0]
    nml_dic['site_params']['lai_coeff'] = x0[1]
    nml_dic['site_params']['omega'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, hv_mod, vv_mod, hv_ob, vv_ob = extract_s1_mod_obs_508(point, sitenml)
    rmse_hv = np.sqrt(np.sum((hv_mod - hv_ob)**2) / len(hv_ob.flatten()))
    rmse_vv = np.sqrt(np.sum((vv_mod - vv_ob)**2) / len(vv_ob.flatten()))
    return rmse_vv #+ rmse_hv


def s1_r2(x0, sitenml='site.nml'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['s'] = x0[0]
    nml_dic['site_params']['lai_coeff'] = x0[1]
    nml_dic['site_params']['omega'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, hv_mod, vv_mod, hv_ob, vv_ob = extract_all_points_s1(sitenml)
    ss_tot_hv = np.sum((hv_ob - np.mean(hv_ob))**2)
    ss_res_hv = np.sum((hv_mod - np.mean(hv_ob))**2)
    hv_r2 = 1 - ss_res_hv / ss_tot_hv
    ss_tot_vv = np.sum((vv_ob - np.mean(vv_ob)) ** 2)
    ss_res_vv = np.sum((vv_mod - np.mean(vv_ob)) ** 2)
    vv_r2 = 1 - ss_res_vv / ss_tot_vv
    return hv_r2, vv_r2


def s1_plot_hv(x0, sitenml='site.nml', point='high'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['s'] = x0[0]
    nml_dic['site_params']['lai_coeff'] = x0[1]
    nml_dic['site_params']['omega'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, hv_mod, vv_mod, hv_ob, vv_ob = extract_s1_mod_obs_508(point, sitenml)
    fig, ax = plt.subplots()
    ax.plot(date, hv_mod, 'o', label='hv simulated')
    ax.plot(date, hv_ob, 'X', label='hv observation')
    ax.set_xlabel('Date')
    ax.set_ylabel(r'Backscatter( m$^{2}$ m$^{-2}$)')
    fig.autofmt_xdate()
    plt.legend()
    plt.show()
    return 'done'


def s1_plot_vv(x0, sitenml='site.nml', point='high'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['s'] = x0[0]
    nml_dic['site_params']['lai_coeff'] = x0[1]
    nml_dic['site_params']['omega'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, hv_mod, vv_mod, hv_ob, vv_ob = extract_s1_mod_obs_508(point, sitenml)
    fig, ax = plt.subplots()
    ax.plot(date, vv_mod, 'o', label='vv simulated')
    ax.plot(date, vv_ob, 'X', label='vv observation')
    ax.set_xlabel('Date')
    ax.set_ylabel(r'Backscatter( m$^{2}$ m$^{-2}$)')
    ax.legend()
    fig.autofmt_xdate()
    plt.show()
    return 'done'


def s2_rmse(x0, sitenml='site.nml', point='high'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['sm_coeff'] = x0[0]
    nml_dic['site_params']['cw'] = x0[1]
    nml_dic['site_params']['rsl1'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, b4_mod, b8a_mod, b4_ob, b8a_ob = extract_s2_mod_obs_508(point, sitenml)
    rmse_b4 = np.sqrt(np.sum((b4_mod - b4_ob)**2) / len(b4_ob.flatten()))
    rmse_b8a = np.sqrt(np.sum((b8a_mod - b8a_ob)**2) / len(b8a_ob.flatten()))
    return rmse_b4 + rmse_b8a


def s2_r2(x0, sitenml='site.nml', point='high'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['sm_coeff'] = x0[0]
    nml_dic['site_params']['cw'] = x0[1]
    nml_dic['site_params']['rsl1'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, b4_mod, b8a_mod, b4_ob, b8a_ob = extract_s2_mod_obs_508(point, sitenml)
    ss_tot_b4 = np.sum((b4_ob - np.mean(b4_ob))**2)
    ss_res_b4 = np.sum((b4_mod - np.mean(b4_ob))**2)
    b4_r2 = 1 - ss_res_b4 / ss_tot_b4
    ss_tot_b8a = np.sum((b8a_ob - np.mean(b8a_ob)) ** 2)
    ss_res_b8a = np.sum((b8a_mod - np.mean(b8a_ob)) ** 2)
    b8a_r2 = 1 - ss_res_b8a / ss_tot_b8a
    return b4_r2, b8a_r2


def s2_plot(x0, sitenml='site.nml', point='high'):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['sm_coeff'] = x0[0]
    nml_dic['site_params']['cw'] = x0[1]
    nml_dic['site_params']['rsl1'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, b4_mod, b8a_mod, b4_ob, b8a_ob = extract_s2_mod_obs_508(point, sitenml)
    fig, ax = plt.subplots()
    ax.plot(date, b4_mod, 'o', label='band 4 simulated')
    ax.plot(date, b4_ob, 'X', label='band 4 observed')
    ax.plot(date, b8a_mod, 'o', label='band 8a simulation')
    ax.plot(date, b8a_ob, 'X', label='band 8a observed')
    ax.set_xlabel('Date')
    ax.set_ylabel('Band Reflectance')
    fig.autofmt_xdate()
    plt.legend()
    plt.show()
    return 'done'


def s2_plot_othbnds(x0, sitenml='site.nml', point='high', bnd_idx1=3, bnd_idx2=7):
    nml_dic = f90nml.read(sitenml)
    nml_dic['site_params']['sm_coeff'] = x0[0]
    nml_dic['site_params']['cw'] = x0[1]
    nml_dic['site_params']['rsl1'] = x0[2]
    nml_dic.write(sitenml, force=True)
    date, b4_mod, b8a_mod, b4_ob, b8a_ob = extract_s2_mod_obs_508(point, sitenml, bnd_indx=(bnd_idx1,bnd_idx2))
    fig, ax = plt.subplots()
    bnds = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8a', 'b8b', 'b9', 'b10', 'b11', 'b12']
    ax.plot(date, b4_mod, 'o', label=bnds[bnd_idx1]+' simulated')
    ax.plot(date, b4_ob, 'X', label=bnds[bnd_idx1]+' observed')
    ax.plot(date, b8a_mod, 'o', label=bnds[bnd_idx2]+' simulation')
    ax.plot(date, b8a_ob, 'X', label=bnds[bnd_idx2]+' observed')
    ax.set_xlabel('Date')
    ax.set_ylabel('Band Reflectance')
    fig.autofmt_xdate()
    plt.legend()
    plt.show()
    return 'done'


def optimize_s1(bnds=[(0.01, 1.5), (0.001, 0.8), (0.08, 0.2)], sitenml='site.nml', point='high'):
    s = 0.015
    lai_coeff = 0.1
    omega = 0.12
    x0 = np.array([s, lai_coeff, omega])
    #res = spop.fmin_l_bfgs_b(s1_rmse, x0, bounds=bnds, approx_grad=True)
    res = spop.differential_evolution(s1_rmse, bnds, args=(sitenml, point))
    return res


def optimize_s2(bnds=[(0.0, 1.0), (0.001, 0.1), (0.01, 0.8)], sitenml='site.nml', point='high'):
    sm_coef = 0.5
    cw = 0.01
    rsl1 = 0.2
    x0 = np.array([sm_coef, cw, rsl1])
    #res = spop.fmin_l_bfgs_b(s2_rmse, x0, bounds=bnds, approx_grad=True)
    #res = spop.minimize(s2_rmse, x0, method='Nelder-Mead')
    res = spop.differential_evolution(s2_rmse, bnds, args=(sitenml, point))
    return res

