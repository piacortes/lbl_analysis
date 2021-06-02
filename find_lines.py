from astropy.table import Table
from astropy.io import fits as pf
import matplotlib.pyplot as plt
import numpy as np
import glob
import time
from astropy.timeseries import LombScargle
import astropy.units as u

files = glob.glob('/net/GSP/users/pcortes/LBL/lblrv_GL205/*.fits')
target = 'GL205'
print('Processing {0} files'.format(len(files)))

fits = pf.open(files[0])
wave_list = np.zeros(len(fits[1].data['RV']))
for i in range(len(fits[1].data['RV'])):
    wave = fits[1].data['WAVE_START'][i]
    wave_list[i] = wave

# JD array, the same for all time series
mjd = np.zeros(len(files))
for i in range(len(files)):
    fits = pf.open(files[i])
    mjd[i] = fits[0].header['MJDMID']

# Remove spectrum zones with high atmospheric absorption
waves = np.unique(wave_list[((wave_list<1800) | (wave_list>2000)) & ((wave_list<1300) | (wave_list>1500))])
print('Analysing {0} lines'.format(len(waves)))

print('Building big table of lines...')
start = time.time()
#files = files[:50]
tbl = np.zeros([len(files),len(waves),8])
for f in range(len(files)):
    fits = pf.open(files[f])
    #mjd = fits[0].header['MJDMID']
    #tbl[:,:,f] = mjd
    for i in range(len(waves)):
        tbl[f,i,1] = waves[i]
        #fits_tbl = fits[1].data[['ORDER','RV','DVRMS','DDV','DDVRMS','DDDV','DDDVRMS'][fits[1].data['WAVE_START'] == waves[i]]
        mask_wave = tuple([fits[1].data['WAVE_START'] == waves[i]])
        rvs = fits[1].data['RV'][mask_wave]
        chis = fits[1].data['DVRMS'][mask_wave]
        orders = fits[1].data['ORDER'][mask_wave]
        dvrmss = fits[1].data['DVRMS'][mask_wave]
        ddvs = fits[1].data['DDV'][mask_wave]
        ddvrmss = fits[1].data['DDVRMS'][mask_wave]
        dddvs = fits[1].data['DDDV'][mask_wave]
        dddvrmss =  fits[1].data['DDDVRMS'][mask_wave]

        rvs_nan = np.isnan(rvs)
        if rvs_nan.all() == True:
            best_rv = order = dvrms = ddv = ddvrms = dddv = np.nan
        else:
            min_indx = np.nanargmin(chis)
            best_rv = rvs[min_indx]
            order = orders[min_indx]
            dvrms = dvrmss[min_indx]
            ddv = ddvs[min_indx]
            ddvrms = ddvrmss[min_indx]
            dddv = ddvs[min_indx]
            dddvrms = ddvrmss[min_indx]
            #print(best_rv)
        tbl[f,i,0] = order
        tbl[f,i,2] = best_rv
        tbl[f,i,3] = dvrms
        tbl[f,i,4] = ddv
        tbl[f,i,5] = ddvrms
        tbl[f,i,6] = dddv
        tbl[f,i,7] = dddvrms

time.sleep(1)
end = time.time()
print('Table of lines ready. It took {0} seconds.'.format(np.round((end-start),2)))

tbl_rv = np.zeros(len(waves),8)
tbl_ddv = np.zeros(len(waves),8)
tbl_dddv = np.zeros(len(waves),8)

print("Let's compute the periodograms...")
t = mjd * u.m/u.s
for i in range(len(waves)):
    rv_line = tbl[:,i,2] * u.m/u.s
    rv_line_err = tbl[:,i,3] * u.m/u.s
    ddv_line = tbl[:,i,4] #* u.m/u.s
    ddv_line_err = tbl[:,i,5] #* u.m/u.s
    dddv_line = tbl[:,i,6] #* u.m/u.s
    dddv_line_err = tbl[:,i,7] #* u.m/u.s

    ls = LombScargle(t,rv_line, rv_line_err)
    frequency, power = ls.autopower()
    ind_peaks = np.argpartition(power, -5)[-5:] #highest five peaks in order
    power_peaks = power[ind_peaks]
    frequency_peaks = frequency[ind_peaks]
    tbl_rv[i,0] = waves[i]
    tbl_rv[i,1] = np.nanmean(rv_line)
    tbl_rv[i,2] = np.nanstd(rv_line)
    tbl_rv[i,3:] = frequency_peaks

    ls = LombScargle(t,ddv_line, ddv_line_err)
    frequency, power = ls.autopower()
    ind_peaks = np.argpartition(power, -5)[-5:] #highest five peaks in order
    power_peaks = power[ind_peaks]
    frequency_peaks = frequency[ind_peaks]
    tbl_ddv[i,0] = waves[i]
    tbl_ddv[i,1] = np.nanmean(ddv_line)
    tbl_ddv[i,2] = np.nanstd(ddv_line)
    tbl_ddv[i,3:] = frequency_peaks

    ls = LombScargle(t,dddv_line, dddv_line_err)
    frequency, power = ls.autopower()
    ind_peaks = np.argpartition(power, -5)[-5:] #highest five peaks in order
    power_peaks = power[ind_peaks]
    frequency_peaks = frequency[ind_peaks]
    tbl_dddv[i,0] = waves[i]
    tbl_dddv[i,1] = np.nanmean(dddv_line)
    tbl_dddv[i,2] = np.nanstd(dddv_line)
    tbl_dddv[i,3:] = frequency_peaks

np.savetxt(target+'_lines_RV.csv',tbl_rv, header='mjd meanRV stdRV peak1 peak2 peak3 peak4 peak5')
np.savetxt(target+'_lines_DDV.csv',tbl_ddv, header='mjd meanDDV stdDDV peak1 peak2 peak3 peak4 peak5')
np.savetxt(target+'_lines_DDDV.csv',tbl_dddv, header='mjd meanDDDV stdDDDV peak1 peak2 peak3 peak4 peak5')

print("DONE.")
