# -*- coding: utf-8 -*-
"""


@author: TomClifford

This file is an initial routine for preprocessing seismic data.
It reads a waveform, filters, removes response, demeans and detrends, finds SNR, FAS, and plots.

"""

#%% import libraries

import obspy
from obspy.clients.fdsn.mass_downloader import CircularDomain, \
    Restrictions, MassDownloader
from obspy.io.xseed import Parser
from obspy.signal import PPSD
from obspy.signal import freqattributes
import os
from scipy.fft import fft, ifft, fftfreq
from scipy.integrate import cumtrapz
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from response_spectrum import *


#%%paths
data_path = r"C:\Users\TomClifford\SlateGeotech\Duynefontyn PSHA - DuynefontynPSHA\05 - GMM\GMM_Scripts\preprocessing"
os.chdir(data_path)


#%% read waveform data into an obspy stream object, st

#origin time
origin_time = obspy.UTCDateTime(2016,10,18,6,25,33)

#https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html
st = obspy.read(os.path.join(data_path, 'waveforms/*'))
print(st)




#%%remove response
inv = obspy.read_inventory(os.path.join(data_path, 'stations/*'))
# response_list = os.listdir(os.path.join(data_path, 'stations'))

#create empty stream object to add once response removed
st_r = obspy.core.stream.Stream()

#prefilter, remove response, and append new trace to stream
#tr is a waveform trace: https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.html
for tr in st:
    #determine nyquist f
    nyq = tr.stats.sampling_rate/2
    #set prefilter according to nyquist
    prefilter = [0.001, 0.005, nyq-5, nyq]
    #find matching response
    tr_response = tr.copy()
    tr_response.remove_response(inventory=inv,
                                pre_filt = prefilter,
                                output = "ACC", 
                                )
    st_r.append(tr_response)
    # print(tr_response)
    

st_rd = st_r.copy()
#https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.detrend.html
#demean
st_rd.detrend('demean')
#detrend
st_rd.detrend('linear')
    

#trim waveform
st_rd.trim(origin_time, origin_time+(1000/2))

for tr in st_rd:
    print(tr)
    tr.plot()
#%%SNR
def snr(trace):
    #input: obspy trace object
    peak = trace.data.max()
    rms =   np.sqrt(np.mean(trace.data**2))
    snr = peak/rms
    return snr
    
    
#%% test snr


for tr in st_rd:
    # tr.plot()
    print(snr(tr))
    
    
#%% FAS


def fas(tr):
    #tr: obspy trace object
    y = fft(tr.data)
    yf = 2.0/tr.stats.npts*np.abs(y[:(tr.stats.npts//2)])
    xf = fftfreq(tr.stats.npts, tr.stats.delta)[:(tr.stats.npts//2)]
    return xf, yf
    


x, y = fas(st_rd[0])


#plot
fig, ax = plt.subplots()
ax.plot(x, y, lw=0.3, color='k')

ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Amplitude")
plt.yscale('log')
plt.show()


#%% response spectra

#get response spectra from trace
#seems number of periods has to be same length as trace?
r = NewmarkBeta([tr.times(), tr.data/100], tr.stats.delta, np.logspace(.1, 10, len(tr.data)))#np.array([.1,1,10])) #convert to cm/s/s

#why returning period instead of time
plt.plot(r.response_spectrum['Period'], r.response_spectrum['Acceleration']) #so this is the waveform


#%% save trace amplitudes, times, and fourier spectra to excel

for tr in st:
    print(tr)

    trace_data = pd.DataFrame({'trace_amplitudes': tr.data,
                               'trace_time' : tr.times()
                               })
    trace_fft = pd.DataFrame({'fftx': fas(tr)[0],
                               'ffty': fas(tr)[1]
                               })
    trace_data.to_csv(os.path.join(data_path, 'raw_traces', str(tr.id)+'_data.csv'))
    trace_fft.to_csv(os.path.join(data_path, 'raw_traces', str(tr.id)+'_fft.csv'))



#%% download data for event - not for final script

#M4.5  2016-10-18 06:25:33.160000
# origin_time = obspy.UTCDateTime(2016,10,18,6,25,33)

# domain = CircularDomain(latitude=-33.676667, longitude=18.431389,
#                                 minradius=0.0, maxradius= 1000/111) #1000 km to deg
# restrictions = Restrictions(
#     #5 minutes before, 30 minutes after origin
#     starttime=origin_time - 5 * 60,
#     endtime=origin_time + 30*60,
#     network='*',
#     station = '*',
#     location='*',
#     # channel='*', #allowing all channels downloads non-seismic data
#     reject_channels_with_gaps=False,
#     )

# mdl = MassDownloader(providers=['IRIS'])

# mdl.download(domain, restrictions, mseed_storage=data_path+"/waveforms",
#               stationxml_storage=data_path+"/stations")

# #downloads 13 waveforms from 3 stations














