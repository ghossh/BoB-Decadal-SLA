import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

## Data variable rename 
def name_ch(input_data):
    #lon=input_data.dims[2]
    #lat=input_data.dims[1]
    #time=input_data.dims[0]
    input_data=input_data.rename({'lat':'latitude'})
    input_data=input_data.rename({'lon':'longitude'})
    input_data=input_data.rename({'time':'time'})
    return input_data

def time_set(data):
    tim=pd.date_range("1-2-1991", periods=len(data.time))
    tim_da = xr.DataArray(tim,[('time',tim)])
    data['time']=tim_da
    return data
def sla_cal(data):
    data['sl']=data.pn/980
    clim=data.sl.sel(time=slice('2000-01-01','2019-12-31')).mean('time')
    data['sla']=data.sl-clim
    data=data.sla
    data=data.where(data<1000)
    data=data.sel(time=slice('1993-01-01','2019-12-31'))
    #data=data.where(data>0.0001)
    return data
def data_proc(data):
    #data=data.pn
    data=name_ch(data)
    data=time_set(data)
    data=sla_cal(data)
    return data


import scipy.signal as signal
def butterworth_lowpass_filter(data, cutoff_time, axis=0):
    cutoff_freq=1/cutoff_time
    nyfreq=cutoff_freq*2
    order=4
    B, A = signal.butter(order, nyfreq, output="ba")
    return signal.filtfilt(B, A, data, axis=0)

def filter_data(data):
    data_low = butterworth_lowpass_filter(data, 400)
    data_low=xr.DataArray(data_low,coords={"time": data.time, "latitude": data.latitude,
                              "longitude": data.longitude},
                      dims=["time", "latitude","longitude"])
    return data_low

BoB_altimeter_raw=xr.open_mfdataset("/home/NCAOR/supriyog/raw_data/SL_no_anom/global_SL_and_current_93_19.nc"
                                 ).sel(latitude=slice(4,25),longitude=slice(78,100)).sla

clim=BoB_altimeter_raw.sel(time=slice('2000-01-01','2019-12-30')).mean('time')


BoB_altimeter=(BoB_altimeter_raw-clim)*100

lcs_path='/home/NCAOR/supriyog/raw_data/LCSCR_model_output/LCSCR_run_by_supi/LCS_ERA5/'
BoB_LCSCR=xr.open_dataset(lcs_path+'LCSCR_final_era5.nc',autoclose=True
                         ).sel(lat=slice(4,25),lon=slice(78,100))
#BoB_LCSCR=BoB_LCSCR.rename({'lon':'longitude'})
#BoB_LCSCR=BoB_LCSCR.rename({'lat':'latitude'})
BoB_LCSCR.dims

BoB_LCSCR=xr.open_mfdataset(lcs_path+'LCSCR_final_era5.nc' ).sel(lat=slice(4,25),lon=slice(78,100))
BoB_LCSBI=xr.open_mfdataset(lcs_path+'LCSBI.nc').sel(lat=slice(4,25),lon=slice(78,100))
BoB_LCSEB=xr.open_mfdataset(lcs_path+'LCSEB.nc').sel(lat=slice(4,25),lon=slice(78,100))
BoB_LCSWB=xr.open_mfdataset(lcs_path+'LCSWB.nc').sel(lat=slice(4,25),lon=slice(78,100))
BoB_LCSEIO=xr.open_mfdataset(lcs_path+'LCSEIO.nc').sel(lat=slice(4,25),lon=slice(78,100))


BoB_LCSCR=data_proc(BoB_LCSCR)
BoB_LCSBI=data_proc(BoB_LCSBI)
BoB_LCSEB=data_proc(BoB_LCSEB)
BoB_LCSWB=data_proc(BoB_LCSWB)
BoB_LCSEIO=data_proc(BoB_LCSEIO)


# BoB_altimeter.to_netcdf('./unfiltered/BoB_altimeter.nc')
# BoB_LCSCR.to_netcdf('./unfiltered/BoB_LCSCR.nc')
# BoB_LCSBI.to_netcdf('./unfiltered/BoB_LCSBI.nc')
# BoB_LCSEB.to_netcdf('./unfiltered/BoB_LCSEB.nc')
# BoB_LCSWB.to_netcdf('./unfiltered/BoB_LCSWB.nc')
# BoB_LCSEIO.to_netcdf('./unfiltered/BoB_LCSEIO.nc')

# BoB_altimeter_low = filter_data(BoB_altimeter)
# BoB_LCSCR_low = filter_data(BoB_LCSCR)
BoB_LCSWB_low = filter_data(BoB_LCSWB)
BoB_LCSEB_low = filter_data(BoB_LCSEB)
BoB_LCSBI_low = filter_data(BoB_LCSBI)
BoB_LCSEIO_low = filter_data(BoB_LCSEIO)

# BoB_altimeter_low.to_netcdf('./filtered/BoB_altimeter_low.nc')
# BoB_LCSCR_low.to_netcdf('./filtered/BoB_LCSCR_low.nc')
BoB_LCSBI_low.to_netcdf('./filtered/BoB_LCSBI_low.nc')
BoB_LCSEB_low.to_netcdf('./filtered/BoB_LCSEB_low.nc')
BoB_LCSWB_low.to_netcdf('./filtered/BoB_LCSWB_low.nc')
BoB_LCSEIO_low.to_netcdf('./filtered/BoB_LCSEIO_low.nc')

