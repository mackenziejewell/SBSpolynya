

# DEPENDENCIES:
import pandas as pd
from datetime import datetime
import xarray as xr
import sys
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np, numpy.ma as ma
import matplotlib.cm as cm
import cartopy, cartopy.crs as ccrs

import metpy.calc as mpcalc
from metpy.units import units

# cartopy plotting
sys.path.append('../')
import plot_simply.geomap as geomap
import geofunc.vectors as vectors
import data_NSIDC.icedrift as icedrift

# standard maps to re-use throughout codes
def makemap(view = 'wide', contours = [], figsize=(8,6), panels=(1,1)):

    map_proj = ccrs.NorthPolarStereo(central_longitude=-140)
#     if view in ['shelf_view', 'small_region']:
#         map_proj = ccrs.NorthPolarStereo(central_longitude=-143)
        
    if panels==(1,1):
        fig, ax = plt.subplots(subplot_kw=dict(projection=map_proj), figsize=figsize)
        axs = [ax]

    else:
        fig, AXS = plt.subplots(*panels, subplot_kw=dict(projection=map_proj), figsize=figsize)
        axs = AXS.ravel()

    for ax in axs:
        ax.set_facecolor('lightgray')

        geomap.land(ax, scale = '10m', color='darkgray', alpha=1, fill_dateline_gap = False, zorder=2)

        if len(contours)>0:
            geomap.gebco_bathymetry(ax, file_path='/Volumes/Seagate_Jewell/KenzieStuff/GEBCO/GEBCO_2024/gebco_2024_n90.0_s55.0_w-180.0_e180.0.nc', 
                                crop_lat=(69, 73.5), crop_lon=(-170, -110), clat=8, clon=15, depth_shade=False, 
                                shade_zorder=0, depth_contours=True, contour_levels=contours, 
                                contour_kwargs={'colors': 'gray', 'linewidths': 0.75, 'linestyles': 'solid', 'zorder': 100},
                                contour_labels=False, text_kwargs={'size': 10, 'color': 'gray', 'weight': 'normal', 'zorder': 100})


        if view == 'wide':
            ax.set_xlim(-700000,500000)
            ax.set_ylim(-2400000,-1900000)
            
        elif view == 'very_tall':
            ax.set_ylim(-2400000,-1300000)
            ax.set_xlim(-600000,450000)


        elif view == 'very_wide':
            ax.set_ylim(-2400000,-1100000)
            ax.set_xlim(-1000000,500000)
            
        elif view == 'circle_view':
            import matplotlib.path as mpath
            ax.set_extent([0, 360,  60, 90], ccrs.PlateCarree())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)

            ax.set_boundary(circle, transform=ax.transAxes)

        elif view == 'large_view':
            ax.set_extent([199, 235,  67.5, 80], ccrs.PlateCarree())



#             ax.set_extent([-180, -105, 68, 90], ccrs.PlateCarree())
            # ax.set_ylim(-2500000,-1000000)
            # ax.set_xlim(-1200000,800000)
        
        elif view == 'shelf_view':
            ax.set_ylim(-2400000,-1850000)
            ax.set_xlim(-780000,350000)
            
#         elif view == 'small_region':
#             ax.set_ylim(-2350000,-2020000)
#             ax.set_xlim(-200000,330000)
            
         
            
        elif view == 'zoom':
            ax.set_ylim(-2400000,-2050000)
            ax.set_xlim(-220000,180000)


        elif view == 'wider_zoom':
            ax.set_ylim(-2400000,-1950000)
            ax.set_xlim(-500000,280000)
            
            
        elif view == 'tall_zoom':
            ax.set_ylim(-2395000,-1940000)
            ax.set_xlim(-420000,185000)

        elif view == 'wider_zoom2':
            ax.set_ylim(-2350000,-2020000)
            ax.set_xlim(-400000,130000)


        elif view == 'reallyzoom':
            ax.set_ylim(-2360000,-2150000)
            ax.set_xlim(-200000,160000)

        elif view == 'narrowzoom':
            ax.set_ylim(-2330000,-2180000)
            ax.set_xlim(-170000,100000)

    if panels==(1,1):
        axs = ax
    return fig, axs
        

def open_daily_winds(year, lat_range, lon_range, time_range = None):
    
    ds = xr.open_dataset(f'/Volumes/Jewell_EasyStore/ECMWF/annual/daily/ERA5_{year}_daily.nc')
    ds.close()
    if time_range == None:
        ds = ds.sel(latitude = lat_range, longitude = lon_range)
    else:
        ds = ds.sel(time = time_range, latitude = lat_range, longitude = lon_range)
        
    return ds

def open_daily_t2m(year, lat_range, lon_range, time_range = None):
    
    ds = xr.open_dataset(f'/Volumes/Seagate_Jewell/KenzieStuff/ERA5/daily_t2m/ERA5_T2m_daily_{year}.nc')
    ds.close()
    
    if time_range == None:
        ds = ds.sel(latitude = lat_range, longitude = lon_range)
    else:
        ds = ds.sel(valid_time = time_range, latitude = lat_range, longitude = lon_range)
    
    return ds

def open_daily_drift(year, lat_range, lon_range, time_range = None):
    
    if type(time_range) == pd.Timestamp or type(time_range) == datetime:
        dates = time_range
    elif time_range == None:
        dates = pd.to_datetime(pd.date_range(datetime(year,1,1),datetime(year,12,31), freq='1D'))
    else:
        dates = pd.to_datetime(pd.date_range(time_range.start, time_range.stop, freq='1D'))
        
    drift = icedrift.open_local_file(dates, crop = [180,310,80,220],include_units = False)
    
    
    drift2 = {}
    drift2['u'] = drift['u']
    drift2['v'] = drift['v']
    drift2['e'] = drift['e']
    drift2['n'] = drift['n']
    
    drift2['proj'] = drift['proj']
    drift2['xx'] = drift['xx']
    drift2['yy'] = drift['yy']
    drift2['lon'] = drift['lon']
    drift2['lat'] = drift['lat']
    
    drift2['s'] = np.sqrt(drift['e']**2+drift['n']**2)

    return drift2



def drift_map_over_time(dates, map_proj, crop = True):

    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
 

    if crop:
        crop = [220,265,110,160]
    else:
        crop = [0,None,0,None]
    drift_map = icedrift.open_local_file(pd.to_datetime(dates),
                        main_path = '/Volumes/Jewell_EasyStore/NSIDC-0116_PPdrift/', 
                        filenametype = 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc', 
                        crop = crop,
                        include_units = False)
    
    E = np.nanmean(drift_map['e'], axis=0)
    N = np.nanmean(drift_map['n'], axis=0)
    # Ep, Np = geomap.fix_cartopy_vectors(E,N,drift_map['lat'])

    E_scaled = np.zeros(E.shape)
    N_scaled = np.zeros(N.shape)

    X_scaled = np.zeros(E.shape)
    Y_scaled = np.zeros(E.shape)
    
    for ii in range(np.shape(E_scaled)[0]):
        for jj in range(np.shape(N_scaled)[0]):

            tail, tip, vec = vectors.project_vectors(map_proj, 
                                                    drift_map['lon'][ii,jj], 
                                                    drift_map['lat'][ii,jj], 
                                                    eastward = E[ii,jj]*units('cm/s'), 
                                                    northward = N[ii,jj]*units('cm/s'), 
                                            final_units = 'm/day')

            E_scaled[ii,jj] = vec[0].magnitude
            N_scaled[ii,jj] = vec[1].magnitude

            X_scaled[ii,jj] = tail[0]
            Y_scaled[ii,jj] = tail[1]
            
        drift_map['E_scaled'] = E_scaled
        drift_map['N_scaled'] = N_scaled
        
        drift_map['X_scaled'] = X_scaled
        drift_map['Y_scaled'] = Y_scaled

    warnings.filterwarnings("default", category=DeprecationWarning) 
    warnings.filterwarnings("default", category=RuntimeWarning)        

    return drift_map


def wind_map_over_time(dates, map_proj, era_lat = slice(75, 68), era_lon = slice(-158,-125)):


    era_map = {}

    counter = 0
    for date in pd.to_datetime(dates):
        try:
            filename = f'/Volumes/Jewell_EasyStore/ECMWF/annual/daily/ERA5_{date.year}_daily.nc'
            with xr.open_dataset(filename) as ds:
                ds_crop = ds.sel(time=date, latitude=era_lat, longitude = era_lon)
                counter+=1
                exists=True
        except:
            print(f'missing {date}')
            exists = False

        if exists:
            if counter == 1:
#                 u10_grid = ds_crop.u10.values
#                 v10_grid = ds_crop.v10.values
                msl_grid = ds_crop.msl.values
#                 vort_grid = mpcalc.vorticity(ds_crop.u10*units('m/s'), ds_crop.v10*units('m/s'), 
#                                              latitude=ds_crop.latitude, longitude=ds_crop.longitude).values

            else:
#                 u10_grid = np.reshape(np.append(u10_grid, ds_crop.u10.values), (counter, *ds_crop.u10.values.shape))
#                 v10_grid = np.reshape(np.append(v10_grid, ds_crop.v10.values), (counter, *ds_crop.u10.values.shape))
                msl_grid = np.reshape(np.append(msl_grid, ds_crop.msl.values), (counter, *ds_crop.u10.values.shape))
#                 vort = mpcalc.vorticity(ds_crop.u10*units('m/s'), ds_crop.v10*units('m/s'), 
#                                              latitude=ds_crop.latitude, longitude=ds_crop.longitude).values
#                 vort_grid = np.reshape(np.append(vort_grid, vort), (counter, *ds_crop.u10.values.shape))

#     era_map['u10'] = u10_grid
#     era_map['v10'] = v10_grid
    era_map['msl'] = msl_grid
#     era_map['vort'] = vort_grid
    
    era_map['lon'], era_map['lat'] = np.meshgrid(ds_crop.longitude.values, ds_crop.latitude.values)

    return era_map


def find_consecutive_values(arr):
    # if not arr:
    #     return []

    result = []
    start_index = 0
    for i in range(1, len(arr)):
        if arr[i] != arr[i - 1] + 1:
            if i - start_index > 1:
                result.append((start_index, i - 1))
            start_index = i
    if len(arr) - start_index > 1:
        result.append((start_index, len(arr) - 1))
    return result