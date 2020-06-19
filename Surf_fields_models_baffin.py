
"""
Created on Thu Jun  4 13:11:41 2020

@author: aristizabal
"""

def GOFS31_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,folder_fig):
    #%% User input
    
    #GOFS3.1 output model location
    url_GOFS_ts = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z'
    
    url_GOFS_uv = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z'
    
    url_GOFS_ssh = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ssh'
    #https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ssh
    
    # Bathymetry file
    #bath_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/nc_files/GEBCO_2014_2D_-100.0_0.0_-60.0_45.0.nc'
    # Bathymetry file
    bath_file = '/home/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    
    #%%
    
    from matplotlib import pyplot as plt
    import numpy as np
    import xarray as xr
    import netCDF4 
    from datetime import datetime, timedelta
    import cmocean
    
    # Do not produce figures on screen
    plt.switch_backend('agg')
    
    # Increase fontsize of labels globally
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('legend',fontsize=14)
    
    #%% Get time bounds for the previous day
    #te = datetime.today()
    #tend = datetime(te.year,te.month,te.day)
    
    ti = datetime.today() - timedelta(1)
    tini = datetime(ti.year,ti.month,ti.day,6)
    
    #%% GOGF 3.1
    
    GOFS_ts = xr.open_dataset(url_GOFS_ts,decode_times=False)
    GOFS_uv = xr.open_dataset(url_GOFS_uv,decode_times=False)
    GOFS_ssh = xr.open_dataset(url_GOFS_ssh,decode_times=False)
    
    lt_GOFS = np.asarray(GOFS_ts['lat'][:])
    ln_GOFS = np.asarray(GOFS_ts['lon'][:])
    tt = GOFS_ts['time']
    t_GOFS = netCDF4.num2date(tt[:],tt.units) 
    
    depth_GOFS = np.asarray(GOFS_ts['depth'][:])
    
    # Conversion from glider longitude and latitude to GOFS convention
    lon_limG = np.empty((len(lon_lim),))
    lon_limG[:] = np.nan
    for i in range(len(lon_lim)):
        if lon_lim[i] < 0: 
            lon_limG[i] = 360 + lon_lim[i]
        else:
            lon_limG[i] = lon_lim[i]
        lat_limG = lat_lim
    
    oklon_GOFS = np.where(np.logical_and(ln_GOFS >= lon_limG[0],ln_GOFS <= lon_limG[1]))[0]
    oklat_GOFS = np.where(np.logical_and(lt_GOFS >= lat_limG[0],lt_GOFS <= lat_lim[1]))[0]
    oktime_GOFS = np.where(t_GOFS == tini)[0][0]
    
    # Conversion from GOFS convention to glider longitude and latitude
    ln_GOFSg= np.empty((len(ln_GOFS),))
    ln_GOFSg[:] = np.nan
    for i in range(len(ln_GOFS)):
        if ln_GOFS[i] > 180: 
            ln_GOFSg[i] = ln_GOFS[i] - 360 
        else:
            ln_GOFSg[i] = ln_GOFS[i]
    lt_GOFSg = lt_GOFS
    
    #lon_GOFS= ln_GOFS[oklon_GOFS]
    lat_GOFS= lt_GOFS[oklat_GOFS]
    lon_GOFSg= ln_GOFSg[oklon_GOFS]
    lat_GOFSg= lt_GOFSg[oklat_GOFS]
    time_GOFS = t_GOFS[oktime_GOFS]
       
    #%% Reading bathymetry data
    
    ncbath = xr.open_dataset(bath_file)
    bath_lat = ncbath.variables['lat'][:]
    bath_lon = ncbath.variables['lon'][:]
    bath_elev = ncbath.variables['elevation'][:]
    
    oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
    oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])
    
    bath_latsub = bath_lat[oklatbath]
    bath_lonsub = bath_lon[oklonbath]
    bath_elevs = bath_elev[oklatbath,:]
    bath_elevsub = bath_elevs[:,oklonbath]  
    
    #%% loading surface temperature and salinity
    
    sst_GOFS = np.asarray(GOFS_ts['water_temp'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])
    sss_GOFS = np.asarray(GOFS_ts['salinity'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])
    ssh_GOFS = np.asarray(GOFS_ssh['surf_el'][oktime_GOFS,oklat_GOFS,oklon_GOFS])
    su_GOFS = np.asarray(GOFS_uv['water_u'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])
    sv_GOFS = np.asarray(GOFS_uv['water_v'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])
    
    #%% Figure sst
    
    kw = dict(levels = np.linspace(24,30,16))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,sst_GOFS[:,:],cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='k')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SST on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig +'GOFS_SST'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure sss
    
    kw = dict(levels = np.linspace(35,37,11))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,sss_GOFS,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='k')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SSS on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_SSS'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure ssh
    
    kw = dict(levels = np.linspace(-0.6,0.6,25))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,ssh_GOFS,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='k')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14) 
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SSH on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_SSH'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure ssh
    
    kw = dict(levels = np.linspace(-0.6,0.6,25))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,ssh_GOFS,cmap=cmocean.cm.curl,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14) 
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SSH on '+str(time_GOFS)[0:13],fontsize=16)
    q=plt.quiver(lon_GOFSg[::5],lat_GOFSg[::5],su_GOFS[::5,::5],sv_GOFS[::5,::5] ,scale=3,scale_units='inches',\
              alpha=0.7)
    plt.quiverkey(q,np.max(lon_GOFSg)+0.7,np.max(lat_GOFSg)+0.4,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file = folder_fig + 'GOFS_SSH_UV'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure sst at 200 meters
    
    kw = dict(levels = np.linspace(10,25,31))
    okdepth = np.where(depth_GOFS >= 200)[0][0]
    temp_200_GOFS = np.asarray(GOFS_ts['water_temp'][oktime_GOFS,okdepth,oklat_GOFS,oklon_GOFS])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,temp_200_GOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS Temperature at 200 m on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_temp_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure salinity at 200 meters
    
    okdepth = np.where(depth_GOFS >= 200)[0][0]
    salt_200_GOFS = np.asarray(GOFS_ts['salinity'][oktime_GOFS,okdepth,oklat_GOFS,oklon_GOFS])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,salt_200_GOFS,cmap=cmocean.cm.haline)#,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS 3.1 Salinity at 200 m on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_salt_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp transect along storm path

    lon_forec_track_interp = np.interp(lat_GOFS,lat_forec_track,lon_forec_track)
    lat_forec_track_interp = lat_GOFS
    
    oklon = np.round(np.interp(lon_forec_track_interp+360,ln_GOFS,np.arange(len(ln_GOFS)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_interp,lt_GOFSg,np.arange(len(lt_GOFSg)))).astype(int)
    
    trans_temp_GOFS = np.empty((len(depth_GOFS),len(lon_forec_track_interp)))
    trans_temp_GOFS[:] = np.nan
    for i in np.arange(len(lon_forec_track_interp)):
        trans_temp_GOFS[:,i] = np.asarray(GOFS_ts['water_temp'][oktime_GOFS,:,oklat[i],oklon[i]])
    
    ok12 = np.where(trans_temp_GOFS <= 12)[0][0]
    max_depth = -depth_GOFS[ok12]
    
    kw = dict(levels = np.linspace(12,32,21))
    
    plt.figure()
    plt.contourf(lt_GOFSg[oklat],-depth_GOFS,trans_temp_GOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(lt_GOFSg[oklat],-depth_GOFS,trans_temp_GOFS,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.title('GOFS Temperature \n along Forecasted Storm Track',fontsize=16)
    plt.ylim([max_depth,0])
    #plt.xlim([20,30])
    
    file = folder_fig + 'GOFS_temp_along_forecasted_track_'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
def RTOFS_oper_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,folder_fig):

    #%% RTOFS files
    folder_RTOFS = '/home/coolgroup/RTOFS/forecasts/domains/hurricanes/RTOFS_6hourly_North_Atlantic/'
    
    nc_files_RTOFS = ['rtofs_glo_3dz_f006_6hrly_hvr_US_east.nc',\
                      'rtofs_glo_3dz_f012_6hrly_hvr_US_east.nc',\
                      'rtofs_glo_3dz_f018_6hrly_hvr_US_east.nc',\
                      'rtofs_glo_3dz_f024_6hrly_hvr_US_east.nc']
    
    # Bathymetry file
    bath_file = '/home/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    
    
    #%% 
    import xarray as xr
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime, timedelta
    import cmocean
    
    # Increase fontsize of labels globally
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('legend',fontsize=14)
    
    #%% Reading bathymetry data
    ncbath = xr.open_dataset(bath_file)
    bath_lat = ncbath.variables['lat'][:]
    bath_lon = ncbath.variables['lon'][:]
    bath_elev = ncbath.variables['elevation'][:]
    
    oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
    oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])
    
    bath_latsub = bath_lat[oklatbath]
    bath_lonsub = bath_lon[oklonbath]
    bath_elevs = bath_elev[oklatbath,:]
    bath_elevsub = bath_elevs[:,oklonbath]  
    
    #%% Get time bounds for the previous day
    #te = datetime.today()
    #tend = datetime(te.year,te.month,te.day)
    
    ti = datetime.today() - timedelta(1)
    tini = datetime(ti.year,ti.month,ti.day)
    
    #name = 'Cristoba;'
    #folder_fig = '/home/aristizabal/Figures/' + name + '/' + ti.strftime('%b-%d') + '/'
    
    #%% Time window
    
    year = int(tini.year)
    month = int(tini.month)
    day = int(tini.day)
    tini = datetime(year, month, day)
    #tend = tini + timedelta(days=1)
    
    #%% Read RTOFS grid and time
    print('Retrieving coordinates from RTOFS')
    
    if tini.month < 10:
        if tini.day < 10:
            fol = 'rtofs.' + str(tini.year) + '0' + str(tini.month) + '0' + str(tini.day)
        else:
            fol = 'rtofs.' + str(tini.year) + '0' + str(tini.month) + str(tini.day)
    else:
        if tini.day < 10:
            fol = 'rtofs.' + str(tini.year) + str(tini.month) + '0' + str(tini.day)
        else:
            fol = 'rtofs.' + str(tini.year) + str(tini.month) + str(tini.day)
    
    ncRTOFS = xr.open_dataset(folder_RTOFS + fol + '/' + nc_files_RTOFS[0])
    latRTOFS = np.asarray(ncRTOFS.Latitude[:])
    lonRTOFS = np.asarray(ncRTOFS.Longitude[:])
    depth_RTOFS = np.asarray(ncRTOFS.Depth[:])
    
    tRTOFS = []
    for t in np.arange(len(nc_files_RTOFS)):
        ncRTOFS = xr.open_dataset(folder_RTOFS + fol + '/' + nc_files_RTOFS[t])
        tRTOFS.append(np.asarray(ncRTOFS.MT[:])[0])
    
    tRTOFS = np.asarray([mdates.num2date(mdates.date2num(tRTOFS[t])) \
              for t in np.arange(len(nc_files_RTOFS))])
    
    #%% Get RTOFS fields
    
    oklonRTOFS = np.where(np.logical_and(lonRTOFS[0,:] >= lon_lim[0],lonRTOFS[0,:] <= lon_lim[1]))[0]
    oklatRTOFS = np.where(np.logical_and(latRTOFS[:,0] >= lat_lim[0],latRTOFS[:,0] <= lat_lim[1]))[0]
    
    i=0
    nc_file = folder_RTOFS + fol + '/' + nc_files_RTOFS[i]
    ncRTOFS = xr.open_dataset(nc_file)
    time_RTOFS = tRTOFS[i]
    lon_RTOFS = lonRTOFS[0,oklonRTOFS]
    lat_RTOFS = latRTOFS[oklatRTOFS,0]
    sst_RTOFS = np.asarray(ncRTOFS.variables['temperature'][i,0,oklatRTOFS,oklonRTOFS])
    sss_RTOFS = np.asarray(ncRTOFS.variables['salinity'][i,0,oklatRTOFS,oklonRTOFS])
    su_RTOFS = np.asarray(ncRTOFS.variables['u'][i,0,oklatRTOFS,oklonRTOFS])
    sv_RTOFS = np.asarray(ncRTOFS.variables['v'][i,0,oklatRTOFS,oklonRTOFS])
    
    #%% SST
    
    kw = dict(levels = np.linspace(24,30,16))
    
    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,sst_RTOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='k')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Oper. SST on '+str(time_RTOFS)[0:13],fontsize=16)
        
    file = folder_fig +'RTOFS_SST'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% SSS
    
    kw = dict(levels = np.linspace(35,37,11))
    
    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,sss_RTOFS,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='k')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Oper. SSS on '+str(time_RTOFS)[0:13],fontsize=16)
        
    file = folder_fig +'RTOFS_SSS'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Surfeace velocity
    
    kw = dict(levels = np.linspace(24,30,16))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,sst_RTOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SST and Surface Vel. on '+str(time_RTOFS)[0:13],fontsize=16)
    q=plt.quiver(lon_RTOFS[::5],lat_RTOFS[::5],su_RTOFS[::5,::5],sv_RTOFS[::5,::5] ,scale=3,scale_units='inches',\
              alpha=0.7)
    plt.quiverkey(q,np.max(lon_RTOFS)+5,np.max(lat_RTOFS)+0.4,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file = folder_fig + 'RTOFS_SST_UV'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure sst at 200 meters
    
    okdepth = np.where(depth_RTOFS >= 200)[0][0]
    temp_200_RTOFS = np.asarray(ncRTOFS.variables['temperature'][i,okdepth,oklatRTOFS,oklonRTOFS])
    
    kw = dict(levels = np.linspace(10,25,31))
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,temp_200_RTOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Temperature at 200 m on '+str(time_RTOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'RTOFS_temp_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure sss at 200 meters
    
    okdepth = np.where(depth_RTOFS >= 200)[0][0]
    salt_200_RTOFS = np.asarray(ncRTOFS.variables['salinity'][i,okdepth,oklatRTOFS,oklonRTOFS])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,salt_200_RTOFS,cmap=cmocean.cm.haline)#,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Salinity at 200 m on '+str(time_RTOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'RTOFS_salt_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp transect along storm path
    
    lon_forec_track_interp = np.interp(lat_RTOFS,lat_forec_track,lon_forec_track)
    lat_forec_track_interp = lat_RTOFS
    
    oklon = np.round(np.interp(lon_forec_track_interp,lonRTOFS[0,:],np.arange(len(lonRTOFS[0,:])))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_interp,latRTOFS[:,0],np.arange(len(latRTOFS[:,0])))).astype(int)
    
    i=0
    trans_temp_RTOFS = np.empty((len(depth_RTOFS),len(lon_forec_track_interp)))
    trans_temp_RTOFS[:] = np.nan
    for x in np.arange(len(lon_forec_track_interp)):
        trans_temp_RTOFS[:,x] = np.asarray(ncRTOFS.variables['temperature'][i,:,oklat[x],oklon[x]])

    ok12 = np.where(trans_temp_RTOFS <= 12)[0][0]
    max_depth = -depth_RTOFS[ok12]
    
    kw = dict(levels = np.linspace(12,32,21))
    
    plt.figure()
    plt.contourf(latRTOFS[oklat,0],-depth_RTOFS,trans_temp_RTOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(latRTOFS[oklat,0],-depth_RTOFS,trans_temp_RTOFS,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.title('RTOFS Temperature \n along Forecasted Storm Track',fontsize=16)
    plt.ylim([max_depth,0])
    #plt.xlim([20,30])
    
    file = folder_fig + 'RTOFS_temp_along_forecasted_track_'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
