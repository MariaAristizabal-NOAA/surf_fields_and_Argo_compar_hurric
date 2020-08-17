
"""
Created on Thu Jun  4 13:11:41 2020

@author: aristizabal
"""

def GOFS31_baffin(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig):
    #%% User input
    
    #GOFS3.1 output model location
    url_GOFS_ts = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z'
    
    url_GOFS_uv = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z'
    
    url_GOFS_ssh = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ssh'
    
    # Bathymetry file
    bath_file = '/home/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    
    #%%
    from matplotlib import pyplot as plt
    import numpy as np
    import xarray as xr
    import netCDF4 
    from datetime import datetime, timedelta
    import cmocean
    import matplotlib.dates as mdates
    
    # Do not produce figures on screen
    plt.switch_backend('agg')
    
    # Increase fontsize of labels globally
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('legend',fontsize=14)
    
    #%% Get time bounds for current day
    #ti = datetime.today() 
    ti = datetime.today() - timedelta(hours=6)
    tini = datetime(ti.year,ti.month,ti.day,12)

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

    ttGOFS = np.asarray([datetime(t_GOFS[i].year,t_GOFS[i].month,t_GOFS[i].day,t_GOFS[i].hour) for i in np.arange(len(t_GOFS))])
    tstamp_GOFS = [mdates.date2num(ttGOFS[i]) for i in np.arange(len(ttGOFS))]
    oktime_GOFS = np.unique(np.round(np.interp(mdates.date2num(tini),tstamp_GOFS,np.arange(len(tstamp_GOFS)))).astype(int))
    time_GOFS = ttGOFS[oktime_GOFS][0]
    #oktime_GOFS = np.where(t_GOFS == tini)[0][0]
    
    # Conversion from GOFS convention to glider longitude and latitude
    ln_GOFSg= np.empty((len(ln_GOFS),))
    ln_GOFSg[:] = np.nan
    for i in range(len(ln_GOFS)):
        if ln_GOFS[i] > 180: 
            ln_GOFSg[i] = ln_GOFS[i] - 360 
        else:
            ln_GOFSg[i] = ln_GOFS[i]
    lt_GOFSg = lt_GOFS
    
    lat_GOFS= lt_GOFS[oklat_GOFS]
    lon_GOFS= ln_GOFS[oklon_GOFS]
    lon_GOFSg= ln_GOFSg[oklon_GOFS]
    lat_GOFSg= lt_GOFSg[oklat_GOFS]
    #time_GOFS = t_GOFS[oktime_GOFS]
       
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
    sst_GOFS = np.asarray(GOFS_ts['water_temp'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])[0,:,:]
    sss_GOFS = np.asarray(GOFS_ts['salinity'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])[0,:,:]
    ssh_GOFS = np.asarray(GOFS_ssh['surf_el'][oktime_GOFS,oklat_GOFS,oklon_GOFS])[0,:,:]
    su_GOFS = np.asarray(GOFS_uv['water_u'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])[0,:,:]
    sv_GOFS = np.asarray(GOFS_uv['water_v'][oktime_GOFS,0,oklat_GOFS,oklon_GOFS])[0,:,:]

    #%% Figure sst
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,sst_GOFS[:,:],cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SST \n on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig +'GOFS_SST'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure sss
    kw = dict(levels = np.arange(salt_lim[0],salt_lim[1],0.5))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,sss_GOFS,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SSS \n on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_SSS'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure ssh
    kw = dict(levels = np.arange(-1.0,1.1,0.1))    

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,ssh_GOFS,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14) 
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SSH \n on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_SSH'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure ssh
    kw = dict(levels = np.arange(-1.0,1.1,0.1))    

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,ssh_GOFS,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14) 
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS SSH \n on '+str(time_GOFS)[0:13],fontsize=16)
    q=plt.quiver(lon_GOFSg[::7],lat_GOFSg[::7],su_GOFS[::7,::7],sv_GOFS[::7,::7] ,scale=3,scale_units='inches',\
              alpha=0.7)
    plt.quiverkey(q,np.max(lon_GOFSg)-0.2,np.max(lat_GOFSg)+0.5,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file = folder_fig + 'GOFS_SSH_UV'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp at 200 meters
    kw = dict(levels = np.arange(temp200_lim[0],temp200_lim[1],1))
    okdepth = np.where(depth_GOFS >= 200)[0][0]
    temp_200_GOFS = np.asarray(GOFS_ts['water_temp'][oktime_GOFS,okdepth,oklat_GOFS,oklon_GOFS])[0,:,:]
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,temp_200_GOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS Temperature at 200 m \n on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_temp_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure salinity at 200 meters
    kw = dict(levels = np.arange(salt200_lim[0],salt200_lim[1],0.3))
    okdepth = np.where(depth_GOFS >= 200)[0][0]
    salt_200_GOFS = np.asarray(GOFS_ts['salinity'][oktime_GOFS,okdepth,oklat_GOFS,oklon_GOFS])[0,:,:]
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,salt_200_GOFS,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS 3.1 Salinity at 200 m \n on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_salt_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 

    #%% Figure bottom temp 
    temp_bott_GOFS = np.asarray(GOFS_ts['water_temp_bottom'][oktime_GOFS,oklat_GOFS,oklon_GOFS])[0,:,:]
    #kw = dict(levels = np.arange(0,np.nanmax(temp_bott_GOFS),1))
    kw = dict(levels = np.arange(tempb_lim[0],tempb_lim[1],1))    

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_GOFSg,lat_GOFSg,temp_bott_GOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('GOFS Bottom Temperature \n on '+str(time_GOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'GOFS_temp_bott'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp transect along storm path
    lon_forec_track_interp = np.interp(lt_GOFS,lat_forec_track,lon_forec_track,left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lt_GOFS)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    
    oklon = np.round(np.interp(lon_forec_track_int+360,ln_GOFS,np.arange(len(ln_GOFS)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lt_GOFS,np.arange(len(lt_GOFS)))).astype(int)
    okdepth = np.where(depth_GOFS <= 350)[0]
    
    #i=0
    trans_temp_GOFS = np.empty((len(depth_GOFS[okdepth]),len(lon_forec_track_int)))
    trans_temp_GOFS[:] = np.nan
    for x in np.arange(len(lon_forec_track_int)):
        trans_temp_GOFS[:,x] = np.asarray(GOFS_ts['water_temp'][oktime_GOFS,okdepth,oklat[x],oklon[x]])
    
    kw = dict(levels = np.arange(tempt_lim[0],tempt_lim[1],1))
    
    plt.figure()
    plt.contourf(lt_GOFS[oklat],-depth_GOFS[okdepth],trans_temp_GOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(lt_GOFSg[oklat],-depth_GOFS[okdepth],trans_temp_GOFS,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Latitude ($^o$)',fontsize=14)
    plt.title('GOFS Temperature \n along Forecasted Storm Track',fontsize=16)
    
    file = folder_fig + 'GOFS_temp_along_forecasted_track_'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
def RTOFS_oper_baffin(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig):

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
    
    #%% Get time bounds for current day
    ti = datetime.today() 
    tini = datetime(ti.year,ti.month,ti.day)
    #tini = datetime(ti.year,ti.month,ti.day,12)
    
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
    
    ncRTOFS = xr.open_dataset(folder_RTOFS + fol + '/' + nc_files_RTOFS[1])
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
    
    t=1
    nc_file = folder_RTOFS + fol + '/' + nc_files_RTOFS[t]
    ncRTOFS = xr.open_dataset(nc_file)
    time_RTOFS = tRTOFS[t]
    lon_RTOFS = lonRTOFS[0,oklonRTOFS]
    lat_RTOFS = latRTOFS[oklatRTOFS,0]
    sst_RTOFS = np.asarray(ncRTOFS.variables['temperature'][0,0,oklatRTOFS,oklonRTOFS])
    sss_RTOFS = np.asarray(ncRTOFS.variables['salinity'][0,0,oklatRTOFS,oklonRTOFS])
    su_RTOFS = np.asarray(ncRTOFS.variables['u'][0,0,oklatRTOFS,oklonRTOFS])
    sv_RTOFS = np.asarray(ncRTOFS.variables['v'][0,0,oklatRTOFS,oklonRTOFS])

    #%% SST
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))
    
    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,sst_RTOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Oper. SST \n on '+str(time_RTOFS)[0:13],fontsize=16)
        
    file = folder_fig +'RTOFS_SST'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% SSS
    kw = dict(levels = np.arange(salt_lim[0],salt_lim[1],0.5))
    
    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,sss_RTOFS,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Oper. SSS \n on '+str(time_RTOFS)[0:13],fontsize=16)
        
    file = folder_fig +'RTOFS_SSS'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Surface velocity
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,sst_RTOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS SST \n on '+str(time_RTOFS)[0:13],fontsize=16)
    q=plt.quiver(lon_RTOFS[::6],lat_RTOFS[::6],su_RTOFS[::6,::6],sv_RTOFS[::6,::6] ,scale=3,scale_units='inches',\
              alpha=0.7)
    plt.quiverkey(q,np.max(lon_RTOFS)-0.2,np.max(lat_RTOFS)+0.5,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file = folder_fig + 'RTOFS_SST_UV'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp at 200 meters
    kw = dict(levels = np.arange(temp200_lim[0],temp200_lim[1],1))
    okdepth = np.where(depth_RTOFS >= 200)[0][0]
    temp_200_RTOFS = np.asarray(ncRTOFS.variables['temperature'][0,okdepth,oklatRTOFS,oklonRTOFS])

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,temp_200_RTOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Temperature at 200 m \n on '+str(time_RTOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'RTOFS_temp_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure salt at 200 meters
    kw = dict(levels = np.arange(salt200_lim[0],salt200_lim[1],0.3))
    okdepth = np.where(depth_RTOFS >= 200)[0][0]
    salt_200_RTOFS = np.asarray(ncRTOFS.variables['salinity'][0,okdepth,oklatRTOFS,oklonRTOFS])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,salt_200_RTOFS,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Salinity at 200 m \n on '+str(time_RTOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'RTOFS_salt_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 

    #%% Figure temp bottom 
    temp_RTOFS = np.asarray(ncRTOFS.variables['temperature'][0,:,oklatRTOFS,oklonRTOFS])
    
    temp_bott_RTOFS = np.empty((temp_RTOFS.shape[1],temp_RTOFS.shape[2]))
    temp_bott_RTOFS[:] = np.nan
    for x in np.arange(temp_RTOFS.shape[1]):
        for y in np.arange(temp_RTOFS.shape[2]):
            okb = np.where(np.isfinite(temp_RTOFS[:,x,y]))[0]
            if len(okb) != 0:
                temp_bott_RTOFS[x,y] = temp_RTOFS[okb[-1],x,y]
            else:
                temp_bott_RTOFS[x,y] = np.nan
                
    kw = dict(levels = np.arange(tempb_lim[0],tempb_lim[1],1))    
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_RTOFS,lat_RTOFS,temp_bott_RTOFS,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('RTOFS Bottom Temperature \n on '+str(time_RTOFS)[0:13],fontsize=16)
    
    file = folder_fig + 'RTOFS_temp_bott'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp transect along storm path
    lon_forec_track_interp = np.interp(latRTOFS[:,0],lat_forec_track,lon_forec_track,left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(latRTOFS[:,0])
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    
    oklon = np.round(np.interp(lon_forec_track_int,lonRTOFS[0,:],np.arange(len(lonRTOFS[0,:])))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,latRTOFS[:,0],np.arange(len(latRTOFS[:,0])))).astype(int)
    okdepth = np.where(depth_RTOFS <= 350)[0]
    
    trans_temp_RTOFS = np.empty((len(depth_RTOFS[okdepth]),len(lon_forec_track_int)))
    trans_temp_RTOFS[:] = np.nan
    for x in np.arange(len(lon_forec_track_int)):
        trans_temp_RTOFS[:,x] = np.asarray(ncRTOFS.variables['temperature'][0,okdepth,oklat[x],oklon[x]])
       
    kw = dict(levels = np.arange(tempt_lim[0],tempt_lim[1],1))
    
    plt.figure()
    plt.contourf(latRTOFS[oklat,0],-depth_RTOFS[okdepth],trans_temp_RTOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(latRTOFS[oklat,0],-depth_RTOFS[okdepth],trans_temp_RTOFS,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Latitude ($^o$)',fontsize=14)
    plt.title('RTOFS Temperature \n along Forecasted Storm Track',fontsize=16)
    
    file = folder_fig + 'RTOFS_temp_along_forecasted_track_'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
   
#%%%%%%%%%%%%%%%%%%%%
def Copernicus_baffin(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig):

    #%%
    # COPERNICUS MARINE ENVIRONMENT MONITORING SERVICE (CMEMS)
    url_cmems = 'http://nrt.cmems-du.eu/motu-web/Motu'
    service_id = 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'
    product_id = 'global-analysis-forecast-phy-001-024'
    depth_min = '0.493'
    out_dir = '/home/aristizabal/crontab_jobs'
        
    # Bathymetry file
    bath_file = '/home/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    
    #%%
    from matplotlib import pyplot as plt
    import numpy as np
    import xarray as xr
    from datetime import datetime, timedelta
    import cmocean
    import os
    import netCDF4
    
    # Do not produce figures on screen
    plt.switch_backend('agg')
    
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
    
    #%% Get time bounds for current day
    ti = datetime.today() 
    tini = datetime(ti.year,ti.month,ti.day,12)
    te= ti + timedelta(1) 
    tend = datetime(te.year,te.month,te.day,12)
    
    #%% Downloading and reading Copernicus output
    motuc = 'python -m motuclient --motu ' + url_cmems + \
        ' --service-id ' + service_id + \
        ' --product-id ' + product_id + \
        ' --longitude-min ' + str(lon_lim[0]-5) + \
        ' --longitude-max ' + str(lon_lim[1]+5) + \
        ' --latitude-min ' + str(lat_lim[0]-5) + \
        ' --latitude-max ' + str(lat_lim[1]+5) + \
        ' --date-min ' + str(tini-timedelta(0.5)) + \
        ' --date-max ' + str(tend+timedelta(0.5)) + \
        ' --depth-min ' + depth_min + \
        ' --depth-max ' + '350' + \
        ' --variable ' + 'thetao' + \
        ' --variable ' + 'so' + \
        ' --variable ' + 'zos' + \
        ' --variable ' + 'uo' + \
        ' --variable ' + 'vo' + \
        ' --variable ' + 'bottomT' + \
        ' --out-dir ' + out_dir + \
        ' --out-name ' + folder_fig.split('/')[-2] + '.nc' + ' ' + \
        ' --user ' + 'maristizabalvar' + ' ' + \
        ' --pwd ' +  'MariaCMEMS2018'
        
    os.system(motuc)
    # Check if file was downloaded

    COP_file = out_dir + '/' + folder_fig.split('/')[-2] + '.nc'
    # Check if file was downloaded
    resp = os.system('ls ' + out_dir +'/' + folder_fig.split('/')[-2] + '.nc')
    if resp == 0:
        COP = xr.open_dataset(COP_file,decode_times=False)

        lat_COP = np.asarray(COP.latitude[:])
        lon_COP = np.asarray(COP.longitude[:])
        depth_COP = np.asarray(COP.depth[:])
        tt = COP.time[:]
        tCOP = netCDF4.num2date(tt[:],tt.units) 
        
        oktimeCOP = np.where(tCOP >= tini)[0][0]    
        time_COP = tCOP[oktimeCOP]
        sst_COP = np.asarray(COP.variables['thetao'][oktimeCOP,0,:,:])
        sss_COP = np.asarray(COP.variables['so'][oktimeCOP,0,:,:])
        ssh_COP = np.asarray(COP.variables['zos'][oktimeCOP,:,:])
        su_COP = np.asarray(COP.variables['uo'][oktimeCOP,0,:,:])
        sv_COP = np.asarray(COP.variables['vo'][oktimeCOP,0,:,:])
        temp_bott_COP = np.asarray(COP.variables['bottomT'][oktimeCOP,:,:])
    else:
        lat_COP = np.empty(1)
        lat_COP[:] = np.nan
        lon_COP = np.empty(1)
        lon_COP[:] = np.nan
        time_COP = np.empty(1)
        time_COP[:] = np.nan
        sst_COP = np.empty(1)
        sst_COP[:] = np.nan
        sss_COP = np.empty(1)
        sss_COP[:] = np.nan
        ssh_COP = np.empty(1)
        ssh_COP[:] = np.nan
        su_COP = np.empty(1)
        su_COP[:] = np.nan
        sv_COP = np.empty(1)
        sv_COP[:] = np.nan
        temp_bott_COP = np.empty(1)
        temp_bott_COP[:] = np.nan
    
    #%% SST
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))
    
    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,sst_COP,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus SST \n on '+str(time_COP)[0:13],fontsize=16)
        
    file = folder_fig +'COP_SST'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% SSS
    kw = dict(levels = np.arange(salt_lim[0],salt_lim[1],0.5))
    
    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,sss_COP,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus SSS \n on '+str(time_COP)[0:13],fontsize=16)
        
    file = folder_fig +'COP_SSS'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure ssh
    kw = dict(levels = np.arange(-1,1.1,0.1))    

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,ssh_COP,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14) 
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus SSH \n on '+str(time_COP)[0:13],fontsize=16)
    
    file = folder_fig + 'COP_SSH'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure ssh
    kw = dict(levels = np.arange(-1,1.1,0.1))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,ssh_COP,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14) 
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus SSH \n on '+str(time_COP)[0:13],fontsize=16)
    q=plt.quiver(lon_COP[::6],lat_COP[::6],su_COP[::6,::6],sv_COP[::6,::6] ,scale=3,scale_units='inches',\
              alpha=0.7)
    plt.quiverkey(q,np.max(lon_COP)-5,np.max(lat_COP)-4,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file = folder_fig + 'COP_SSH_UV'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure temp at 200 meters
    kw = dict(levels = np.arange(temp200_lim[0],temp200_lim[1],1))
    okdepth = np.where(depth_COP >= 200)[0][0]
    temp_200_COP = np.asarray(COP.variables['thetao'][oktimeCOP,okdepth,:,:])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,temp_200_COP,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus Temperature at 200 m \n on '+str(time_COP)[0:13],fontsize=16)
    
    file = folder_fig + 'COP_temp_at_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
    
    #%% Figure salinity at 200 meters
    kw = dict(levels = np.arange(salt200_lim[0],salt200_lim[1],0.3))
    okdepth = np.where(depth_COP >= 200)[0][0]
    salt_200_COP = np.asarray(COP.variables['so'][oktimeCOP,okdepth,:,:])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,salt_200_COP,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus Salinity at 200 m \n on '+str(time_COP)[0:13],fontsize=16)
    
    file = folder_fig + 'COP_salt_200m'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 

    #%% Bottom temperature
    kw = dict(levels = np.arange(tempb_lim[0],tempb_lim[1],1))

    fig, ax = plt.subplots()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_COP,lat_COP,temp_bott_COP,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('Copernicus Bottom Temperature \n on '+str(time_COP)[0:13],fontsize=16)

    file = folder_fig +'COP_temp_bott'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Figure temp transect along storm path
    lon_forec_track_interp = np.interp(lat_COP,lat_forec_track,lon_forec_track,left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_COP)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    
    oklon = np.round(np.interp(lon_forec_track_int,lon_COP,np.arange(len(lon_COP)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_COP,np.arange(len(lat_COP)))).astype(int)
    okdepth = np.where(depth_COP <= 350)[0]
    
    trans_temp_GOFS = np.empty((len(depth_COP[okdepth]),len(lon_forec_track_int)))
    trans_temp_GOFS[:] = np.nan
    for i in np.arange(len(lon_forec_track_int)):
        trans_temp_GOFS[:,i] = np.asarray(COP['thetao'][oktimeCOP,okdepth,oklat[i],oklon[i]])
    
    kw = dict(levels = np.arange(tempt_lim[0],tempt_lim[1],1))
    
    plt.figure()
    plt.contourf(lat_COP[oklat],-depth_COP[okdepth],trans_temp_GOFS,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(lat_COP[oklat],-depth_COP[okdepth],trans_temp_GOFS,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Latitude ($^o$)',fontsize=14)
    plt.title('Copernicus Temperature \n along Forecasted Storm Track',fontsize=16)

    file = folder_fig + 'COP_temp_along_forecasted_track'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

    
