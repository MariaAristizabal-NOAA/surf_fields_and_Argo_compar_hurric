
"""
Created on Thu Jun 11 13:32:41 2020

@author: aristizabal
"""
def GOFS_RTOFS_vs_Argo_floats(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,folder_fig):
    #%% User input
    
    #GOFS3.1 output model location
    url_GOFS_ts = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z'
    
    
    # RTOFS files
    folder_RTOFS = '/home/coolgroup/RTOFS/forecasts/domains/hurricanes/RTOFS_6hourly_North_Atlantic/'
    
    nc_files_RTOFS = ['rtofs_glo_3dz_f006_6hrly_hvr_US_east.nc',\
                      'rtofs_glo_3dz_f012_6hrly_hvr_US_east.nc',\
                      'rtofs_glo_3dz_f018_6hrly_hvr_US_east.nc',\
                      'rtofs_glo_3dz_f024_6hrly_hvr_US_east.nc']
        
    # COPERNICUS MARINE ENVIRONMENT MONITORING SERVICE (CMEMS)
    url_cmems = 'http://nrt.cmems-du.eu/motu-web/Motu'
    service_id = 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'
    product_id = 'global-analysis-forecast-phy-001-024'
    depth_min = '0.493'
    out_dir = '/home/aristizabal/crontab_jobs'
        
    # Bathymetry file
    #bath_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/nc_files/GEBCO_2014_2D_-100.0_0.0_-60.0_45.0.nc'    
    bath_file = '/home/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    
    # Argo floats
    url_Argo = 'http://www.ifremer.fr/erddap'
    
    #%%
    
    from matplotlib import pyplot as plt
    import numpy as np
    import xarray as xr
    import netCDF4 
    from datetime import datetime, timedelta
    import cmocean
    import matplotlib.dates as mdates
    from erddapy import ERDDAP
    import pandas as pd
    import os
    
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
    #ti = datetime.today()
    ti = datetime.today() - timedelta(1) - timedelta(hours=6)
    tini = datetime(ti.year,ti.month,ti.day)
    te = ti + timedelta(2)
    tend = datetime(te.year,te.month,te.day)

    #%% Look for Argo datasets 
    
    e = ERDDAP(server = url_Argo)
    
    # Grab every dataset available
    #datasets = pd.read_csv(e.get_search_url(response='csv', search_for='all'))
    
    kw = {
        'min_lon': lon_lim[0],
        'max_lon': lon_lim[1],
        'min_lat': lat_lim[0],
        'max_lat': lat_lim[1],
        'min_time': str(tini),
        'max_time': str(tend),
    }
    
    search_url = e.get_search_url(response='csv', **kw)
    
    # Grab the results
    search = pd.read_csv(search_url)
    
    # Extract the IDs
    dataset = search['Dataset ID'].values
    
    msg = 'Found {} Datasets:\n\n{}'.format
    print(msg(len(dataset), '\n'.join(dataset)))
    
    dataset_type = dataset[0]
    
    constraints = {
        'time>=': str(tini),
        'time<=': str(tend),
        'latitude>=': lat_lim[0],
        'latitude<=': lat_lim[1],
        'longitude>=':lon_lim[0],
        'longitude<=': lon_lim[1],
    }
    
    variables = [
     'platform_number',  
     'time',
     'pres',
     'longitude',
     'latitude', 
     'temp',
     'psal',
    ]
    
    e = ERDDAP(
        server = url_Argo,
        protocol = 'tabledap',
        response = 'nc'
    )
    
    e.dataset_id = dataset_type
    e.constraints=constraints
    e.variables=variables
    
    print(e.get_download_url())
    
    df = e.to_pandas(
        parse_dates=True,
        skiprows=(1,)  # units information can be dropped.
    ).dropna()
    
    argo_ids = np.asarray(df['platform_number'])
    argo_times = np.asarray(df['time (UTC)'])
    argo_press = np.asarray(df['pres (decibar)'])
    argo_lons = np.asarray(df['longitude (degrees_east)'])
    argo_lats = np.asarray(df['latitude (degrees_north)'])
    argo_temps = np.asarray(df['temp (degree_Celsius)'])
    argo_salts = np.asarray(df['psal (PSU)'])
    
    #%% GOGF 3.1
        
    try:
        GOFS_ts = xr.open_dataset(url_GOFS_ts,decode_times=False)
        
        lt_GOFS = np.asarray(GOFS_ts['lat'][:])
        ln_GOFS = np.asarray(GOFS_ts['lon'][:])
        tt = GOFS_ts['time']
        t_GOFS = netCDF4.num2date(tt[:],tt.units) 
        depth_GOFS = np.asarray(GOFS_ts['depth'][:])
    except Exception as err:
        print(err)
        GOFS_ts = np.nan
        lt_GOFS = np.nan
        ln_GOFS = np.nan
        depth_GOFS = np.nan
        t_GOFS = ti
    
    #%% Map Argo floats
     
    lev = np.arange(-9000,9100,100)
    plt.figure()
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo) 
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='gold')
    plt.plot(lon_forec_cone,lat_forec_cone,'.-b',markersize=1)
    plt.plot(lon_best_track,lat_best_track,'or',markersize=3)
    
    argo_idd = np.unique(argo_ids)
    for i,id in enumerate(argo_idd): 
        okind = np.where(argo_ids == id)[0]
        plt.plot(np.unique(argo_lons[okind]),np.unique(argo_lats[okind]),'s',color='darkorange',markersize=5,markeredgecolor='k')
    
    plt.title('Argo Floats ' + str(tini)[0:13]+'-'+str(tend)[0:13],fontsize=16)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    
    file = folder_fig + 'ARGO_lat_lon'
    #file = folder_fig + 'ARGO_lat_lon_' + str(np.unique(argo_times)[0])[0:10]
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 

    #%% Figure argo float vs GOFS and vs RTOFS
    
    argo_idd = np.unique(argo_ids)
    
    for i,id in enumerate(argo_idd): 
        print(id)
        okind = np.where(argo_ids == id)[0]
        argo_time = np.asarray([datetime.strptime(t,'%Y-%m-%dT%H:%M:%SZ') for t in argo_times[okind]])
        
        argo_lon = argo_lons[okind]
        argo_lat = argo_lats[okind]
        argo_pres = argo_press[okind]
        argo_temp = argo_temps[okind]
        argo_salt = argo_salts[okind]
        
        # GOFS
        print('Retrieving variables from GOFS')
        if isinstance(GOFS_ts,float): 
            temp_GOFS = np.nan
            salt_GOFS = np.nan
        else:
            #oktt_GOFS = np.where(t_GOFS >= argo_time[0])[0][0]
            ttGOFS = np.asarray([datetime(t_GOFS[i].year,t_GOFS[i].month,t_GOFS[i].day,t_GOFS[i].hour) for i in np.arange(len(t_GOFS))])
            tstamp_GOFS = [mdates.date2num(ttGOFS[i]) for i in np.arange(len(ttGOFS))]
            oktt_GOFS = np.unique(np.round(np.interp(mdates.date2num(argo_time[0]),tstamp_GOFS,np.arange(len(tstamp_GOFS)))).astype(int))[0]
            oklat_GOFS = np.where(lt_GOFS >= argo_lat[0])[0][0]
            oklon_GOFS = np.where(ln_GOFS >= argo_lon[0]+360)[0][0]
            temp_GOFS = np.asarray(GOFS_ts['water_temp'][oktt_GOFS,:,oklat_GOFS,oklon_GOFS])
            salt_GOFS = np.asarray(GOFS_ts['salinity'][oktt_GOFS,:,oklat_GOFS,oklon_GOFS])
        
        # RTOFS 
        #Time window
        year = int(argo_time[0].year)
        month = int(argo_time[0].month)
        day = int(argo_time[0].day)
        tini = datetime(year, month, day)
        tend = tini + timedelta(days=1)
        
        # Read RTOFS grid and time
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
        
        oktt_RTOFS = np.where(mdates.date2num(tRTOFS) >= mdates.date2num(argo_time[0]))[0][0]
        oklat_RTOFS = np.where(latRTOFS[:,0] >= argo_lat[0])[0][0]
        oklon_RTOFS = np.where(lonRTOFS[0,:]  >= argo_lon[0])[0][0]
            
        nc_file = folder_RTOFS + fol + '/' + nc_files_RTOFS[oktt_RTOFS]
        ncRTOFS = xr.open_dataset(nc_file)
        #time_RTOFS = tRTOFS[oktt_RTOFS]
        temp_RTOFS = np.asarray(ncRTOFS.variables['temperature'][0,:,oklat_RTOFS,oklon_RTOFS])
        salt_RTOFS = np.asarray(ncRTOFS.variables['salinity'][0,:,oklat_RTOFS,oklon_RTOFS])
        #lon_RTOFS = lonRTOFS[0,oklon_RTOFS]
        #lat_RTOFS = latRTOFS[oklat_RTOFS,0]
        
        # Downloading and reading Copernicus output
        motuc = 'python -m motuclient --motu ' + url_cmems + \
        ' --service-id ' + service_id + \
        ' --product-id ' + product_id + \
        ' --longitude-min ' + str(argo_lon[0]-2/12) + \
        ' --longitude-max ' + str(argo_lon[0]+2/12) + \
        ' --latitude-min ' + str(argo_lat[0]-2/12) + \
        ' --latitude-max ' + str(argo_lat[0]+2/12) + \
        ' --date-min ' + '"' + str(tini-timedelta(0.5)) + '"' + \
        ' --date-max ' + '"' + str(tend+timedelta(0.5)) + '"' + \
        ' --depth-min ' + depth_min + \
        ' --depth-max ' + str(np.nanmax(argo_pres)+1000) + \
        ' --variable ' + 'thetao' + ' ' + \
        ' --variable ' + 'so'  + ' ' + \
        ' --out-dir ' + out_dir + \
        ' --out-name ' + str(id) + '.nc' + ' ' + \
        ' --user ' + 'maristizabalvar' + ' ' + \
        ' --pwd ' +  'MariaCMEMS2018'
        
        os.system(motuc)
        # Check if file was downloaded

        COP_file = out_dir + '/' + str(id) + '.nc'
        # Check if file was downloaded
        resp = os.system('ls ' + out_dir +'/' + str(id) + '.nc')
        if resp == 0:
            COP = xr.open_dataset(COP_file)

            latCOP = np.asarray(COP.latitude[:])
            lonCOP = np.asarray(COP.longitude[:])
            depth_COP = np.asarray(COP.depth[:])
            tCOP = np.asarray(mdates.num2date(mdates.date2num(COP.time[:])))
        else:
            latCOP = np.empty(1)
            latCOP[:] = np.nan
            lonCOP = np.empty(1)
            lonCOP[:] = np.nan
            tCOP = np.empty(1)
            tCOP[:] = np.nan

        oktimeCOP = np.where(mdates.date2num(tCOP) >= mdates.date2num(tini))[0][0]
        oklonCOP = np.where(lonCOP >= argo_lon[0])[0][0]
        oklatCOP = np.where(latCOP >= argo_lat[0])[0][0]
        
        temp_COP = np.asarray(COP.variables['thetao'][oktimeCOP,:,oklatCOP,oklonCOP])
        salt_COP = np.asarray(COP.variables['so'][oktimeCOP,:,oklatCOP,oklonCOP])
        
        # Figure temp
        plt.figure(figsize=(5,6))
        plt.plot(argo_temp,-argo_pres,'.-',linewidth=2,label='ARGO Float id '+str(id))
        plt.plot(temp_GOFS,-depth_GOFS,'.-',linewidth=2,label='GOFS 3.1',color='red')
        plt.plot(temp_RTOFS,-depth_RTOFS,'.-',linewidth=2,label='RTOFS',color='g')
        plt.plot(temp_COP,-depth_COP,'.-',linewidth=2,label='Copernicus',color='darkorchid')
        plt.ylim([-1000,0])
        plt.title('Temperature Profile on '+ str(argo_time[0])[0:13] +
                  '\n [lon,lat] = [' \
                  + str(np.round(argo_lon[0],3)) +',' +\
                      str(np.round(argo_lat[0],3))+']',\
                      fontsize=16)
        plt.ylabel('Depth (m)',fontsize=14)
        plt.xlabel('$^oC$',fontsize=14)
        plt.legend(loc='lower right',fontsize=14)
        
        file = folder_fig + 'ARGO_vs_GOFS_RTOFS_COP_temp_' + str(id) 
        plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
        
        # Figure salt
        plt.figure(figsize=(5,6))
        plt.plot(argo_salt,-argo_pres,'.-',linewidth=2,label='ARGO Float id '+str(id))
        plt.plot(salt_GOFS,-depth_GOFS,'.-',linewidth=2,label='GOFS 3.1',color='red')
        plt.plot(salt_RTOFS,-depth_RTOFS,'.-',linewidth=2,label='RTOFS',color='g')
        plt.plot(salt_COP,-depth_COP,'.-',linewidth=2,label='Copernicus',color='darkorchid')
        plt.ylim([-1000,0])
        plt.title('Salinity Profile on '+ str(argo_time[0])[0:13] +
                  '\n [lon,lat] = [' \
                  + str(np.round(argo_lon[0],3)) +',' +\
                      str(np.round(argo_lat[0],3))+']',\
                      fontsize=16)
        plt.ylabel('Depth (m)',fontsize=14)
        plt.legend(loc='lower right',fontsize=14)
        
        file = folder_fig + 'ARGO_vs_GOFS_RTOFS_COP_salt_' + str(id)
        plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 



