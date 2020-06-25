
"""
Created on Mon Jun 15 14:31:20 2020

@author: aristizabal
"""
    
#%% url NHC    

url_nhc = 'https://www.nhc.noaa.gov/gis/'

# GoMex
lon_lim_GoMex = [-98,-80]
lat_lim_GoMex = [15,32.5]

# Caribbean
lon_lim_Car = [-87,-60]
lat_lim_Car = [8,30]

# SAB
lon_lim_SAB = [-82,-64]
lat_lim_SAB = [25,36]

# MAB
#lon_lim_MAB = [-77,-64]
#lat_lim_MAB = [35,42]
lon_lim_MAB = [-77,-55]
lat_lim_MAB = [35,46]

lon_lim = []
lat_lim = []

#%%

from datetime import datetime, timedelta
import numpy as np
import os
import requests
import urllib.request
from bs4 import BeautifulSoup
import glob
from zipfile import ZipFile
import sys

sys.path.append('/home/aristizabal/Code/surf_fields_and_Argo_compar_hurric')

from Surf_fields_models_baffin import GOFS31_baffin, RTOFS_oper_baffin, Copernicus_baffin
from GOFS_RTOFS_oper_vs_Argo_floats_baffin import GOFS_RTOFS_vs_Argo_floats

#%% Get time bounds for the previous day
#te = datetime.today()
#tend = datetime(te.year,te.month,te.day)

#ti = datetime.today() - timedelta(1)
ti = datetime.today() - timedelta(1)
tini = datetime(ti.year,ti.month,ti.day)

#%% Download kmz files
os.system('rm -rf *best_track*')
os.system('rm -rf *TRACK_latest*')
os.system('rm -rf *CONE_latest*')

r = requests.get(url_nhc)
data = r.text

soup = BeautifulSoup(data,"html.parser")

for i,s in enumerate(soup.find_all("a")):
    ff = s.get('href')
    if type(ff) == str:
        if np.logical_and('kmz' in ff, str(tini.year) in ff):
            if 'CONE_latest' in ff:
                file_name = ff.split('/')[3]
                print(ff, file_name)
                urllib.request.urlretrieve(url_nhc[:-4] + ff , file_name)
            if 'TRACK_latest' in ff:
                file_name = ff.split('/')[3]
                print(ff, file_name)
                urllib.request.urlretrieve(url_nhc[:-4] + ff ,file_name)
            if 'best_track' in ff:
                file_name = ff.split('/')[1]
                print(ff,file_name)
                urllib.request.urlretrieve(url_nhc + ff ,file_name)

#%%
kmz_files = glob.glob('*.kmz')

if len(kmz_files) == 0:
    name = ''
    lon_forec_track = []
    lat_forec_track = []
    lon_best_track = []
    lat_best_track = []

# NOTE: UNTAR  the .kmz FILES AND THEN RUN FOLLOWING CODE
for f in kmz_files:
    os.system('cp ' + f + ' ' + f[:-3] + 'zip')
    #os.system('mkdir ' + f[:-4])
    os.system('unzip -o ' + f + ' -d ' + f[:-4])

#%% get names zip and kml files
zip_files = glob.glob('*.zip')
zip_files = [f for f in zip_files if np.logical_or('al' in f,'AL' in f)]

#%%
for i,f in enumerate(zip_files):
    kmz = ZipFile(f, 'r')
    if 'TRACK' in f:
        kml_f = glob.glob(f[:-4]+'/*.kml')
        kml_track = kmz.open(kml_f[0].split('/')[1], 'r').read()

        # Get TRACK coordinates
        soup = BeautifulSoup(kml_track,'html.parser')

        lon_forec_track = np.empty(len(soup.find_all("point")))
        lon_forec_track[:] = np.nan
        lat_forec_track = np.empty(len(soup.find_all("point")))
        lat_forec_track[:] = np.nan
        for i,s in enumerate(soup.find_all("point")):
            print(s.get_text("coordinates"))
            lon_forec_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[0])
            lat_forec_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[1])


    else:
        if 'CONE' in f:
            kml_f = glob.glob(f[:-4]+'/*.kml')
            kml_cone = kmz.open(kml_f[0].split('/')[1], 'r').read()

            # CONE coordinates
            soup = BeautifulSoup(kml_cone,'html.parser')

            lon_forec_cone = []
            lat_forec_cone = []
            for i,s in enumerate(soup.find_all("coordinates")):
                coor = s.get_text('coordinates').split(',0')
                for st in coor[1:-1]:
                    lon_forec_cone.append(st.split(',')[0])
                    lat_forec_cone.append(st.split(',')[1])

            lon_forec_cone = np.asarray(lon_forec_cone).astype(float)
            lat_forec_cone = np.asarray(lat_forec_cone).astype(float)

        else:
            kml_f = glob.glob(f[:-4]+'/*.kml')
            kml_best_track = kmz.open(kml_f[0].split('/')[1], 'r').read()

            # best track coordinates
            soup = BeautifulSoup(kml_best_track,'html.parser')

            lon_best_track = np.empty(len(soup.find_all("point")))
            lon_best_track[:] = np.nan
            lat_best_track = np.empty(len(soup.find_all("point")))
            lat_best_track[:] = np.nan
            for i,s in enumerate(soup.find_all("point")):
                print(s.get_text("coordinates"))
                lon_best_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[0])
                lat_best_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[1])

            #get name
            for f in soup.find_all('name'):
                if 'AL' in f.get_text('name'):
                    name = f.get_text('name')
                    
#%% 
name = '03l2020 Cristobal'

lon_forec_track = np.asarray([-90,-91.1,-90,-91.3,-91.1,-91.2,-91.5,-92.0,-92.5,-93,\
                   -93.3,-92.5,-92.5,-92.4,-92.1,-91.1,-90.0])
lat_forec_track = np.arange(16,33,1)


if np.logical_and(np.mean(lon_forec_track) >= lon_lim_GoMex[0],\
                  np.mean(lon_forec_track) <= lon_lim_GoMex[1]):
    if np.logical_and(np.mean(lat_forec_track) >= lat_lim_GoMex[0],\
                  np.mean(lat_forec_track) <= lat_lim_GoMex[1]):
        lon_lim = lon_lim_GoMex
        lat_lim = lat_lim_GoMex
        
if np.logical_and(np.mean(lon_forec_track) >= lon_lim_Car[0],\
                  np.mean(lon_forec_track) <= lon_lim_Car[1]):
    if np.logical_and(np.mean(lat_forec_track) >= lat_lim_Car[0],\
                  np.mean(lat_forec_track) <= lat_lim_Car[1]):
        lon_lim = lon_lim_Car
        lat_lim = lat_lim_Car  
        
if np.logical_and(np.mean(lon_forec_track) >= lon_lim_SAB[0],\
                  np.mean(lon_forec_track) <= lon_lim_SAB[1]):
    if np.logical_and(np.mean(lat_forec_track) >= lat_lim_SAB[0],\
                  np.mean(lat_forec_track) <= lat_lim_SAB[1]):
        lon_lim = lon_lim_SAB
        lat_lim = lat_lim_SAB 
        
if np.logical_and(np.mean(lon_forec_track) >= lon_lim_MAB[0],\
                  np.mean(lon_forec_track) <= lon_lim_MAB[1]):
    if np.logical_and(np.mean(lat_forec_track) >= lat_lim_MAB[0],\
                  np.mean(lat_forec_track) <= lat_lim_MAB[1]):
        lon_lim = lon_lim_MAB
        lat_lim = lat_lim_MAB 

#%% 

if np.logical_and(len(name) != 0,len(lon_lim) != 0):
    os.chdir('/www/web/rucool/hurricane/Hurricane_season_' + str(tini.year))
    os.system('mkdir ' + ti.strftime('%b-%d')  ) 
    os.system('mkdir ' + ti.strftime('%b-%d') + '/' + name.split(' ')[1] )

    folder_fig = '/www/web/rucool/hurricane/Hurricane_season_' \
        + str(tini.year)+ '/' + ti.strftime('%b-%d') + '/' + name.split(' ')[1] + '/'

    try:
        print('Reading Argo floats')
        GOFS_RTOFS_vs_Argo_floats(lon_lim,lat_lim,folder_fig)
    except Exception as err:
        print(err)
   
    try: 
        print('Reading RTOFS')
        RTOFS_oper_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,folder_fig)
    except Exception as err:
        print(err)    

    try:
        print('Reading GOFS 3.1')
        GOFS31_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,folder_fig)
    except Exception as err:
        print(err) 

    try:
        print('Reading Copernicus')
        Copernicus_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,folder_fig)
    except Exception as err:
        print(err) 
