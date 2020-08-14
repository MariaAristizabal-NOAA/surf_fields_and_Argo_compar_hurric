
"""
Created on Mon Jun 15 14:31:20 2020

@author: aristizabal
"""
    
#%% url NHC    

url_nhc = 'https://www.nhc.noaa.gov/gis/'

# GoMex
lon_lim_GoMex = [-100,-78]
lat_lim_GoMex = [15,32.5]

# Caribbean
#lon_lim_Car = [-87,-60]
#lat_lim_Car = [8,30]
lon_lim_Car = [-87,-60]
lat_lim_Car = [8,45]

# SAB
lon_lim_SAB = [-82,-64]
lat_lim_SAB = [25,36]

# MAB
lon_lim_MAB = [-77,-52]
lat_lim_MAB = [35,46]

# Atlantic
lon_lim_ATL = [-60,-10]
lat_lim_ATL = [8,45]

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

#%% Get time bounds for current day
ti = datetime.today()
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
zip_files_track_latest = [f for f in zip_files if 'TRACK' in f]

#%%
for i,f in enumerate(zip_files_track_latest):
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

        zip_file_cone_latest = [fl for fl in zip_files if np.logical_and(f.split('_')[0][2:] in fl,'CONE_latest' in fl)][0]
        kmz = ZipFile(zip_file_cone_latest, 'r')
        kml_ff = glob.glob(zip_file_cone_latest[:-4]+'/*.kml')
        kml_cone_latest = kmz.open(kml_ff[0].split('/')[1], 'r').read()
        soup = BeautifulSoup(kml_cone_latest,'html.parser')

        lon_forec_cone = []
        lat_forec_cone = []
        for i,s in enumerate(soup.find_all("coordinates")):
            coor = s.get_text('coordinates').split(',0')
            for st in coor[1:-1]:
                lon_forec_cone.append(st.split(',')[0])
                lat_forec_cone.append(st.split(',')[1])

        lon_forec_cone = np.asarray(lon_forec_cone).astype(float)
        lat_forec_cone = np.asarray(lat_forec_cone).astype(float)        

        zip_file_best_track = [fl for fl in zip_files if np.logical_and(f.split('_')[0][2:] in fl,'best_track' in fl)][0]
        kmz = ZipFile(zip_file_best_track, 'r')
        kml_ff = glob.glob(zip_file_best_track[:-4]+'/*.kml')
        kml_best_track = kmz.open(kml_ff[0].split('/')[1], 'r').read()
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
    if np.logical_or(np.min(lon_forec_track) >= lon_lim_GoMex[0],\
                  np.min(lon_forec_track) <= lon_lim_GoMex[1]):
        if np.logical_and(np.min(lat_forec_track) >= lat_lim_GoMex[0],\
                  np.min(lat_forec_track) <= lat_lim_GoMex[1]):
            lon_lim = lon_lim_GoMex
            lat_lim = lat_lim_GoMex
            temp_lim = [25,31.1]
            salt_lim = [31,37.1]
            temp200_lim = [5,24.6]
            salt200_lim = [35.5,37.6]
            tempb_lim = [0,25.6]
            tempt_lim = [6,30.6]
        
    if np.logical_and(np.min(lon_forec_track) >= lon_lim_Car[0],\
                  np.min(lon_forec_track) <= lon_lim_Car[1]):
        if np.logical_and(np.min(lat_forec_track) >= lat_lim_Car[0],\
                  np.min(lat_forec_track) <= lat_lim_Car[1]):
            lon_lim = lon_lim_Car
            lat_lim = lat_lim_Car  
            temp_lim = [25,31.1]
            salt_lim = [31,37.1]
            temp200_lim = [5,24.6]
            salt200_lim = [35.5,37.6]
            tempb_lim = [0,25.6]        
            tempt_lim = [6,30.6]

    if np.logical_and(np.min(lon_forec_track) >= lon_lim_SAB[0],\
                  np.min(lon_forec_track) <= lon_lim_SAB[1]):
        if np.logical_and(np.min(lat_forec_track) >= lat_lim_SAB[0],\
                  np.min(lat_forec_track) <= lat_lim_SAB[1]):
            lon_lim = lon_lim_SAB
            lat_lim = lat_lim_SAB 
            temp_lim = [25,30.6]
            salt_lim = [35,37.6]
            temp200_lim = [10,22.6]
            salt200_lim = [35.5,37.6]
            tempb_lim = [0,28.6]
            tempt_lim = [6,30.6]
        
    if np.logical_and(np.min(lon_forec_track) >= lon_lim_MAB[0],\
                  np.min(lon_forec_track) <= lon_lim_MAB[1]):
        if np.logical_and(np.min(lat_forec_track) >= lat_lim_MAB[0],\
                  np.min(lat_forec_track) <= lat_lim_MAB[1]):
            lon_lim = lon_lim_MAB
            lat_lim = lat_lim_MAB 
            temp_lim = [10,30.6]
            salt_lim = [30,36.6]
            temp200_lim = [5,20.6]
            salt200_lim = [35.0,37.6]
            tempb_lim = [0,26.6]
            tempt_lim = [6,28.6]

    if np.max(lon_forec_track) > lon_lim_MAB[1]:
        lon_lim = lon_lim_ATL
        lat_lim = lat_lim_ATL
        temp_lim = [22,31.1]
        salt_lim = [31,38.1]
        temp200_lim = [4,24.6]
        salt200_lim = [35.0,37.1]
        tempb_lim = [0,25.6]
        tempt_lim = [6,30.6]

    lon_lim = [np.min(lon_forec_track)-5,np.max(lon_forec_track)+5]
    lat_lim = [np.min(lat_forec_track)-5,np.max(lat_forec_track)+5]

#%% 
    if np.logical_and(len(name) != 0,len(lon_lim) != 0):
        os.chdir('/www/web/rucool/hurricane/Hurricane_season_' + str(tini.year))
        os.system('mkdir ' + ti.strftime('%b-%d')  ) 
        os.system('mkdir ' + ti.strftime('%b-%d') + '/' + name.split(' ')[1] )

        folder_fig = '/www/web/rucool/hurricane/Hurricane_season_' \
             + str(tini.year)+ '/' + ti.strftime('%b-%d') + '/' + name.split(' ')[1] + '/'

        try:
            print('Reading Argo floats')
            GOFS_RTOFS_vs_Argo_floats(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,folder_fig)
        except Exception as err:
            print(err)
   
        try: 
            print('Reading RTOFS')
            RTOFS_oper_baffin(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig)
        except Exception as err:
            print(err)    

        try:
            print('Reading GOFS 3.1')
            GOFS31_baffin(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig)
        except Exception as err:
            print(err) 

        try:
            print('Reading Copernicus')
            Copernicus_baffin(lon_forec_track,lat_forec_track,lon_forec_cone,lat_forec_cone,lon_best_track,lat_best_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig)
        except Exception as err:
            print(err) 

        os.chdir('/home/aristizabal/Code/Surf_fields_and_Argo_compar_hurric')

#%%
'''
name = '03l2020 Cristobal'

lon_forec_track = np.asarray([-90,-91.1,-90,-91.3,-91.1,-91.2,-91.5,-92.0,-92.5,-93,\
                   -93.3,-92.5,-92.5,-92.4,-92.1,-91.1,-90.0])
lat_forec_track = np.arange(16,33,1)

os.chdir('/www/web/rucool/hurricane/Hurricane_season_' + str(tini.year))
os.system('mkdir ' + ti.strftime('%b-%d')  )
os.system('mkdir ' + ti.strftime('%b-%d') + '/' + name.split(' ')[1] )
folder_fig = '/www/web/rucool/hurricane/Hurricane_season_' \
             + str(tini.year)+ '/' + ti.strftime('%b-%d') + '/' + name.split(' ')[1] + '/'

if np.logical_or(np.min(lon_forec_track) >= lon_lim_GoMex[0],\
              np.min(lon_forec_track) <= lon_lim_GoMex[1]):
    if np.logical_and(np.min(lat_forec_track) >= lat_lim_GoMex[0],\
              np.min(lat_forec_track) <= lat_lim_GoMex[1]):
        lon_lim = lon_lim_GoMex
        lat_lim = lat_lim_GoMex
        temp_lim = [25,31.6]
        salt_lim = [31,37.6]
        temp200_lim = [5,24.6]
        salt200_lim = [35.5,37.0]
        tempb_lim = [0,24.6]
        tempt_lim = [6,30.6]

if np.logical_and(np.min(lon_forec_track) >= lon_lim_Car[0],\
              np.min(lon_forec_track) <= lon_lim_Car[1]):
    if np.logical_and(np.min(lat_forec_track) >= lat_lim_Car[0],\
              np.min(lat_forec_track) <= lat_lim_Car[1]):
        lon_lim = lon_lim_Car
        lat_lim = lat_lim_Car
        temp_lim = [25,31.6]
        salt_lim = [31,37.6]
        temp200_lim = [5,24.6]
        salt200_lim = [35.5,37.0]
        tempb_lim = [0,24.6]
        tempt_lim = [6,30.6]
        
if np.logical_and(np.min(lon_forec_track) >= lon_lim_SAB[0],\
              np.min(lon_forec_track) <= lon_lim_SAB[1]):
    if np.logical_and(np.min(lat_forec_track) >= lat_lim_SAB[0],\
              np.min(lat_forec_track) <= lat_lim_SAB[1]):
        lon_lim = lon_lim_SAB
        lat_lim = lat_lim_SAB
        temp_lim = [10,30.6]
        salt_lim = [30,36.6]
        temp200_lim = [5,20.6]
        salt200_lim = [35.5,36.6]
        tempb_lim = [0,26.6]
        tempt_lim = [6,28.6]

if np.logical_and(np.min(lon_forec_track) >= lon_lim_MAB[0],\
              np.min(lon_forec_track) <= lon_lim_MAB[1]):
    if np.logical_and(np.min(lat_forec_track) >= lat_lim_MAB[0],\
              np.min(lat_forec_track) <= lat_lim_MAB[1]):
        lon_lim = lon_lim_MAB
        lat_lim = lat_lim_MAB
        temp_lim = [10,30.6]
        salt_lim = [30,36.6]
        temp200_lim = [5,20.6]
        salt200_lim = [35.5,36.6]
        tempb_lim = [0,26.6]
        tempt_lim = [6,28.6]

print('Reading RTOFS')
RTOFS_oper_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig)

print('Reading GOFS 3.1')
GOFS31_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig)

print('Reading Copernicus')
Copernicus_baffin(lon_forec_track,lat_forec_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig)

os.chdir('/home/aristizabal/Code/surf_fields_and_Argo_compar_hurric')
'''
