# -*- coding: utf-8 -*-
""" Script to download both ATL08 or ATL06 data
Created on Feb 11, 2020
Last Run/Update: Mar 30, 2022 using v4.1; Aug 26, 2021

Modification2: origianl name: icesat2_search_and_download_ATL
TODO: on calling read_atl06 function, following exception
    ERROR:root:	Exception in reading hdf group (ground track), df length = 17591
    ERROR:root:	Exception in reading hdf group (ground track), df length = 17591
    
@author: Bidhya N Yadav
"""
#%%
from doctest import master
import math
import os, sys
import shutil
import time
from datetime import timedelta

import numpy as np
import pandas as pd

import geopandas as gpd
from shapely.geometry import Point #, Polygon, mapping
# from shapely.geometry.polygon import orient
from statistics import mean
import h5py
# # To read KML files with geopandas, we will need to enable KML support in fiona (disabled by default)
# fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
import json
import zipfile
import io

import requests

from xml.etree import ElementTree as ET
import logging
# from rema_dem import get_rema_strip_polygon
from create_dem_metadata import get_strip_polygon

import warnings  # March 30, 2022  
# warnings.simplefilter  (action='ignore', category=FutureWarning)  # FutureWarning: pandas.Int64Index is deprecated and will be removed from pandas in a future version. Use pandas.Index with the appropriate dtype instead.
# warnings.simplefilter  (action='ignore', category=UserWarning)    # UserWarning: Column names longer than 10 characters will be truncated when saved to ESRI Shapefile.
warnings.filterwarnings(action="ignore", category=FutureWarning, message="pandas.Int64Index")  ##FutureWarning: pandas.Int64Index is deprecated and will be removed from pandas in a future version
# warnings.filterwarnings(action='ignore', category=UserWarning, message="Column names longer than 10 characters will be truncated") # UserWarning: Column names longer than 10 characters will be truncated when saved to ESRI Shapefile.
# warnings.filterwarnings(action='ignore', category=UserWarning, message="Geometry is in a geographic CRS")
warnings.filterwarnings(action='ignore', category=DeprecationWarning)  # , message="The array interface is deprecated and will no longer work in Shapely 2.0"
warnings.filterwarnings(action="ignore", category=UserWarning)
warnings.filterwarnings(action="ignore", category=Warning)

# One more: ShapelyDeprecationWarning: The array interface is deprecated and will no longer work in Shapely 2.0. Convert the '.coords' to a numpy array instead.

# Some global variables and definition [must choose one]
short_name = 'ATL06'  # 'ATL08' 'ATL06' 'ATL03'
# ===================================================================================


def move_files_from_order(icesat2_path):
    ''' Extract files from downloaded subfoder [ie one by orderid] 
        and move all hdf files in one location
    '''
    hdf_path = icesat2_path
    for root, dirs, files in os.walk(f'{icesat2_path}/downloads', topdown=False):
        for f in files:
            if f.endswith('h5') or f.endswith('.hdf'):  # .hdf was added for modis data
                try:
                    shutil.move(os.path.join(root, f), hdf_path)
                except:
                    logging.error("Extraction from downloads Error (SHUTIL): perhaps file already exist")
    # os.rmdir(f'{icesat2_path}/downloads') #for empty
    shutil.rmtree(f'{icesat2_path}/downloads')

#%% To Parse hdf file and convert csv and shapefile for analysis/visualization
from astropy.time import Time
def gps2dyr(time, offset = 0):
    """ Converte GPS time to decimal years. Helper function"""
    time = time + offset
    gps_time = Time(time, format='gps') #.decimalyear
    iso_time = Time(gps_time, format='iso')
    iso_time = iso_time.value
    # Conver to pandas datetime [not sure if it is utc]
    dt = pd.to_datetime(iso_time)
    return dt


# Convert points to line (for visualization and maybe even for cross-track analysis)
def convert_points_to_line(icesat2_path, short_name = 'ATL06'):
    ''' Conbine all the csv files of ATL03 ground tracks into one shapefile of polylines
        This is initially conceptualized for easy visualization
        In future, this can be used for cross-track analysis
    '''
    files = os.listdir(icesat2_path)
    # csv_files = [f for f in files if f.endswith('.csv') and 'ATL03' in f]
    csv_files = [f for f in files if f.endswith('.csv') and short_name in f]

    res_dict = {}
    for idx, fname in enumerate(csv_files):
        df  = pd.read_csv(f'{icesat2_path}/{fname}')
        line_dict = {}
        for strip in df.strip.unique():
            gt = df[df.strip==strip]
            if len(gt)>=2:#else line cannot be constructed
                geom = LineString(gt[['lon', 'lat', 'h_li']].values)
                line_dict[strip] = geom
        line_df = pd.DataFrame.from_dict(line_dict, orient='index', columns=['geometry'])
        line_df['fname'] = fname.split('.csv')[0] #to drop csv suffix
        res_dict[idx] = line_df
    for k in res_dict.keys():
        if k == 0:
            dfx = res_dict[k]
            dfx['idx'] = k
        else:
            dfx1 = res_dict[k]
            dfx1['idx'] = k
            dfx = pd.concat([dfx, dfx1], axis=0)
    dfx['strip'] = dfx.index
    dfx['date'] = dfx.fname.apply(lambda x: pd.to_datetime(x.split('_')[2])) #[:8]
    dfx['rgt'] = dfx.fname.apply(lambda x: x.split('_')[3][:4])
    dfx['cycle'] = dfx.fname.apply(lambda x: x.split('_')[3][4:6])
    dfx['segment'] = dfx.fname.apply(lambda x: x.split('_')[3][6:])


    line_gdf = gpd.GeoDataFrame(dfx, geometry='geometry', crs = 'EPSG:4326')
    proj4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' # TODO : check how this is used
    line_gdf = line_gdf.to_crs(proj4)  # so we can calculate length
    line_gdf['date'] = line_gdf['date'].dt.strftime('%Y-%m-%d')  # To prevent DriverSupportError: ESRI Shapefile does not support datetime fields
    # line_gdf.to_file(f'{icesat2_path}/lines.gpkg', driver='GPKG')
    return line_gdf


""" 
    files = os.listdir(icesat2_path)
    hdf_files = [f for f in files if f.endswith('.h5') and 'ATL06' in f]
    logging.info(f'Number of HDF files : {len(hdf_files)}')
    for f in hdf_files:
        logging.info(f'  Parsing : {f}')
        hdf_path = f'{icesat2_path}/{f}'

"""


def read_atl03(icesat2_path, append_attrs=True, output_shapefile=False, output_folder=None, gis_output='shp', overwrite=False):
    """ Read ATL03 file and output 6 reduced files.
        Slightly different from ATL06 and 08 due to data volume
        We will save each of six beams in a separate shapefile
        hdf_path: Full path to ATL03 hdf file [this is different from ATL06 08
        output_path : If you want to save outputs (csv, shp) to different folder 

        NB: Here, we don't have fill values for height, so perhaps not required to put nan then drop rows!
            Watch out of ERROR is size of Shapefile > 2GB
    """
    files = os.listdir(icesat2_path)
    hdf_files = [f for f in files if f.endswith('.h5') and 'ATL03' in f]
    logging.info(f'Number of HDF files : {len(hdf_files)}')
    for f in hdf_files:
        logging.info(f'  Parsing : {f}')
        hdf_path = f'{icesat2_path}/{f}'
        # This is the old default behaviour: extract the original hdf filename without the extension
        out_path = f'{os.path.splitext(hdf_path)[0]}'  # fullpath to hdf filename, but without the .h5 extension
        # logging.info(f"out_path : {out_path}")
        # if output_folder:
        #     # If output folder is specified for csv and shapefile, extract just the filename (without .h5) and append to the specified path
        #     hdf_fname = out_path.split('/')[-1]  # extract hdf file's filename only
        #     out_path = f'{output_folder}/{hdf_fname}'
        # print(out_path)
        # if os.path.exists(f'{out_path}.csv') and overwrite == False:
        #     print('Already Processed')
        #     return
        #     # elif overwrite == True:
        #     #     print('Overwriting csv outputs ...')
        #     #     pass  # dummy pass to prevent nesting below  [WHY?]
        logging.info(f'Processing HDF file: {hdf_path}')
        res_dict = {}  #not used because each ground-track will be saved separately
        # meta_dict = {}  # These will hold metadata required for scalars per ground-track
        # qual_str_count = ''
        group = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
        # Loop trough beams
        # Perhaps read file first, then loop through groups; should be faster
        with h5py.File(hdf_path, 'r') as fi:
            # subset group based on data
            group = [g for g in list(fi.keys()) if g in group]
            group1 = len(group)
            group = [g for g in group if 'heights' in fi[f'/{g}']]
            group2 = len(group)
            if len(group) == 0:
                logging.info('No heights information in any group')
                # return  # can't use continue here as no loop yet
                continue  # now that it is inside for loop perhaps continue should work
            if group2 < group1:
                logging.info(f'Some Non-empty groups: {group2}/{group1}')
            # NB: Assert if at least one group present else may be error due to enumeration
            for k,g in enumerate(group):
                # if 'heights' in fi[f'/{g}']:
                # 1) Read in data for a single beam
                lat = fi[f'/{g}/heights/lat_ph']
                lon = fi[f'/{g}/heights/lon_ph']
                h_ph = fi[f'/{g}/heights/h_ph']
                
                dist_ph_along = fi[f'/{g}/heights/dist_ph_along']
                quality_ph = fi[f'/{g}/heights/quality_ph']  #new march 31, 2022 Photon quality
                # Indicates the quality of the associated photon. 0=nominal, 1=possible_afterpulse, 2=possible_impulse_response_effect, 3=possible_tep. Use this flag in conjunction with signal_conf_ph to identify those photons that are likely noise or likely signal
                t_dt = fi[f'/{g}/heights/delta_time']
                t_ref = fi['/ancillary_data/atlas_sdp_gps_epoch'] #scalar 1 value; required for offset
                segment_length = fi[f'/{g}/geolocation/segment_length']
                segment_ph_cnt = fi[f'/{g}/geolocation/segment_ph_cnt']

                # Confidence level associated with each photon event selected as signal.
                # 0=noise. 1=added to allow for buffer but algorithm classifies as background; 2=low; 3=med; 4=high
                # This parameter is a 5xN array where N is the number of photons in the granule, and the 5 rows indicate signal finding for each surface type (in order: land, ocean, sea ice, land ice and inland water).
                # Events not associated with a specific surface type have a confidence level of -1. Events evaluated as TEP returns have a confidence level of -2
                signal_conf_ph = fi[f'/{g}/heights/signal_conf_ph'] # q_flag (old name) This is N by 5 matrix so can't be converted to dataframe directly            
                # 2) Make Pandas dataframe
                # Error: , 'segment_length':segment_length, 'segment_ph_cnt':segment_ph_cnt
                df = pd.DataFrame({'lon':lon, 'lat':lat, 't_dt':t_dt, 'h_ph': h_ph, 'quality_ph':quality_ph, 'dist_ph_along':dist_ph_along})
                # df = pd.DataFrame({'lon':lon, 'lat':lat, 't_dt':t_dt, 'h_ph': h_ph, 'quality_ph':quality_ph, 'dist_ph_along':dist_ph_along, 'segment_length':segment_length, 'segment_ph_cnt':segment_ph_cnt})
                # , 'signal_conf_ph':signal_conf_ph 'q_flag':q_flag, 't_dt':t_dt
                # Convert GPS time to actual time using function
                df['t_dt'] = df['t_dt'].apply(gps2dyr, offset=t_ref[0])

                if append_attrs:
                    # Append the Flag values for photon classification (this is N by 5 matrix)
                    # TODO: investigate this more
                    flag_cols = ['land', 'ocean', 'sea_ice', 'land_ice', 'inld_water']  # inland_water is more than 10 chars hence problem with shapefile
                    # flag_df = pd.DataFrame(q_flag, columns=flag_cols)
                    flag_df = pd.DataFrame(signal_conf_ph, columns=flag_cols)
                    df = pd.concat([df, flag_df], axis=1)
                # To save separately by ground-track
                df.to_csv(f'{out_path}_{g}.csv', index=False)
                # Collect Dataframe in a Dict to combine all ground tracks for 1 h5 file in single csv
                res_dict[g] = df
            # Combine Dataframes for each of 6 ground-tracks into single Dataframe
            count = 0
            for k in res_dict.keys():
                if count == 0:
                    df = res_dict[k]
                    df['strip'] = k
                    count += 1
                else:
                    df1 = res_dict[k]
                    df1['strip'] = k
                    df = pd.concat([df, df1], axis=0)
            # Choose filename for csv and shapefile
            df.to_csv(f'{out_path}.csv', index=False)

            # if output_shapefile:
            #     # 2. Convert to Geopandas
            #     # Create shapefile for a subset (1 or 0.1% of all points to get the overview of points
            #     #step = int(len(df)/10000)
            #     #if step < 1:
            #     #    # to guard against very small datasets
            #     #    step = 1
            #     #df = df.iloc[::step]
            #     df['geometry'] = df[['lon', 'lat']].apply(lambda x: Point(x), axis=1) #because this takes a long time
            #     gdf = gpd.GeoDataFrame(df.drop(columns=['lon', 'lat']), geometry='geometry', crs = 'epsg:4326')
            #     # gdf = gpd.GeoDataFrame(df[['h_li', 'geometry']], geometry='geometry', crs = {'init': 'epsg:4326'})
            #     # gdf.crs = {'init': 'epsg:4326'} #not yet verified or checked with what ICESAT-2 metadata provides
            #     if gis_output == 'shp':
            #         gdf.to_file(f'{out_path}.shp')
            #     elif gis_output == 'gpkg':
            #         gdf.to_file(f'{out_path}.gpkg', driver='GPKG')


def read_atl06(icesat2_path, gis_output='shp'):
    """ Read 1 ATL06 file and output 6 reduced files.     
        Extract variables of interest and separate the ATL06 file 
        into each beam (ground track) and ascending/descending orbits."""
    files = os.listdir(icesat2_path)
    hdf_files = [f for f in files if f.endswith('.h5') and 'ATL06' in f]
    logging.info(f'Number of HDF files : {len(hdf_files)}')
    for f in hdf_files:
        logging.info(f'  Parsing : {f}')
        hdf_path = f'{icesat2_path}/{f}'
        res_dict = {}
        meta_dict = {} #These will hold metadata required for scalars per ground-track
        group = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
        qual_str_count = ''
        # Loop trough beams
        # Perhaps read file first, then loop through groups; should be faster
        with h5py.File(hdf_path, 'r') as fi:
            # subset group based on data
            group = [g for g in list(fi.keys()) if g in group]
            group1 = len(group)
            group = [g for g in group if 'land_ice_segments' in fi[f'/{g}']]
            group2 = len(group)
            if group2<group1:
                logging.info(f'\tNon-empty groups: {group2}/{group1}')
            # NB: Assert if at least one group present else may be error due to enumeration
            if len(group) == 0:
                logging.info('No any Group!')
                continue
            for k,g in enumerate(group):
                #if 'land_ice_segments' in fi[f'/{g}']:
                # 1) Read in data for a single beam #
                lat = fi[f'/{g}/land_ice_segments/latitude'][:]
                lon = fi[f'/{g}/land_ice_segments/longitude'][:]
                h_li = fi[f'/{g}/land_ice_segments/h_li'][:] #nan
                #s_li = fi[f'/{g}/land_ice_segments/h_li_sigma'][:] #nan
                t_dt = fi[f'/{g}/land_ice_segments/delta_time'][:]
                q_flag = fi[f'/{g}/land_ice_segments/atl06_quality_summary'][:]
                t_ref = fi['/ancillary_data/atlas_sdp_gps_epoch'][:] #scalar 1 value; required for offset
                #rgt = fi['/orbit_info/rgt'][:] * np.ones(len(lat)) #scalar 1 value
                #orb = np.full_like(h_li, k) #scalar 1 value
                
                meta_dict['t_ref'] = t_ref #dictionary of metadata (will be used in future to investigate data)
                # 2) Make Pandas dataframe
                #t_dt = gps2dyr(t_dt, offset=t_ref) #time consuming 
                df = pd.DataFrame({'lon':lon, 'lat':lat, 'h_li': h_li, 'q_flag':q_flag, 't_dt':t_dt})
                #Convert GPS time to actual time using function
                df['t_dt'] = df['t_dt'].apply(gps2dyr, offset=t_ref[0])
                #df['fname'] = f.split('.')[0] #to save hdf filename; could important when multiple hdf files over a rema tile
                #df.index = df.t_dt #creates difficultly in plotting
                # Fill Nans for na-data and drop
                df.loc[df.h_li>3e38, 'h_li'] = np.nan
                df = df.dropna()
                all_points = len(df)
                #df = df[df.q_flag==0] #select only flags of zero (good values) can be empty sometimes
                good_quality_points = len(df) #len(df[df.q_flag==0])                
                qual_str_count = qual_str_count + f'{g}={good_quality_points}/{all_points}; ' #Appending to save string for all groups within a granule
                if len(df)>0:
                    # Assemble ground track into a dictionary, later we convert to csv and shp through df
                    res_dict[g] = df
            #----------------------------------------------------------------------------------------------
            # Now that ATL06 data from separate ground tracks are in one dict, merge it to df and save to csv/shp
            if len(res_dict)>0:
                if good_quality_points < all_points:
                    logging.info(f'\t\tGood/Total Points: {qual_str_count}   Total GTs = {len(res_dict)}')
                # To guard againt empty result dictionary created with no icesat2 passing the quality control above
                # 1. Combine Dataframes for each of 6 ground-tracks into single Dataframe
                count = 0
                for k in res_dict.keys():
                    # k = 'gt1l', 'gt1r' etc
                    if count == 0:
                        df = res_dict[k]
                        df['strip'] = k
                        count += 1
                    else:
                        df1 = res_dict[k]
                        df1['strip'] = k
                        df = pd.concat([df, df1], axis=0)
                # Log the time range of icesat2 data (could be useful for understanding why some data is large)
                # This may be creating exception when empty
                time_range = df.t_dt.max() - df.t_dt.min()
                # logging.info(f"\tRows = {len(df)} \t Time Range = {time_range.total_seconds()} seconds")
                # Choose filename for csv and shapefile
                atl_fname = os.path.splitext(hdf_path)[0].split('/')[-1]
                #df = df[df.q_flag==0] # Already done above for each ground track
                df.to_csv(f'{icesat2_path}/{atl_fname}.csv', index=False)
                
                # 2. Convert to Geopandas
                df['geometry'] = df[['lon', 'lat']].apply(lambda x: Point(x), axis=1)
                gdf = gpd.GeoDataFrame(df[['t_dt', 'h_li', 'q_flag', 'strip', 'geometry']], geometry='geometry', crs = 'EPSG:4326')
                #gdf.crs = {'init': 'epsg:4326'} # Verify with ICESAT-2 metadata
                if gis_output=='shp':
                    gdf['t_dt'] = gdf['t_dt'].dt.strftime('%Y-%m-%d %H:%M:%S.%f') #To prevent DriverSupportError: ESRI Shapefile does not support datetime fields
                    gdf.to_file(f'{icesat2_path}/{atl_fname}.shp')
                elif gis_output=='gpkg':
                    gdf.to_file(f'{icesat2_path}/{atl_fname}.gpkg', driver='GPKG')
            else:
                logging.info(f"\t\tNo good quality Ground Track in this HDF file; csv or shp not created")
    # Perhaps merge if more than one hdf file present per tile (or write 2 separate functions)
#----------------------------------------------------------------------------------------------------

def read_atl08(icesat2_path, gis_output='shp'):
    """ Read 1 ATL06 file and output 6 reduced files.     
        Extract variables of interest and separate the ATL06 file 
        into each beam (ground track) and ascending/descending orbits.
    """
    files = os.listdir(icesat2_path)
    hdf_files = [f for f in files if f.endswith('.h5') and 'ATL08' in f]
    logging.info(f'Number of HDF files: {len(hdf_files)}')
    for f in hdf_files:
        logging.info(f'Processing HDF file: {f}')
        hdf_path = f'{icesat2_path}/{f}'
        res_dict = {}
        meta_dict = {} #These will hold metadata required for scalars per ground-track
        group = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
        qual_str_count = ''
        # Loop trough beams
        # Perhaps read file first, then loop through groups; should be faster
        with h5py.File(hdf_path, 'r') as fi:
            # subset group based on data
            group = [g for g in list(fi.keys()) if g in group]
            group1 = len(group)
            group = [g for g in group if 'land_segments' in fi[f'/{g}']]
            group2 = len(group)
            if group2<group1:
                logging.info(f'Non-empty groups: {group2}/{group1}')
            # NB: Assert if at least one group present else may be error due to enumeration
            if len(group) == 0:
                logging.info('No any Group!')
                continue

            t_ref = fi['/ancillary_data/atlas_sdp_gps_epoch'][:] #scalar 1 value
            #surf_type = fi['/ds_surf_type'][:] #land water inland water etc; but probably does not exist, giving error
            for k,g in enumerate(group):
                #try:
                # 1) Read in data for a single beam #
                lat = fi[f'/{g}/land_segments/latitude'][:]
                lon = fi[f'/{g}/land_segments/longitude'][:]
                #h_li = fi[f'/{g}/land_segments/terrain/h_te_median'][:] #nan
                #s_li = fi[f'/{g}/land_ice_segments/h_li_sigma'][:] #nan
                t_dt = fi[f'/{g}/land_segments/delta_time'][:]  #mean_pass_time: Mean time for the segment in number of GPS seconds since the ATLAS SDP epoch. 
                # The ATLAS Standard Data Products (SDP) epoch offset is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and 
                # the ATLAS SDP epoch. By adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the time in gps_seconds relative to the GPS epoch can be computed.";

                layer_flag = fi[f'/{g}/land_segments/layer_flag'][:]  #Not required perhaps; likely_clear likely_cloudy; 1 means clouds or blowing snow are likely present. A value of 0 indicates the likely absence of clouds or blowing snow
                # Extract everything for terrain and canopy variables
                terrain_keys = fi[f'{g}/land_segments/terrain'].keys()
                terrain_keys = list(terrain_keys)
                terrain_keys.remove('h_te_best_fit_20m') #2D data, can't be saved in table/shp
                terrain_keys.remove('subset_te_flag') #2D data, can't be saved in table/shp

                terrain_dict = {}
                for tk in terrain_keys:
                    terrain_dict[tk] = fi[f'{g}/land_segments/terrain/{tk}'] #no data yet, we are just mapping
                #terrain = pd.DataFrame.from_dict(terrain_dict)
                # Do the same for Canopy
                canopy_keys = fi[f'{g}/land_segments/canopy'].keys()
                canopy_keys = list(canopy_keys)
                # Remove these two keys as they contain some multivalue tuples which are just no-data in version 002; updated for version 5 [Feb 01, 2022]
                canopy_keys.remove('canopy_h_metrics')
                canopy_keys.remove('canopy_h_metrics_abs')
                canopy_keys.remove('h_canopy_20m')
                canopy_keys.remove('subset_can_flag')
                canopy_dict = {}
                for ck in canopy_keys:
                    canopy_dict[ck] = fi[f'{g}/land_segments/canopy/{ck}']
                    
                # Remove the canopy subset flag that seems added in version 3 because this is an array and we can't save array to shapefile
                # Not required anymore as these are removed from list above [Feb 01, 2022]
                # if 'subset_can_flag' in canopy_dict:
                #     del canopy_dict['subset_can_flag']
                # if 'subset_te_flag' in terrain_dict:
                #     del terrain_dict['subset_te_flag']

                #Merge two dictionaries (order should be retained; verify when running again with Python 3.7)
                #https://thispointer.com/how-to-merge-two-or-more-dictionaries-in-python/
                terrain_dict.update(canopy_dict)

                #rgt = fi['/orbit_info/rgt'][:] * np.ones(len(lat)) #scalar 1 value
                #orb = np.full_like(h_li, k) #scalar 1 value
                
                meta_dict['t_ref'] = t_ref #dictionary of metadata (will be used in future to investigate data)
                
                # 2) To Make Pandas dataframe
                #t_dt = gps2dyr(t_dt, offset=t_ref) #time consuming 
                # Collect everythin into one dictionary
                gt_dict = {'lon':lon, 'lat':lat, 't_dt':t_dt, 'layer_flag':layer_flag}
                gt_dict.update(terrain_dict)
                #df = pd.DataFrame({'lon':lon, 'lat':lat, 'h_li': h_li, 'q_flag':q_flag, 't_dt':t_dt})
                df = pd.DataFrame.from_dict(gt_dict)
                # Jan 31, 2021: attrs Fillvalues does not exist in version, or may have been moved to new location
                nan_value =  fi[f'/{g}/land_segments/terrain/h_te_mean'].attrs['_FillValue']  #prior we needed [0] to extract value; all terrain keys have same fill values, so can use one
                nan_value = nan_value.item(0)  #Feb 02, 2022: now this is array, extract the value
                #nan_value = df.max().max() #Not Guaranteed to work for small segments, so don't rely on this
                # df= df.replace(nan_value, np.nan)  ##error
                df = df[df.apply(lambda x: x!=nan_value)]  #new Jan 31, 2021; but not verified

                #Convert GPS time to actual time using function
                df['t_dt'] = df['t_dt'].apply(gps2dyr, offset=t_ref[0])
                #df['fname'] = f.split('.')[0] #to save hdf filename; could important when multiple hdf files over a rema tile
                #df.index = df.t_dt #creates difficultly in plotting
                # Fill Nans for na-data and drop
                #df.loc[df.h_li>3e38, 'h_li'] = np.nan     OLD               
                #df = df.dropna() OLD

                all_points = len(df)
                #df = df[df.q_flag==0] #select only flags of zero (good values) can be empty sometimes
                good_quality_points = len(df) #len(df[df.q_flag==0])
                qual_str_count = qual_str_count + f'{g}={good_quality_points}/{all_points}; '
                if len(df)>0:
                    # Assemble ground track into a dictionary, later we convert to csv and shp through df
                    res_dict[g] = df
                # except:
                #     # Most like this error is due to empty dataframe 
                #     # may not exist anymore since we are no dropping the bad quality data
                #     logging.error(f'\tError parsing {g}')
            #----------------------------------------------------------------------------------------------
            # Now that ATL08 data from separate ground tracks are in one dict, merge it to df and save to csv/shp
            if len(res_dict)>0:
                logging.info(f'\tGood/Total Points: {qual_str_count}   Total GTs = {len(res_dict)}')
                # To guard againt empty result dictionary created with no icesat2 passing the quality control above
                # 1. Combine Dataframes for each of 6 ground-tracks into single Dataframe
                count = 0
                for k in res_dict.keys():
                    # k = 'gt1l', 'gt1r' etc
                    if count == 0:
                        df = res_dict[k]
                        df['strip'] = k
                        count += 1
                    else:
                        df1 = res_dict[k]
                        df1['strip'] = k
                        df = pd.concat([df, df1], axis=0)
                # Log the time range of icesat2 data (could be useful for understanding why some data is large)
                # This may be creating exception when empty
                time_range = df.t_dt.max() - df.t_dt.min()
                #logging.info(f"\tRows = {len(df)} \t Time Range = {time_range.total_seconds()} seconds")
                # Choose filename for csv and shapefile
                atl_fname = os.path.splitext(hdf_path)[0].split('/')[-1]
                #df = df[df.q_flag==0] # Already done above for each ground track
                df.to_csv(f'{icesat2_path}/{atl_fname}.csv', index=False)
                
                # 2. Convert to Geopandas
                df['geometry'] = df[['lon', 'lat']].apply(lambda x: Point(x), axis=1)
                #gdf = gpd.GeoDataFrame(df[['t_dt', 'h_li', 'q_flag', 'strip', 'geometry']], geometry='geometry')
                gdf = gpd.GeoDataFrame(df, geometry='geometry', crs = 'EPSG:4326')
                # gdf.crs = {'init': 'epsg:4326'} #not yet verified or checked with what ICESAT-2 metadata provides
                if gis_output=='shp':
                    gdf['t_dt'] = gdf['t_dt'].dt.strftime('%Y-%m-%d %H:%M:%S.%f') #To prevent DriverSupportError: ESRI Shapefile does not support datetime fields
                    gdf.to_file(f'{icesat2_path}/{atl_fname}.shp')
                elif gis_output=='gpkg':
                    gdf.to_file(f'{icesat2_path}/{atl_fname}.gpkg', driver='GPKG')
                # gdf_json = gdf.to_json()
                # with open(f'{icesat2_path}/{atl_fname}.json', 'w') as json_file:
                #     json.dump(gdf_json, json_file)
            else:
                logging.info(f"\tNo good quality Ground Track in this HDF file; csv or shp not created")

# Read the shapefile: Either directly of from Rema Tile
def download_process_icesat2(strip_folder, strip, icesat2_path):
    ''' Main code to ownload IS-2 data by DEM strip
        Inputs:
        strip_folder: 
    icesat2_path: where files will be downloaded and processed; created inside this function
    Update Feb 11, 2020: Now we will just order the files with this function
    '''
    # This part is already taken care in the main script; this gives error anyway if the folder does not exist
    # # don't download if hdf (one or more already exist) here
    # files = os.listdir(icesat2_path)
    # hdf_files = [f for f in files if f.endswith('.h5')]
    # if len(hdf_files) > 0:
    #     logging.info(f"HDF files already exist, hdf_count = {len(hdf_files)}, skipping download")
    #     return
    
    # Get ouline of strip
    gdf = get_strip_polygon(strip_folder) #, OLD: get_rema_strip_polygon 'arctic' option used only for detour of hma
    gdf = gdf.to_crs('EPSG:4326')
    lat = gdf.centroid.y.values[0]  # to adapt timedelta based on latitude
    logging.info(f'Latitude = {lat:.2f}')
    shp_json = gdf.to_json()  # use json as shapefile, or kml too, shp file not yet tried

    # Create directories
    if not os.path.exists(icesat2_path):
        logging.info(f'Output folder: {icesat2_path}')
        os.makedirs(icesat2_path) #exist_ok=True to prevent complaint if directory exist
        os.mkdir(f'{icesat2_path}/downloads')
        #os.mkdir(f'{icesat2_path}/hdf')
        #os.mkdir(f'{icesat2_path}/gis')
    
    #Input temporal range     
    # Read datetime from shapefile : to extract the minimum and maximum dates for image acquisition
    min1 = pd.to_datetime(gdf.time1).min()
    min2 = pd.to_datetime(gdf.time2).min()
    max1 = pd.to_datetime(gdf.time1).max()
    max2 = pd.to_datetime(gdf.time2).max()
    tmin = min(min1, min2)
    tmax = max(max1, max2)
    img_time_separation = tmax - tmin
    # New (Aug 21, 2021): making time-delta a function of latitude
    if abs(lat) <= 65:
        days_offset = 30  # end_date = '2020-08-30'     # 3 months, farthest from pole TODO generate end month by adding number of months/days
    elif abs(lat) > 65 and abs(lat) <= 75:
        days_offset = 20
    elif abs(lat) > 75 and abs(lat) <= 80:
        days_offset = 10  # Closest to pole, select smaller end time
    elif abs(lat) > 80 and abs(lat) <= 85:
        days_offset = 5  # Closest to pole, select smaller end time
    elif abs(lat) > 85:
        days_offset = 3  # Closest to pole, select smaller end time    
    days_offset = timedelta(days = days_offset) # d originally 30 days
    tmin = tmin - days_offset #d
    tmax = tmax + days_offset #d
    
    #make sure dem is collected after icesat2 data began, else no use of further processing
    cutoff_time = pd.to_datetime('2018/10/15') #official '2018-10-14T00:00:00.000Z'
    if (tmax < cutoff_time):
        logging.info(f"NOT Processed: older than 2018/10/15")
        return
    # Saving shapefile is not required to download data; only for future ananlysis
    gdf.to_file(f'{icesat2_path}/strip_outline.shp')
    
    '''
    #% Use geojson, KML, or shapefile for POST request
    fiona.supported_drivers['KML'] = 'rw'
    kml_filepath = 'rema3.kml' #antar_test.kml
    # Save to kml to use with pre-existing workflow from summer school
    gdf.to_file(kml_filepath, driver='KML')
    #kml_filepath = f'{strip_folder}/outline/outline.kml'
    '''
    
    start_date, start_time = str(tmin).split()
    end_date, end_time = str(tmax).split()
    # Comment these two lines if not running for hma
    #start_date = '2018-10-15' #yyyy-MM-dd format '2018-10-14' #'2019-04-01' # start of ATL06 2018/10/14
    #end_date = '2019-12-01' #'2019-09-30'
    
    #start_date = '2019-04-01' #yyyy-MM-dd format '2018-10-14' #'2019-04-01' # start of ATL06 2018/10/14
    start_time = '00:00:00' #HH:mm:ss format
    #end_date = '2019-10-28' #'2019-09-30'
    end_time = '23:59:59'
    temporal = start_date + 'T' + start_time + 'Z' + ',' + end_date + 'T' + end_time + 'Z'
    logging.info(f'Image time separation = {img_time_separation}.  Days_offset added = {days_offset.days}. Total Days = {(tmax - tmin).days} Search Date: {start_date} to {end_date}.')
    
    # Get json response from CMR collection metadata . This is high-level metadata on a dataset or "collection".
    # In JSON format
    search_params = {'short_name': short_name}
    cmr_collections_url = 'https://cmr.earthdata.nasa.gov/search/collections.json'
    response = requests.get(cmr_collections_url, params=search_params)
    results = json.loads(response.content)
    #pprint.pprint(results)
    
    # Get the lastest version of available data
    versions = [i['version_id'] for i in results['feed']['entry']]
    latest_version = max(versions)
    # To use Bounding Box for subsetting
    bounding_box = ','.join(map(str, list(gdf.total_bounds))) #map will convert all float to str; then we do the join on string
    #aoi = '1'
    # Bounding box subsetting (bbox) in same format as bounding_box
    #bbox = bounding_box #not used in post  
    # bounding box input:
    search_params = {
            'short_name': short_name,
            'version': latest_version,
            'temporal': temporal,
            'page_size': 100,
            'page_num': 1,
            'bounding_box': bounding_box
            }
    
    # To use polygon for subsetting
    '''
    #Integer position based indexing of GeoDataFrame object to get it into a shapeply geometry object.
    poly = gdf.iloc[0].geometry
    # Simplify polygon. The larger the tolerance value, the more simplified the polygon.
    poly = poly.simplify(0.01, preserve_topology=False)
    # Orient counter-clockwise
    poly = orient(poly, sign=1.0)
    print(poly)
    #Format dictionary to polygon coordinate pairs for CMR polygon filtering
    polygon = ','.join([str(c) for xy in zip(*poly.exterior.coords.xy) for c in xy])
    
    # If polygon input (either via coordinate pairs or shapefile/KML/KMZ):
    search_params = {
            'short_name': short_name,
            'version': latest_version,
            'temporal': temporal,
            'page_size': 100,
            'page_num': 1,
            'polygon': polygon
            }
    # NB: bbox new parameter for spatial subsetting
    '''
    #print('CMR search parameters: ', search_params)    
    # Query number of granules using our (paging over results)
    granule_search_url = 'https://cmr.earthdata.nasa.gov/search/granules'
    headers={'Accept': 'application/json'} #if not defined in first block to get the token
    granules = []
    while True:
        response = requests.get(granule_search_url, params=search_params, headers=headers)
        results = json.loads(response.content)
        if len(results['feed']['entry']) == 0:
            # Out of results, so break out of loop
            break
        # Collect results and increment page_num
        granules.extend(results['feed']['entry'])
        search_params['page_num'] += 1
        #print(count)
        #count +=1
    # Get number of granules over my area and time of interest
    if len(granules)==0:
        logging.error("Stopping because Number of Granules = {len(granules)}")
        return    
    granule_sizes = [float(granule['granule_size']) for granule in granules]
    # Average size of granules in MB
    logging.info(f'Granules: Count = {len(granules)}, Avg Size = {mean(granule_sizes):.1f} MB, Total Size = {sum(granule_sizes):.2f} MB')
    
    # =============================================================================
    # Create a token
    '''
    We will generate a token needed in order to access data using your Earthdata Login credentials, and we will apply that token to the following queries. If you do not already have an Earthdata Login account, 
    go to http://urs.earthdata.nasa.gov to register. Your password will be prompted for privacy.
    
    # Earthdata Login credentials
    # Enter your Earthdata Login user name
    uid = 'bssxxya'
    # Enter your email address associated with your Earthdata Login account
    email = ''
    #pswd = getpass.getpass('Earthdata Login password: ')
    
    # Request token from Common Metadata Repository using Earthdata credentials
    token_api_url = 'https://cmr.earthdata.nasa.gov/legacy-services/rest/tokens'
    # hostname and ip may be not used!
    hostname = socket.gethostname() # get my personal computer name, eg staff-xxx-xxx; 'STAFF-BY-xx'
    ip = socket.gethostbyname(hostname)
    data = {
        'token': {
            'username': uid,
            'password': pswd,
            'client_id': 'NSIDC_client_id',
            'user_ip_address': ip
        }
    }    
    headers={'Accept': 'application/json'}
    #response = requests.post(token_api_url, json=data, headers=headers)
    #token = json.loads(response.content)['token']['id']
    '''
    # subsetting and reformatting services    
    # Query service capability URL 
    #from xml.etree import ElementTree as ET
    capability_url = f'https://n5eil02u.ecs.nsidc.org/egi/capabilities/{short_name}.{latest_version}.xml'
    #print(capability_url)
    
    #Create session to store cookie and pass credentials to capabilities url
    uid, pswd, token = get_api_key()
    # If above line does not work, enter your credentials manually as below
    #uid = 'enter '
    #email = ''
    #pswd = ''
    #pswd = getpass.getpass('Earthdata Login password: ') #B2f
    #token = os.environ['EARTHDATA_TOKEN'] #'WSADG-GEGC-sE887-E43F-2DS9SD-49KDF'   
    session = requests.session()
    s = session.get(capability_url)
    response = session.get(s.url,auth=(uid,pswd))
    logging.info(f'Server Response = {response}')
    # New check (Dec 19, 2020) but not checked yet in script; to exit here if server response is not 200
    if response.status_code != 200:
        logging.error(f'Exiting Script because server response of {response} is not 200')
        sys.exit(0)

    #logging.info('Before Assert')
    #assert(response.status_code == 200) #This could trigeer unknown error
    #logging.info('After Assert')
    root = ET.fromstring(response.content)    
    # collect lists with each service option
    subagent = [subset_agent.attrib for subset_agent in root.iter('SubsetAgent')]
    # variable subsetting
    variables = [SubsetVariable.attrib for SubsetVariable in root.iter('SubsetVariable')]
    variables_raw = [variables[i]['value'] for i in range(len(variables))]
    variables_join = [''.join(('/',v)) if v.startswith('/') == False else v for v in variables_raw] 
    variable_vals = [v.replace(':', '/') for v in variables_join]
    
    # reformatting
    formats = [Format.attrib for Format in root.iter('Format')]
    format_vals = [formats[i]['value'] for i in range(len(formats))]
    format_vals.remove('')
    
    #print(subagent)
    if len(subagent) < 1 :
        agent = 'NO'
        print(agent)
    #print(format_vals)
    #print(len(variable_vals))
    
    # Coverage is a required field for api: also check if possible without coverage field; so it downloads everything for that Ground reference track
    coverage = ''
#    if short_name == 'ATL06':
#        coverage = '/ancillary_data/atlas_sdp_gps_epoch,/gt1l/land_ice_segments/atl06_quality_summary,/gt1l/land_ice_segments/delta_time,/gt1l/land_ice_segments/h_li,/gt1l/land_ice_segments/h_li_sigma,/gt1l/land_ice_segments/latitude,/gt1l/land_ice_segments/longitude,/gt1l/land_ice_segments/segment_id,/gt1l/land_ice_segments/sigma_geo_h,/gt1r/land_ice_segments/atl06_quality_summary,/gt1r/land_ice_segments/delta_time,/gt1r/land_ice_segments/h_li,/gt1r/land_ice_segments/h_li_sigma,/gt1r/land_ice_segments/latitude,/gt1r/land_ice_segments/longitude,/gt1r/land_ice_segments/segment_id,/gt1r/land_ice_segments/sigma_geo_h,/gt2l/land_ice_segments/atl06_quality_summary,/gt2l/land_ice_segments/delta_time,/gt2l/land_ice_segments/h_li,/gt2l/land_ice_segments/h_li_sigma,/gt2l/land_ice_segments/latitude,/gt2l/land_ice_segments/longitude,/gt2l/land_ice_segments/segment_id,/gt2l/land_ice_segments/sigma_geo_h,/gt2r/land_ice_segments/atl06_quality_summary,/gt2r/land_ice_segments/delta_time,/gt2r/land_ice_segments/h_li,/gt2r/land_ice_segments/h_li_sigma,/gt2r/land_ice_segments/latitude,/gt2r/land_ice_segments/longitude,/gt2r/land_ice_segments/segment_id,/gt2r/land_ice_segments/sigma_geo_h,/gt3l/land_ice_segments/atl06_quality_summary,/gt3l/land_ice_segments/delta_time,/gt3l/land_ice_segments/h_li,/gt3l/land_ice_segments/h_li_sigma,/gt3l/land_ice_segments/latitude,/gt3l/land_ice_segments/longitude,/gt3l/land_ice_segments/segment_id,/gt3l/land_ice_segments/sigma_geo_h,/gt3r/land_ice_segments/atl06_quality_summary,/gt3r/land_ice_segments/delta_time,/gt3r/land_ice_segments/h_li,/gt3r/land_ice_segments/h_li_sigma,/gt3r/land_ice_segments/latitude,/gt3r/land_ice_segments/longitude,/gt3r/land_ice_segments/segment_id,/gt3r/land_ice_segments/sigma_geo_h,/orbit_info/cycle_number,/orbit_info/rgt,/orbit_info/orbit_number'
    
    # Temporal subsetting KVP
    timevar = start_date + 'T' + start_time + ',' + end_date + 'T' + end_time
    #timevar = temporal #see if this works [not yet]
    # Request data from the NSIDC data access API.
    ''' The API is structured as a URL with a base plus individual key-value-pairs (KVPs) separated by ‘&’. 
        The base URL of the NSIDC API is: '''
    base_url = 'https://n5eil02u.ecs.nsidc.org/egi/request' ##Set NSIDC data access base URL
    
    # Set number of granules requested per order, which we will initially set to 10.
    page_size = 10
    #Determine number of pages basd on page_size and total granules. Loop requests by this value
    page_num = math.ceil(len(granules)/page_size)
    #Set request mode. 
    request_mode = 'async' #with synchronous, there is possibility of timeout errors
    
    #token = os.environ['EARTHDATA_TOKEN']
    email = ''
    
    #Create config dictionary
    config_params = {
        'request_mode': request_mode, 
        'page_size': page_size,  
        'token': token, 
        'email': email,   
    }
    
    # Determine how many individual orders we will request based on the number of granules requested
    #print(f'Page numbers = {page_num}')
    ##Print API base URL + request parameters : But this is not used anywhere
    #API_request = f'{base_url}?short_name={short_name}&version={latest_version}&temporal={temporal}&time={timevar}&polygon={polygon}&bbox={polygon}&Coverage={coverage}&request_mode={request_mode}&page_size={page_size}&page_num={page_num}&token={token}&email={email}'
    #print(API_request)
    # create a new dictionary of NSIDC API KVPs to be used in our subset request
    # Adding customization parameter dictionary 
    custom_params = {
        'time': timevar,
        'Coverage': coverage, 
        }
    # Creating final request parameter dictionary with search, config, and customization parameters.
    subset_request_params = {**search_params, **config_params, **custom_params}
    #print(subset_request_params)
    #session = requests.session() #New here, may not work
    #s = session.get(capability_url)
    #response = session.get(s.url,auth=(uid,pswd))
    #root = ET.fromstring(response.content)
    # Request data service for each page number, and unzip outputs
    no_processing_params = {
        'time': timevar,
        'Coverage': coverage,
        'agent' : 'NO',
        'include_meta' : 'Y',
    }

    
    # Creating additional final request parameter dictionary with search, config, and no processing parameters.
    native_request_params = {**search_params, **config_params, **no_processing_params}

    # Choose method to make a request
    method = 'post' #'get' 'post'
    #print(native_request_params)
    #native_request_params = subset_request_params #only temporary
    # Order data
    orderIDs = [] #To save order numbers
    for i in range(page_num):
        page_val = i + 1
        #print('Order: ', page_val)
        if method=='get':
            # shapefile is not sent to the server
            native_request_params.update({'page_num': page_val}) 
            #native_request_params.update({'bbox': polygon}) #BY: new
            # For all requests other than spatial file upload, use get function    
            request = session.get(base_url, params=native_request_params)
        #for Post
        if method == 'post':
            # We are sending shapefile (in KML format) to the server
            subset_request_params.update( {'page_num': page_val})
            # Post polygon to API endpoint for polygon subsetting to subset based on original, non-simplified KML file
            #shape_post = {'shapefile': open(kml_filepath, 'rb')}
            shape_post = {'shapefile': shp_json}
            request = session.post(base_url, params=subset_request_params, files=shape_post) 
        logging.info(f'HTTP response code: {request.status_code}')
    
        # Raise bad request: Loop will stop for bad response code.
        request.raise_for_status()
        #print('Order request URL: ', request.url)
        esir_root = ET.fromstring(request.content)
        #print('Order request response XML content: ', request.content)
    
        #Look up order ID
        orderlist = []
        for order in esir_root.findall("./order/"):
            orderlist.append(order.text)
        orderID = orderlist[0]
        logging.info(f'order ID: {orderID}')
        # This should be the end of ordering        
        orderIDs.append(orderID) # seems we are missing some order
        
        #'''
        #Create status URL
        statusURL = base_url + '/' + orderID
        #print('status URL: ', statusURL)
    
        #Find order status
        request_response = session.get(statusURL)
        logging.info(f'Order Status HTTP response: {request_response.status_code}')
        
        # Raise bad request: Loop will stop for bad response code.
        request_response.raise_for_status()
        request_root = ET.fromstring(request_response.content)
        statuslist = []
        for status in request_root.findall("./requestStatus/"):
            statuslist.append(status.text)
        status = statuslist[0]
        #print('Data request ', page_val, ' is submitting...')
        #print('Initial request status is ', status)
    
        sleep_sec = max(5, len(granules))
        time.sleep(sleep_sec)
        #Continue loop while request is still processing
        while status == 'pending' or status == 'processing': 
            #print('Status is not complete. Trying again.')
            # sleep_sec = max(5, len(granules))
            time.sleep(10)
            loop_response = session.get(statusURL)
            # Raise bad request: Loop will stop for bad response code.
            loop_response.raise_for_status()
            loop_root = ET.fromstring(loop_response.content)
            #find status
            statuslist = []
            for status in loop_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            #print('Retry request status is: ', status)
            if status == 'pending' or status == 'processing':
                continue
        #Order can either complete, complete_with_errors, or fail:
        # Provide complete_with_errors error message:
        if status == 'complete_with_errors' or status == 'failed':
            messagelist = []
            for message in loop_root.findall("./processInfo/"):
                messagelist.append(message.text)
            logging.error(f'error messages: {messagelist}')
            #pprint.pprint(messagelist)
    
        # Download zipped order if status is complete or complete_with_errors
        #'https://n5eil02u.ecs.nsidc.org/esir/5000000402535/166238470/
        #processed_ATL06_20181102070512_05290110_002_01.h5
        if status == 'complete' or status == 'complete_with_errors':
            downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
            #print('Beginning download of zipped output...')
            zip_response = session.get(downloadURL)
            # Raise bad request: Loop will stop for bad response code.
            zip_response.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                z.extractall(f'{icesat2_path}/downloads')
            #print('Data request', page_val, 'is complete.')
        else: logging.info('Request failed.')
    logging.info('Download of icesat2 complete')
    # Cleanup : move the hdf files to base folder and deleted the "downloads" folder
    move_files_from_order(icesat2_path)
    #read_atl06(icesat2_path) #moved outside to work independent of download
    #'''
    return orderIDs
    #time.sleep(1) #To delay another request to server

#%%
def get_api_key():
#    import os 
#    dirpath = os.getcwd()
#    print("current directory is : " + dirpath)
#    foldername = os.path.basename(dirpath)
#    print("Directory name is : " + foldername)
    #return os.environ['PL_API_KEY']
    #return os.environ['EARTHDATA_TOKEN']
    #from pathlib import Path
    #current_dir = Path('__file__').parents[0]
    #os.path.dirname(os.path.abspath('__file__'))
    #print(current_dir)
    #config_file = os.path.join(current_dir, 'earthdata_config.ini')
    # C:/Github/giuh/scripts/earthdata_config.ini
##    with open('earthdata_config.ini') as f:
##        for line in f.readlines():
##            line = line.strip()
##            if 'uid' in line:
##                uid = line.split('=')[-1]
##            if 'pswd' in line:
##                pswd = line.split('=')[-1]
##            if 'token' in line:
##                token = line.split('=')[-1]
    #return uid, pswd, token
    return 'earthdata_usermane', 'earthdata_pwd', 'token_than_need_to_be_replace_every_few_months'

#%%
if __name__ == "__main__":
    ''' This generic script should work for both ATL06 and ATL08
    At a minimum:
        1. Define the DEM strips location
            strips_folder
        2. Define the output location to save downloaded icesat2 files
            base_icesat2_path
    '''
    from giuh_helpers import tic, toc
    tic()
    import argparse

    parser = argparse.ArgumentParser(description='Download Icesat2 Data for DEM strips.')
    parser.add_argument('region', help='region folder name', type=str)
    parser.add_argument('dem_type', help='REMA ArcticDEM EarthDEM', type=str)
    parser.add_argument('--short_name', help='Icesat2 product shortname', type=str, default='ATL06')

    args = parser.parse_args()
    # dem_region = args.dem_region  # rema or arcticdem #antarctica or arctic
    # start_date = args.start_date  # '2020-06-01'  # Same for all latitudes
    short_name = args.short_name  # ATL06 or ATL08 
    region = args.region  # sys.argv[1] #pass region from outside
    dem_type =args.dem_type # dem_type = ['REMA', 'ArcticDEM', 'EarthDEM']
    print(args)    
    logging.basicConfig(filename=f'{short_name}_{region}.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')
    # Location of where DEMs are staged
    # region = arcticdem_05_greenland_northeast arcticdem_06_greenland_northwest arcticdem_05_greenland_northeast arcticdem_04_greenland_central arcticdem_05_greenland_northeast  
    # region = rema_01_subantarctic_islands region_03_peninsula_south hma_2019jun26 region_20_ross_shelf
    if dem_type == 'ArcticDEM':
        dir_prefix = '/fs/byo/howat-data5'
        strip_version = "v4.1"
    elif dem_type == 'REMA':
        dir_prefix = '/fs/byo/howat-data4' # /fs/project/howat.4 : OLD
        strip_version = "v4"
    elif dem_type == "EarthDEM":
        dir_prefix = "/fs/project/howat.4"  #/EarthDEM/region_31_alaska_south/strips_v4/2m
        strip_version = "v4"
    else:
        print('Need either ArcticDEM or REMA as input for demtype, EXITING script')
        sys.exit(0)
    # strip_version = ['v4', 'unf', 'v4.1']
    # strip_version = strip_version[0] #choose based on how folder names look
    #f'{dir_prefix}/EarthDEM/{region}/strips_unf/2m' #for Alaska but not tried yet as ATL06 does not exist for this region
    # strips_folder = f'{dir_prefix}/{dem_type}/{region}/strips_v4.1/2m'  #ArcticDEM/Greenland new data 
    # strips_folder = f'{dir_prefix}/{dem_type}/{region}/strips_v4/2m'  #REMA but still old data 
    strips_folder = f'{dir_prefix}/{dem_type}/{region}/strips_{strip_version}/2m'  # Can be used for all EarthDEM, REMA, and ArcticDEM 

    strips = os.listdir(f'{strips_folder}')
    strips = [f for f in strips if not f.endswith('.json')]  # Mar 30, 2022 (for new version of DEMs)

    logging.info(f'Total Strips for {region}= {len(strips)}')

    # Set the icesat2_path : where files will be downloaded and processed (created inside code)
    base_icesat2_path = f'/fs/project/howat.4/icesat2/{dem_type}/{region}_{short_name}'
    # base_icesat2_path = f'/fs/project/howat.4/icesat2/prototyping/temp/{dem_type}/{region}_{short_name}'
    
    # Get only WV3 files for now
    cutoff_dt = pd.to_datetime('2018-10-16') #This is the date when ICESAT started collecting data
    strips = [strip for strip in strips if pd.to_datetime(strip.split('_')[1])>=cutoff_dt]
    # Aside: next two lines for testing only (March 31, 2022)
    # strips = [strip for strip in strips if strip.startswith('W3')]
    # strips = strips[20:120]

    logging.info(f'Strips after 2018-10-16 for {region}= {len(strips)}')
    # To procwss just one or multiple named strips uncomment this line below
    #strips = ['W1W1_20190410_1020010081378B00_10200100850FB100_2m_lsf']
    download_error = []
    hdf_error = []

    strip_index = 0
    strip_order_dict = {} #To save order numbers for each Tile

    # # Aside [May 12, 2022]: Only download data for DEMs that has only one segment per DEM (for Ians proposal)
    # #--------------------------------------------------------------------------
    # one_seg_strip_list = []
    # for strip in strips:
    #     strip_folder = f'{strips_folder}/{strip}'
    #     files = os.listdir(strip_folder)
    #     meta_files = [f for f in files if f.endswith("_meta.txt")]
    #     if len(meta_files) == 1:
    #         one_seg_strip_list.append(strip)
    # sensor_list = [s.split("_")[0] for s in one_seg_strip_list]
    # sensor_list = list(set(sensor_list))
    # master_list = []
    # for sensor in sensor_list:
    #     subset = [s for s in one_seg_strip_list if sensor in s]
    #     # Adaptive subset based on how much data we have per each sensor type
    #     if len(subset)<10:
    #         master_list = master_list + subset
    #     # elif len(subset)<30:
    #     #     master_list = master_list + subset[::3]
    #     # elif len(subset)<50:
    #     #     master_list = master_list + subset[::3]
    #     # elif len(subset)<100:
    #     #     master_list = master_list + subset[::5]
    #     else:
    #         master_list = master_list + subset[::5]  #Get every 10th strip
    #     logging.info(f"sensor: {sensor}, count: {len(subset)}")
    # logging.info(f"master_list: {len(master_list)}, one_seg_strip_list: {len(one_seg_strip_list)}")
    # strips = master_list # one_seg_strip_list #
    # logging.info(f'Strips with one segment DEM for {region}= {len(strips)}')
    # #--------------------------------------------------------------------------
    # Aside: these two stips requested by Ian for region_34 only
    # strips.append("WV02_20190226_103001008ED5C100_103001008B691F00_2m_lsf_v030403")
    strips = ["WV03_20190204_1040010045934B00_10400100484A3900_2m_lsf_v030403"]  #region_31

    for strip in strips:
        logging.info('==============================================================================================================')
        strip_index += 1 # just to keep count of strips
        logging.info(f'{strip_index}/{len(strips)}   {region}:: {strip}')
        strip_folder = f'{strips_folder}/{strip}'
        # # Aside [May 12, 2022]: Only download data for DEMs that has only one segment per DEM (for Ians proposal)
        # #--------------------------------------------------------------------------
        # files = os.listdir(strip_folder)
        # meta_files = [f for f in files if f.endswith("_meta.txt")]
        # if len(meta_files) > 1:
        #     # Don't process if more than one segment for DEM 
        #     continue
        # #--------------------------------------------------------------------------
        
        icesat2_path = f'{base_icesat2_path}/{strip}'
        if os.path.exists(icesat2_path):
            logging.info(f'Skipping Download and Processing {region}: {strip}')
            # TODO: Also check if downloads subfolder exit: this implies 
            # read_atl03(icesat2_path)  # ASIDE: May 12, 2022
            continue
        try:
            # Download icesat2 data for a REMA strip using API
            # Remove trailing _v0xxxx for new ArcticDEM files
            #strip = strip.split('_v0')[0]
            # In the newer version of ArcticDEM for Greenland, there are some missing strip even though the folder exist
            if len(os.listdir(f'{strip_folder}'))<1:
                # Check if no DEM exist; otherwise the script will crash when downloading
                print(f'MISSING STRIP: {strip_folder}')
                logging.error(f'MISSING STRIP: {strip_folder}')
                continue
            orderdict = download_process_icesat2(strip_folder, strip, icesat2_path)
            strip_order_dict[strip] = orderdict
            # time.sleep(1) #To delay another request to server
        except:
            logging.error(f'Strip_Index: {strip_index}  Exception in data download: {region} {strip}')
            download_error.append(f'{strip}')
        #Parsing-----------------------------------------------------
        # Think if its better to parse separately
        try:
            # Parse HDF file: convert to csv and shapefile/geopackage
            if short_name == 'ATL06':
                read_atl06(icesat2_path)
            if short_name == 'ATL08':
                read_atl08(icesat2_path)
            if short_name == 'ATL03':
                pass  # Dont parse yet
                # read_atl03(icesat2_path)
        except:
            logging.error(f' Strip_Index: {strip_index}  Exception in processing hdf file: {region} {strip}')
            hdf_error.append(f'{strip}')
        logging.info(f'End for strip_index = {strip_index}\n')
        #Parsing----------------------------------------------------------------------------------------------------------------------------------        
    if len(download_error)>0:
        logging.error('The following list contains Download Errors. Diagnose them separately')
        logging.info(download_error) #save to list of error strips to diagnose later
    if len(hdf_error)>0:
        logging.error('The following list contains HDF error. Diagnose them separately')
        logging.info(hdf_error) #save to list of error strips to diagnose later
    running_time = toc()
    logging.info(f'Total Running Time: {running_time}')    
    # Save the orders returned by the nsidc server as csv file
    #df = pd.DataFrame.from_dict(strip_order_dict, orient='index')
    #df.to_csv(f'{base_icesat2_path}/orders.csv')
    print(f"Completed: Search and Download icesat for {short_name}")
    #curl -O -J --dump-header response-header.txt "https://n5eil02u.ecs.nsidc.org/egi/request?short_name=ATL06&version=001&bounding_box=-115,-75.3,-111,-74&
    #bbox=-115,-75.3,-111,-74&time=2019-01-01T00:00:00,2019-01-31T23:59:59&token=4E3042A0-E3B6-FEE0-417B-F69CA0C481C6"
