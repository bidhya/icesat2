#!/usr/bin/env python
# coding: utf-8

""" Download ATL06/ATL08 data for summer months by 100 KM PGC tiles for either antarctica or greenland
    
    The scrip will run for just one summer. So, for different year summertime 
    data, call this scipt separately for each year.

    Updates:
    May 08, 2022: Downloading latest v5 for Antarctica
    May 09, 2022: renamed "icesat2_summer_download.py"

    Known Issues
    ------------
    July 01, 2021 : Updated to solve problem for tiles around +-180 longitude (using -ve buffer)
                    Only for Antarctia

    Inputs:
        - dem region
        - start date
        - shapefile of pgc tiles (hardcoded)
        Hardcoded:
            - 100 km grid shapefile for REMA and Greenland (from Ian Howat)
        End date will be adaptively determined by script based on latitude

    Recommended to use the following time for calling the script
        Jun-01 : Summer start date for Greenland
        Dec-01  : Summer start data for Antarctica

    Known Issue: Basically, a server error. Restart the script, it will start where the script failed 
                    The folder where it failed can also be deleted to be sure that the folder is processed
                    For post processing: check if any folder still has /downloads subfolder

        Traceback (most recent call last):
        File "/home/yadav.111/Github/giuh/scripts/icesat_summer_by_tile.py", line 146, in <module>
            # bbox = bounding_box #not used in post
        File "/home/yadav.111/miniconda3/lib/python3.9/json/__init__.py", line 346, in loads
            return _default_decoder.decode(s)
        File "/home/yadav.111/miniconda3/lib/python3.9/json/decoder.py", line 337, in decode
            obj, end = self.raw_decode(s, idx=_w(s, 0).end())
        File "/home/yadav.111/miniconda3/lib/python3.9/json/decoder.py", line 355, in raw_decode
            raise JSONDecodeError("Expecting value", s, err.value) from None
        json.decoder.JSONDecodeError: Expecting value: line 1 column 1 (char 0)
"""

import math
import os
import shutil
import sys
import time
from datetime import timedelta
# import numpy as np
import pandas as pd
import geopandas as gpd
# from pandas.core.accessor import delegate_names
import json
import zipfile
import io
import requests
import logging
# Other custom modules specific to icesat-2
from icesat2_download import move_files_from_order, read_atl06, read_atl08, get_api_key
import argparse
import warnings
warnings.filterwarnings("ignore", category=UserWarning)  # geopandas creating centroid based on lat/lon rather than projected gdf

parser = argparse.ArgumentParser(description='Download Icesat2 Data from Summer by fixed size tiles.')
parser.add_argument('dem_region', help='antarctica or greenland', type=str)
parser.add_argument('start_date', help='Starting Date', type=str)
parser.add_argument('--short_name', help='Icesat2 product shortname', type=str, default='ATL06')

args = parser.parse_args()
dem_region = args.dem_region  # antarctica or arctic
start_date = args.start_date  # '2020-06-01'  # Same for all latitudes
short_name = args.short_name  # ATL06 or ATL08

# start_date = '2020-06-01'  # This will remain same for all latitudes
logging.basicConfig(filename=f'{dem_region}_summer_{start_date[:4]}.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')
logging.info(f'{dem_region} {start_date} {short_name} new UPDATE')

# 1. Choose output Directory
if dem_region == 'antarctica':
    tile_gdf = gpd.read_file('/fs/project/howat.4/icesat2/ancillary/REMA_Tile_Index/REMA_Tile_Index_Rel1_1.shp')
    tile_gdf.name = tile_gdf.tile  # trick to get rid of _8m suffix so that data is reconsided with greenland tiles
elif dem_region == 'greenland':
    tile_gdf = gpd.read_file('/fs/project/howat.4/icesat2/ancillary/ArcticDEMTiles_Greenland/ArcticDEMTiles_Greenland.shp')
else:
    print('Need either antarctica or greenland as input for demtype, exiting script')
    sys.exit(0)
# tile_gdf.crs = 'EPSG:3413'  # if prj file was not supplied
# tile_gdf = tile_gdf.to_crs('EPSG:4326')

# tile_gdf['lat'] = tile_gdf.geometry.centroid.y  # latitude of centroid, but not useful anymore as data on unprojected to lat/lon
# tile_gdf['row'] = tile_gdf.name.apply(lambda x: x.split('_')[0])
# tile_gdf['col'] = tile_gdf.name.apply(lambda x: x.split('_')[1]) #we will use col to subset problem tiles around -180, +180 longitude

# Select one region from Antarctica
regions = list(tile_gdf.name)  # slightly different naming of tiles for Greenland
# regions.remove("30_14")  # due to unkown error May 10, 2023  
# region = [reg for reg in regions if reg.startswith(f'rema_{reg_number}_')][0]
# July 1, 2021: Some tiles around -180, +180 degree is confused by geopandas, so process these separately by using -ve buffer (say 100m for pgc tiles)
# valid only for Antarctica
# problem_regions = ['17_30_8m', '18_30_8m', '19_30_8m', '20_30_8m', '21_30_8m', '22_30_8m', '23_30_8m', '24_30_8m', '25_30_8m', '26_30_8m', '27_30_8m', '28_30_8m', '29_30_8m']
problem_regions = ['17_30', '18_30', '19_30', '20_30', '21_30', '22_30', '23_30', '24_30', '25_30', '26_30', '27_30', '28_30', '29_30']
# start_date = '2020-06-01'  # This will remain same for all latitudes

total_tiles = len(regions)
tile_cnt = 0
# START LOOP TO PULL DATA FOR EACH TILE
# for region in regions[start_idx:end_idx]:
for region in regions:
    gdf = tile_gdf[tile_gdf.name == region]  # select one polygon
    if region in problem_regions and dem_region == 'antarctica':
        logging.info('Data along +-180 degrees seam, using -100m buffer')
        gdf = gdf.buffer(-100)  # WARN:  only for problematic antarctica tiles
    gdf = gdf.to_crs('EPSG:4326')
    lat = gdf.centroid.y.values[0]
    logging.info(f'Latitude = {lat}')
    # Use absolute value of latitude for consistency along North and South Poles
    # if abs(lat) > 80:
    #     # download these files later, perhaps because it takes longer
    #     logging.info("Skipping download below 80 degree lat")
    #     continue
    tile_cnt += 1
    logging.info(f'Starting : {region} : {tile_cnt}/{total_tiles}----------------------------------------------------------')
    # icesat2_path = f'/fs/project/howat.4/icesat2/{dem_region}_tiles_{short_name}/{start_date[:4]}/{region}'
    icesat2_path = f'/fs/project/howat.4/icesat2/{dem_region}_summer_{short_name}/{start_date[:4]}/{region}'
    if not os.path.exists(icesat2_path):
        logging.info(f'Create Output Directory : {icesat2_path}')
        os.makedirs(icesat2_path)  # exist_ok=True to prevent complaint if directory exist
        os.mkdir(f'{icesat2_path}/downloads')
    else:
        # Added Dec 24, 2020 to ignore data for PGC tiles are already downloaded
        if os.path.exists(f'{icesat2_path}/downloads'):
            # Updated May 07, 2022: Clean download again where previous script failed by first deleting the downloads folder
            logging.info('Aside: DOWNLOAD AGAIN WHERE IT FAILED DUE TO SERVER')
            shutil.rmtree(f'{icesat2_path}/downloads')
            logging.info('Aside: Deleting and creating downloads leaf-folder')
            os.mkdir(f'{icesat2_path}/downloads')
        else:
            logging.info('Previously downloaded')
            continue
    if abs(lat) <= 65:
        days_offset = 91  # end_date = '2020-08-30'     # 3 months, farthest from pole TODO generate end month by adding number of months/days
    elif abs(lat) > 65 and abs(lat) <= 75:
        days_offset = 60
    elif abs(lat) > 75 and abs(lat) <= 80:  # change 80 to 82
        days_offset = 30  # Closest to pole, select smaller end time
    elif abs(lat) > 80 and abs(lat) <= 85:
        # TODO: May 07, 2022: Increase offset by 2 to 4 days because one in greenland did not have any data for 2019, and 2020
        days_offset = 7  # Closest to pole, select smaller end time
    elif abs(lat) > 85:
        days_offset = 3  # Closest to pole, select smaller end time

    # Prototyping different time-offset schemes (but has to be implemented uniformly across all regions, try when v6 data is available)
    # if abs(lat) <= 65:
    #     days_offset = 91  # end_date = '2020-08-30'     # 3 months, farthest from pole TODO generate end month by adding number of months/days
    # elif abs(lat) > 65 and abs(lat) <= 70:
    #     days_offset = 60
    # elif abs(lat) > 70 and abs(lat) <= 75:
    #     days_offset = 50
    # elif abs(lat) > 75 and abs(lat) <= 82:  # change 80 to 82
    #     days_offset = 30  # Closest to pole, select smaller end time
    # elif abs(lat) > 80 and abs(lat) <= 85:
    #     # TODO: May 07, 2022: Increase offset by 2 to 4 days because one in greenland did not have any data for 2019, and 2020
    #     days_offset = 7  # Closest to pole, select smaller end time
    # elif abs(lat) > 85:
    #     days_offset = 3  # Closest to pole, select smaller end time

    end_date = pd.to_datetime(start_date) + timedelta(days=days_offset)  # '2020-02-28' directly
    end_date = end_date.strftime('%Y-%m-%d')
    temporal = f'{start_date} , {end_date}'
    logging.info(f'Date Range: {temporal}')
    shp_json = gdf.to_json()

    # To use Bounding Box for subsetting
    bounding_box = ','.join(map(str, list(gdf.total_bounds)))
    # Bounding box subsetting (bbox) in same format as bounding_box
    # bbox = bounding_box #not used in post
    # Run this occasionally to check the latest version
    search_params = {'short_name': short_name}
    cmr_collections_url = 'https://cmr.earthdata.nasa.gov/search/collections.json'
    response = requests.get(cmr_collections_url, params=search_params)
    results = json.loads(response.content)

    # Dates for REMA: Oct 2018 to Apr 2019 Austral Summer?
    # Dates for Greenland 2019-08-01 to 2019-10-01 
    # --------------------------------------------------------------------------------------------------------------------------------------
    # Define the Granules search parameters the NSIDC Server
    latest_version = '005'
    search_params = {
            'short_name': short_name,
            'version': latest_version,
            'temporal': temporal,
            'page_size': 100,
            'page_num': 1,
            'bounding_box': bounding_box
            }
    granule_search_url = 'https://cmr.earthdata.nasa.gov/search/granules'
    headers = {'Accept': 'application/json'}  # if not defined in first block to get the token
    granules = []
    count = 0
    while True:
        #logging.info(count) #seems to be used for debugging
        response = requests.get(granule_search_url, params=search_params, headers=headers)
        results = json.loads(response.content)
        if len(results['feed']['entry']) == 0:
            # Out of results, so break out of loop
            break
        # Collect results and increment page_num
        granules.extend(results['feed']['entry'])
        search_params['page_num'] += 1
        count += 1
    # Get number of granules over my area and time of interest
    logging.info(f'Number of Granules = {len(granules)}')
    granule_sizes = [float(granule['granule_size']) for granule in granules]
    logging.info(f'{len(granule_sizes)} Total Granules of Size: {sum(granule_sizes):.2f} MB')

    # ## 5. Choose the authentication mechanism
    ## 6. Request token from Common Metadata Repository using Earthdata credentials
    # Query service capability URL : seems only to see what service is availabe
    from xml.etree import ElementTree as ET
    capability_url = f'https://n5eil02u.ecs.nsidc.org/egi/capabilities/{short_name}.{latest_version}.xml'
    email = ''
    # uid = 'earthdata_uname'
    # pswd = 'earthdata_password'
    # token = 'xxx3CB87-39DA-CC37-DCE5-89A3D4A60sss'
    uid, pswd, token = get_api_key()

    # Create session to store cookie and pass credentials to capabilities url
    session = requests.session()
    s = session.get(capability_url)

    response = session.get(s.url,auth=(uid,pswd))
    logging.info(f'Server Response = {response}')
    # TODO: The script should exit here if server response is not 200
    if response.status_code != 200:
        logging.error(f'Exiting Script because server response of {response} is not 200')
        sys.exit(0)

    coverage = ''
    # Temporal subsetting KVP
    timevar = temporal  # this seemed to work as well, when used temporal below directly
    # Request data from the NSIDC data access API.
    '''the API is structured as a URL with a base plus individual key-value-pairs (KVPs) separated by ‘&’. 
    The base URL of the NSIDC API is: '''
    base_url = 'https://n5eil02u.ecs.nsidc.org/egi/request'  # Set NSIDC data access base URL

    # Set number of granules requested per order, which we will initially set to 10.
    page_size = 20  # 10
    # Determine number of pages basd on page_size and total granules. Loop requests by this value
    page_num = math.ceil(len(granules)/page_size)
    logging.info(f'page_num = {page_num}')
    #Set request mode. 
    request_mode = 'async' #with synchronous, there is possibility of timeout errors

    #Create config dictionary
    config_params = {
        'request_mode': request_mode, 
        'page_size': page_size,  
        'token': token, 
        'email': email,   
    }

    #timevar replaced with temporal: and it seems to work fine
    custom_params = {
        'time': temporal,
        'Coverage': coverage, 
        }
    # Creating final request parameter dictionary with search, config, and customization parameters.
    subset_request_params = {**search_params, **config_params, **custom_params}

    # Choose method to make a request
    method = 'post' #'get' 'post'
    for i in range(page_num):
        page_val = i + 1
        logging.info(f'Order: {page_val}')
        if method == 'post':
            # We are sending shapefile (in KML format) to the server
            subset_request_params.update( {'page_num': page_val})
            # Post polygon to API endpoint for polygon subsetting to subset based on original, non-simplified KML file
            #shape_post = {'shapefile': open(kml_filepath, 'rb')}
            shape_post = {'shapefile': shp_json}        
            request = session.post(base_url, params=subset_request_params, files=shape_post, auth=(uid,pswd))
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

        #Create status URL
        statusURL = base_url + '/' + orderID
        logging.info(f'status URL: {statusURL}')

        #Find order status
        request_response = session.get(statusURL)
        #logging.info(f'Order Status HTTP response: {request_response.status_code}')

        # Raise bad request: Loop will stop for bad response code.
        request_response.raise_for_status()
        request_root = ET.fromstring(request_response.content)
        statuslist = []
        for status in request_root.findall("./requestStatus/"):
            statuslist.append(status.text)
        status = statuslist[0]
        logging.info(f'Data request {page_val} is submitting...')
        logging.info(f'Initial request status is {status}')
        # Continue loop while request is still processing
        while status == 'pending' or status == 'processing': 
            logging.info('Status is not complete. Trying again.')
            sleep_sec = 10
            time.sleep(sleep_sec)
            loop_response = session.get(statusURL)
            # Raise bad request: Loop will stop for bad response code.
            loop_response.raise_for_status()
            loop_root = ET.fromstring(loop_response.content)
            #find status
            statuslist = []
            for status in loop_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            logging.info(f'Retry request status is: {status}')
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
            logging.info('Beginning download of zipped output...')
            zip_response = session.get(downloadURL)
            # Raise bad request: Loop will stop for bad response code.
            try:
                # Added Dec 24, 2020 because of error due to completing the order but having no zip file to download (perhaps no data due to spatial/temporal subsetting)
                zip_response.raise_for_status()
                with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                    z.extractall(f'{icesat2_path}/downloads')
                logging.info(f'Data request {page_val} is complete.')
            except:
                logging.error('Zip file download error for orderID: {downloadURL}')
        else:
            logging.info('Request failed.')

    # Further processing
    # Cleanup : move the hdf files to base folder and deleted the "downloads" folder
    move_files_from_order(icesat2_path)
    logging.info('End Download and Extraction')
    time.sleep(3) #Before downloading next series of pbs command for download

    # # Convert/Parse the HDF files to csv file/shapefile
    if short_name == 'ATL08':
        logging.info(f'Parsing {short_name} hdf file--------------------------------------------')
        atl_gdf = read_atl08(icesat2_path)
    if short_name == 'ATL06':
        logging.info(f'Parsing {short_name} hdf file---------------------------------------------')
        atl_gdf = read_atl06(icesat2_path)
    logging.info('---------------------------------------------------------------------------------------\n')
logging.info('=======================================================================================\n')
