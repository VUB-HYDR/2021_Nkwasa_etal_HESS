#!/usr/bin/env python
'''this script generates the plant/harvest dates + irrigation + fertilization 
        for rainfed & irrigated hrus for any specific period of time from global datasets

Author  : albert nkwasa
Contact : albert.nkwasa@vub.be
Date    : 2021.04.13

'''

# Import modules

import pandas as pd
import xarray as xr
import numpy as np
import os
import shutil
import itertools
import glob
import datetime
import sys
import numpy as np

print('\t >Rememebr to have the correct crop maps in the correct folder')
now = datetime.datetime.now()

# model file name
datazone = 'africa'  # subject to change

# define the hru land use codes for rainfed and irrigated
non_irrigated_lum = 'agrl_lum'  # subject to change
irrigated_lum = 'agrr_lum'  # subject to change

# setting the working directory
wkdir = '/{}/Scenarios/Default'.format(datazone)  # subject to change
try:
    os.chdir(wkdir)
except:
    print("\t> the specified model at {0} is invalid".format(wkdir))


# check if the management folder exists and delete (folder for the new management files)
try:
    shutil.rmtree('management_tables')
except:
    pass

# make a new folder for the new mangement files
try:
    os.makedirs("management_tables")
except:
    pass

os.chdir("management_tables")

# setting the hru paths
# hru landuses specifically agrl
# subject to change
file_1 = '/{}/Scenarios/Default/TxtInOut/hru-data.hru'.format(datazone)
df_1 = pd.read_csv(file_1, delimiter='\s+', skiprows=1)
df_1 = df_1.drop(['topo', 'hydro', 'soil', 'soil_plant_init',
                  'surf_stor', 'name', 'id', 'snow', 'field'], axis=1)
# hru coordinates
# subject to change
file_2 = '/{}/Scenarios/Default/TxtInOut/hru.con'.format(datazone)
df_2 = pd.read_csv(file_2, delimiter='\s+', skiprows=1)
df_2 = df_2.drop(['id', 'name', 'gis_id', 'area', 'elev',
                  'wst', 'cst', 'ovfl', 'rule', 'out_tot'], axis=1)

df = pd.concat([df_1, df_2], axis=1)

# extract the hru coordinates


def get_coordinates(lum_use):
    '''
    this function returns hrus and the corresponding coordinates
    '''
    df_lum = df[df['lu_mgt'] == lum_use]
    df_lum = df_lum.reset_index(drop=True)
    lats = df_lum['lat'].tolist()
    lons = df_lum['lon'].tolist()
    col_hrus = df_lum['hru'].tolist()
    lats_lons = list(zip(lats, lons))
    coordinates = dict(zip(col_hrus, lats_lons))
    return coordinates, col_hrus


hru_coords, col_hru = get_coordinates(non_irrigated_lum)
hru_coords_irr, col_hru_irr = get_coordinates(irrigated_lum)

# extract the hru coordinates
# subject to change; use the landuse maps generated
crop_path = '/landuse_map_folder'


def crop_extraction(crop_coord, path_agrl_lu):
    '''
    this extracts the specific crops from the defined landuse_year map
    '''
    print('\t >extracting crops from landuse map')
    crop_dic = {}
    os.chdir(path_agrl_lu)
    for i in crop_coord:
        crop_list = []
        for k in os.listdir():
            crop_array = xr.open_rasterio(k)
            crop_values = crop_array.sel(
                x=crop_coord[i][1], y=crop_coord[i][0], method="nearest")
            crop_list.append(crop_values.values)
        merged_crops = list(itertools.chain(*crop_list))
        # caters for no data values (e.g -9999 or 3.48e+38 or nan values)  # no data values are assigned bare cropland
        merged_crops_clean = list(map(lambda x: 6.0 if str(
            x) == 'nan' or x < 0 else x, merged_crops))
        crop_dic[i] = merged_crops_clean
    return crop_dic


crops = crop_extraction(hru_coords, crop_path)
crops_irr = crop_extraction(hru_coords_irr, crop_path)

# Paths to the plant and harvest data maps folder for the selected crops (maize, wheat, soy, banana, sugarcane)
plant_date_dir = '/plant_rainfed'  # subject to change
plant_date_dir_irr = '/plant_irrigated'  # subject to change
harvest_date_dir = '/harvest_rainfed'  # subject to change
harvest_date_dir_irr = '/harvest_irrigated'  # subject to change

# check to see if 5 maps exist in each folder for the 5 crops


def check_tiffs(path_tiffs, number_tiffs):
    if len(glob.glob1(path_tiffs, "*.tif")) != number_tiffs:
        print("\t> the number of tiffs is invalid, check plant & harvest folders")
        sys.exit()
    else:
        return path_tiffs


plant_rainfed_tiffs = check_tiffs(plant_date_dir, 5)
plant_irrigated_tiffs = check_tiffs(plant_date_dir_irr, 5)
harvest_rainfed_tiffs = check_tiffs(harvest_date_dir, 5)
harvest_irrigated_tiffs = check_tiffs(harvest_date_dir_irr, 5)

# extract the plant and harvest days for each hru from the plant nad harvest maps


def extract_plant_harv(crop_coord, date_dir, operation, hru_type):
    os.chdir(date_dir)
    date_files = {}
    print('\t >extracting {0} dates for {1} HRUs from rasters'.format(
        operation, hru_type))
    for i in crop_coord:
        date_list = []
        for k in os.listdir(date_dir):
            if k.endswith('.tif'):
                date_array = xr.open_rasterio(k)
                date_values = date_array.sel(
                    x=crop_coord[i][1], y=crop_coord[i][0], method="nearest")
                date_list.append(date_values.values)
        merged_dates = list(itertools.chain(*date_list))
        merged_dates_clean = list(
            map(lambda x: 0.0 if str(x) == 'nan' or x < 0 else x, merged_dates))
        date_files[i] = merged_dates_clean
    return date_files


plant_date_files = extract_plant_harv(
    hru_coords, plant_date_dir, 'plant', 'rainfed')
plant_date_files_irr = extract_plant_harv(
    hru_coords_irr, plant_date_dir_irr, 'plant', 'irrigated')
harvest_date_files = extract_plant_harv(
    hru_coords, harvest_date_dir, 'harvest', 'rainfed')
harvest_date_files_irr = extract_plant_harv(
    hru_coords_irr, harvest_date_dir_irr, 'harvest', 'irrigated')

# path to the fertilizer files for N and P
n_file_dir = '/nitrogen_fertilizer'
p_file_dir = '/phosphorus_fertilizer'

# extract the fertilizer for N and P for each hru


def check_fert(crop_directory, fertilizer_directory):
    crop_files = []
    fertilizer_files = []
    for k in os.listdir(crop_directory):
        crop_files.append(k)
    for h in os.listdir(fertilizer_directory):
        fertilizer_files.append(h)

    if len(crop_files) != len(fertilizer_files):
        print("\t> The number of crop maps is different from the number of fertilizer maps. CHECK!!")
        sys.exit()
    else:
        pass


check_fert(crop_path, n_file_dir)
check_fert(crop_path, p_file_dir)


def extract_fertilizer(crop_coord, fert_dir, file_name, type_fert):
    os.chdir(fert_dir)
    fery_files = {}
    print('\t >extracting {} fertilizer from rasters'.format(type_fert))
    for i in crop_coord:
        n_list = []
        for k in os.listdir(fert_dir):
            if k.startswith(file_name):
                n_fert_array = xr.open_rasterio(k)
                n_values = n_fert_array.sel(
                    x=crop_coord[i][1], y=crop_coord[i][0], method="nearest")
                n_list.append(n_values.values)
        merged = list(itertools.chain(*n_list))
        merged_clean = list(map(lambda x: 0.0 if str(x) ==
                                'nan' or x < 0 else x, merged))
        fery_files[i] = merged_clean
    return fery_files


nfery_files = extract_fertilizer(hru_coords, n_file_dir, 'nfery', 'N')
nfery_files_irr = extract_fertilizer(hru_coords_irr, n_file_dir, 'nfery', 'N')
pfery_files = extract_fertilizer(hru_coords, p_file_dir, 'pfery', 'P')
pfery_files_irr = extract_fertilizer(hru_coords_irr, p_file_dir, 'pfery', 'P')

# combining the dictionaries-rainfed
plant_harv_dict = {key: (crops[key], plant_date_files[key],
                         harvest_date_files[key]) for key in plant_date_files}
nfery_pfery_dict = {
    key: (nfery_files[key], pfery_files[key]) for key in nfery_files}
# combining the dictionaries-irrigated
plant_harv_dict_irr = {key: (crops_irr[key], plant_date_files_irr[key],
                             harvest_date_files_irr[key]) for key in plant_date_files_irr}
nfery_pfery_dict_irr = {
    key: (nfery_files_irr[key], pfery_files_irr[key]) for key in nfery_files_irr}

# functions for indexing
# subject to change
os.chdir('/{}/Scenarios/Default/management_tables'.format(datazone))


def l_padding(var, padding):
    l_padded = var.ljust(padding)
    return l_padded


def r_padding(var, padding):
    r_padded = var.rjust(padding)
    return r_padded


print('\t >writing to tables')

# start with the management file

mgt_file = 'management.sch'
mf = open(mgt_file, 'w+')
mf1 = 'management.sch: -{0}- version nkwasa 2021 {1}'.format(
    datazone,  now.strftime("%Y-%m-%d %H:%M:%S")) + '\n'
mf.write(mf1)
mf2 = l_padding('name', 27) + l_padding('numb_ops', 12) + l_padding('numb_auto', 21) + l_padding('op_typ', 12) + l_padding('mon', 12) + \
    l_padding('day', 12) + l_padding('hu_sch', 12) + l_padding('op_data1',
                                                               12) + l_padding('op_data2', 12) + l_padding('op_data3', 10) + '\n'
mf.write(mf2)


# with fertilizer and planting in the same dtl,
number_ops = 1

'''landuse.lum file'''
# landuse.lum file - for adding the new landuse mgt operation
# get a copy of the file in the wrkdir
# subject to change
original_lum_file = '/{}/Scenarios/Default/TxtInOut/landuse.lum'.format(
    datazone)
# subject to change
copy_lum_file = '/{}/Scenarios/Default/management_tables/landuse_copy.lum'.format(
    datazone)
shutil.copy2(original_lum_file, copy_lum_file)

# subject to change
lum_file = '/{}/Scenarios/Default/management_tables/landuse_copy.lum'.format(
    datazone)
with open(lum_file, 'r') as f:
    for i in range(2):
        f.readline()
    agrl_line_1 = f.readline().split()
    agrl_line_2 = f.readline().split()
    agrl_line_3 = f.readline().split()
    agrl_line_4 = f.readline().split()
    agrl_line_5 = f.readline().split()
    agrl_line_6 = f.readline().split()

# decision table for planting and harvesting
no_years = os.listdir(crop_path)
conditions = len(no_years)*3 + 1
conditions_irr = len(no_years)*3 + 2
alternatives = len(no_years)*2 + 1
alternatives_irr = len(no_years)*4 + 1
actions = len(no_years)*2 + 1 + len(no_years)*2
actions_irr = len(no_years)*2 + 2 + len(no_years)*2

# start the decision table writing
lum_dtl = 'lum.dtl'
g = open(lum_dtl, 'w+')
g.write('lum.dtl: -{0}- version nkwasa 2021 {1}'.format(datazone,
                                                        now.strftime("%Y-%m-%d %H:%M:%S")) + '\n')
# number of tables
numb_tables = str(len(col_hru) + len(col_hru_irr)) + '\n' + '\n'
g.write(numb_tables)

journal_crops = {}
community_plant = {}
list_mgt = []
list_mgt_irr = []

# reading the dictionary data and wrriting tables
# rainfed hrus
for h, n in zip(nfery_pfery_dict, plant_harv_dict):
    line1 = l_padding('name', 25) + l_padding('conds', 11) + \
        l_padding('alts', 10) + l_padding('acts', 10) + '\n'
    line2 = l_padding('fh{0}'.format(h), 29) + l_padding(str(conditions), 10) + \
        l_padding(str(alternatives), 10) + l_padding(str(actions), 10) + '\n'
    line3 = l_padding('var', 27) + l_padding('obj', 6) + l_padding('obj_num', 18) + \
        l_padding('lim_var', 19) + l_padding('lim_op', 11) + \
        l_padding('lim_const', 15)
    for j in range(1, alternatives+1):
        line3 = line3 + l_padding('alt' + str(j), 10)
    line3 = line3 + '\n'
    g.write(line1),  g.write(line2),  g.write(line3)
    for count, item in enumerate(plant_harv_dict[h][0], start=0):
        if item == 1.0:  # wheat
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(
                '{}'.format(plant_harv_dict[h][1][0])+'0000', 16) + ('   -   '*count + '   =   ') + ('   -   '*(alternatives - count - 1)) + '\n'
            community_plant['agr_{}'.format(h)] = 'swht_comm'
            journal_crops[h] = 'swht'
        elif item == 2.0:  # maize
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(
                '{}'.format(plant_harv_dict[h][1][1])+'0000', 16) + ('   -   '*count + '   =   ') + ('   -   '*(alternatives - count - 1)) + '\n'
            community_plant['agr_{}'.format(h)] = 'corn_comm'
            journal_crops[h] = 'corn'
        elif item == 3.0:  # soyabean
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(
                '{}'.format(plant_harv_dict[h][1][2])+'0000', 16) + ('   -   '*count + '   =   ') + ('   -   '*(alternatives - count - 1)) + '\n'
            community_plant['agr_{}'.format(h)] = 'soyb_comm'
            journal_crops[h] = 'soyb'
        elif item == 4.0:
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(
                '{}'.format(plant_harv_dict[h][1][4])+'0000', 16) + ('   -   '*count + '   =   ') + ('   -   '*(alternatives - count - 1)) + '\n'
            community_plant['agr_{}'.format(h)] = 'bana_comm'
            journal_crops[h] = 'bana'
        elif item == 5.0:  # sugarcane
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(
                '{}'.format(plant_harv_dict[h][1][4])+'0000', 16) + ('   -   '*count + '   =   ') + ('   -   '*(alternatives - count - 1)) + '\n'
            community_plant['agr_{}'.format(h)] = 'sugc_comm'
            journal_crops[h] = 'sugc'
        # incase there is a pixel mismatch-and bare soil is gotten, use agrl but put the plant/harvest date as (plant 100, harvest 349 - common days for unknown)
        elif item == 6.0:
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding(
                '-', 8) + l_padding('100.0000', 16) + ('   -   '*count + '   =   ') + ('   -   '*(alternatives - count - 1)) + '\n'
            community_plant['agr_{}'.format(h)] = 'agrl_comm'
            journal_crops[h] = 'agrl'
        g.write(line4)

    for count_h, item_h in enumerate(plant_harv_dict[h][0], start=0):
        if item == 1.0:  # wheat
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(
                plant_harv_dict[h][2][0])+'0000', 16) + ('   -   '*(len(plant_harv_dict[h][0])+count_h) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count_h)) + '\n'
        elif item == 2.0:  # maize
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(
                plant_harv_dict[h][2][1])+'0000', 16) + ('   -   '*(len(plant_harv_dict[h][0])+count_h) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count_h)) + '\n'
        elif item == 3.0:  # soyabean
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(
                plant_harv_dict[h][2][2])+'0000', 16) + ('   -   '*(len(plant_harv_dict[h][0])+count_h) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count_h)) + '\n'
        elif item == 4.0:
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(
                plant_harv_dict[h][2][4])+'0000', 16) + ('   -   '*(len(plant_harv_dict[h][0])+count_h) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count_h)) + '\n'
        elif item == 5.0:  # sugarcane
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(
                plant_harv_dict[h][2][4])+'0000', 16) + ('   -   '*(len(plant_harv_dict[h][0])+count_h) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count_h)) + '\n'
        # incase there is a pixel mismatch-and bare soil is gotten, use maize but put the plant/harvest date as (plant 100, harvest 349 - common days of corn)
        elif item == 6.0:
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(
                '345.0000', 16) + ('   -   '*(len(plant_harv_dict[h][0])+count) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count)) + '\n'
        g.write(line5)

    for count_yr, item_yr in enumerate(plant_harv_dict[h][0]):
        line6 = l_padding('year_rot', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(str(1 + count_yr)+'.00000',
                                                                                                                                              16) + ('   -   '*count_yr + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0])-1) + '   =   ') + ('   -   '*(len(plant_harv_dict[h][0]) - count_yr)) + '\n'
        g.write(line6)

    line7 = l_padding('year_rot', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding(
        '-', 8) + l_padding(str(len(plant_harv_dict[h][0]))+'.00000', 16) + ('   -   '*(alternatives-1) + '   >   ') + '\n'
    line8 = l_padding('act_typ', 27) + l_padding('obj', 6) + l_padding('obj_num', 21) + l_padding('name', 16) + l_padding(
        'option', 15) + l_padding('const', 10) + l_padding('const2', 14) + l_padding('fp', 9) + l_padding('outcome', 9) + '\n'

    g.write(line7), g.write(line8)

    for count_p, item_p in enumerate(plant_harv_dict[h][0], start=0):
        if item_p == 1.0:  # wheat
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_swht', 22) + l_padding('swht', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives - count_p - 1)) + '\n'
        elif item_p == 2.0:  # maize
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_corn', 22) + l_padding('corn', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives - count_p - 1)) + '\n'
        elif item_p == 3.0:  # soyb
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_soyb', 22) + l_padding('soyb', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives - count_p - 1)) + '\n'
         # remember that since i dont have dates for banana, i am using dates for sugarcane {rice dates(3) not being used[saved for asia]}
        elif item_p == 4.0:
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_bana', 22) + l_padding('bana', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives - count_p - 1)) + '\n'
        elif item_p == 5.0:  # sugarcane
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_sugc', 22) + l_padding('sugc', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives - count_p - 1)) + '\n'
        elif item_p == 6.0:  # maize-bare
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_agrl', 22) + l_padding('agrl', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives - count_p - 1)) + '\n'
        g.write(line9)

    for count_k, item_k in enumerate(plant_harv_dict[h][0], start=0):
        if item_k == 1.0:  # wheat
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('swht', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('grain', 14) + ('   n   '*(len(plant_harv_dict[h][0]) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict[h][0]) - count_k)) + '\n'
        elif item_k == 2.0:  # maize
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('corn', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('grain', 14) + ('   n   '*(len(plant_harv_dict[h][0]) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict[h][0]) - count_k)) + '\n'
        elif item_k == 3.0:  # soyb
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('soyb', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('grain', 14) + ('   n   '*(len(plant_harv_dict[h][0]) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict[h][0]) - count_k)) + '\n'
         # remember that since i dont have dates for banana, i am using dates for sugarcane {rice dates(3) not being used[saved for asia]}
        elif item_k == 4.0:
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('bana', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('grain', 14) + ('   n   '*(len(plant_harv_dict[h][0]) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict[h][0]) - count_k)) + '\n'
        elif item_k == 5.0:  # sugarcane
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('sugc', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('grain', 14) + ('   n   '*(len(plant_harv_dict[h][0]) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict[h][0]) - count_k)) + '\n'
        elif item_k == 6.0:  # maize-bare
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('agrl', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('grain', 14) + ('   n   '*(len(plant_harv_dict[h][0]) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict[h][0]) - count_k)) + '\n'
        g.write(line10)

    for count_m, item_m in enumerate(nfery_pfery_dict[h][0], start=0):
        line_m = l_padding('fertilize', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('ElementalN', 22) + l_padding('elem_n', 12) + l_padding(str(round(
            item_m*10, 0))+'0000', 12) + l_padding('1.00000', 11) + l_padding('broadcast', 14) + ('   n   '*count_m + '   y   ') + ('   n   '*(alternatives - count_m - 1)) + '\n'
        g.write(line_m)

    for count_n, item_n in enumerate(nfery_pfery_dict[h][1], start=0):
        line_n = l_padding('fertilize', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('ElementalP', 22) + l_padding('elem_p', 12) + l_padding(str(round(
            item_n*10, 0))+'0000', 12) + l_padding('1.00000', 11) + l_padding('broadcast', 14) + ('   n   '*count_n + '   y   ') + ('   n   '*(alternatives - count_n - 1)) + '\n'
        g.write(line_n)

    line11 = l_padding('rot_reset', 27) + l_padding('hru', 12) + l_padding('0', 12) + l_padding('reset_1', 21) + l_padding('null', 11) + \
        l_padding('1.00000', 11) + l_padding('0.00000', 11) + l_padding('null',
                                                                        14) + ('   n   '*(alternatives-1) + '   y   ') + '\n'+'\n'
    g.write(line11)


# managament file
    mf3 = l_padding('agrl_{}'.format(h), 34) + \
        l_padding('0', 13) + l_padding('1', 21) + '\n'
    mf4 = r_padding('fh{0}'.format(h), 66) + '\n'
    mf.write(mf3), mf.write(mf4)
# landuse file
    with open(lum_file, 'a') as l_file:
        if agrl_line_1[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_1[1], 13) + l_padding(agrl_line_1[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_1[4], 17) + l_padding(agrl_line_1[5], 25) + l_padding(
                agrl_line_1[6], 18) + l_padding(agrl_line_1[7], 10) + l_padding(agrl_line_1[8], 26) + l_padding(agrl_line_1[9], 18) + l_padding(agrl_line_1[10], 18) + l_padding(agrl_line_1[11], 18) + l_padding(agrl_line_1[12], 18) + l_padding(agrl_line_1[13], 18) + '\n'
            l_file.write(lum)
            list_mgt.append('agr_{0}'.format(h))
        elif agrl_line_2[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_2[1], 13) + l_padding(agrl_line_2[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_2[4], 17) + l_padding(agrl_line_2[5], 25) + l_padding(
                agrl_line_2[6], 18) + l_padding(agrl_line_2[7], 10) + l_padding(agrl_line_2[8], 26) + l_padding(agrl_line_2[9], 18) + l_padding(agrl_line_2[10], 18) + l_padding(agrl_line_2[11], 18) + l_padding(agrl_line_2[12], 18) + l_padding(agrl_line_2[13], 18) + '\n'
            l_file.write(lum)
            list_mgt.append('agr_{0}'.format(h))
        elif agrl_line_3[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_3[1], 13) + l_padding(agrl_line_3[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_3[4], 17) + l_padding(agrl_line_3[5], 25) + l_padding(
                agrl_line_3[6], 18) + l_padding(agrl_line_3[7], 10) + l_padding(agrl_line_3[8], 26) + l_padding(agrl_line_3[9], 18) + l_padding(agrl_line_3[10], 18) + l_padding(agrl_line_3[11], 18) + l_padding(agrl_line_3[12], 18) + l_padding(agrl_line_3[13], 18) + '\n'
            l_file.write(lum)
            list_mgt.append('agr_{0}'.format(h))
        elif agrl_line_4[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_4[1], 13) + l_padding(agrl_line_4[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_4[4], 17) + l_padding(agrl_line_4[5], 25) + l_padding(
                agrl_line_4[6], 18) + l_padding(agrl_line_4[7], 10) + l_padding(agrl_line_4[8], 26) + l_padding(agrl_line_4[9], 18) + l_padding(agrl_line_4[10], 18) + l_padding(agrl_line_4[11], 18) + l_padding(agrl_line_4[12], 18) + l_padding(agrl_line_4[13], 18) + '\n'
            l_file.write(lum)
            list_mgt.append('agr_{0}'.format(h))
        elif agrl_line_5[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_5[1], 13) + l_padding(agrl_line_5[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_5[4], 17) + l_padding(agrl_line_5[5], 25) + l_padding(
                agrl_line_5[6], 18) + l_padding(agrl_line_5[7], 10) + l_padding(agrl_line_5[8], 26) + l_padding(agrl_line_5[9], 18) + l_padding(agrl_line_5[10], 18) + l_padding(agrl_line_5[11], 18) + l_padding(agrl_line_5[12], 18) + l_padding(agrl_line_5[13], 18) + '\n'
            l_file.write(lum)
            list_mgt.append('agr_{0}'.format(h))
        elif agrl_line_6[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_6[1], 13) + l_padding(agrl_line_6[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_6[4], 17) + l_padding(agrl_line_6[5], 25) + l_padding(
                agrl_line_6[6], 18) + l_padding(agrl_line_6[7], 10) + l_padding(agrl_line_6[8], 26) + l_padding(agrl_line_6[9], 18) + l_padding(agrl_line_6[10], 18) + l_padding(agrl_line_6[11], 18) + l_padding(agrl_line_6[12], 18) + l_padding(agrl_line_6[13], 18) + '\n'
            l_file.write(lum)
            list_mgt.append('agr_{0}'.format(h))


# irrigated hrus

# irrigation water to apply to irrigated hrus; drip irrigation is used but it can be changed
irr_water = 60
for h, n in zip(nfery_pfery_dict_irr, plant_harv_dict_irr):
    line1 = l_padding('name', 25) + l_padding('conds', 11) + \
        l_padding('alts', 10) + l_padding('acts', 10) + '\n'
    line2 = l_padding('fh{0}'.format(h), 29) + l_padding(str(conditions_irr), 10) + \
        l_padding(str(alternatives_irr), 10) + \
        l_padding(str(actions_irr), 10) + '\n'
    line3 = l_padding('var', 27) + l_padding('obj', 6) + l_padding('obj_num', 18) + \
        l_padding('lim_var', 19) + l_padding('lim_op', 11) + \
        l_padding('lim_const', 15)
    for j in range(1, alternatives_irr + 1):
        line3 = line3 + l_padding('alt' + str(j), 10)
    line3 = line3 + '\n'
    g.write(line1),  g.write(line2),  g.write(line3)

    line_xx = l_padding('w_stress', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('0.90000', 16) + (
        '   -   '*(len(plant_harv_dict_irr[h][0])) + '   <   '*(len(plant_harv_dict_irr[h][0])*2)) + ('   -   '*(len(plant_harv_dict_irr[h][0])+1)) + '\n'
    g.write(line_xx)

    for count, item in enumerate(plant_harv_dict_irr[h][0], start=0):
        if item == 1.0:  # wheat
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][1][0])+'0000', 16) + (
                '   -   '*count + '   =   ') + ('   -   '*((len(plant_harv_dict_irr[h][0]) - 1) + count)) + '   >   ' + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - count)*3) + count)) + '\n'
            community_plant['agr_{}'.format(h)] = 'swht_comm'
            journal_crops[h] = 'swht_irr'
        elif item == 2.0:  # maize
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][1][1])+'0000', 16) + (
                '   -   '*count + '   =   ') + ('   -   '*((len(plant_harv_dict_irr[h][0]) - 1) + count)) + '   >   ' + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - count)*3) + count)) + '\n'
            community_plant['agr_{}'.format(h)] = 'corn_comm'
            journal_crops[h] = 'corn_irr'
        elif item == 3.0:  # soyabean
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][1][2])+'0000', 16) + (
                '   -   '*count + '   =   ') + ('   -   '*((len(plant_harv_dict_irr[h][0]) - 1) + count)) + '   >   ' + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - count)*3) + count)) + '\n'
            community_plant['agr_{}'.format(h)] = 'soyb_comm'
            journal_crops[h] = 'soyb_irr'
        elif item == 4.0:
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][1][4])+'0000', 16) + (
                '   -   '*count + '   =   ') + ('   -   '*((len(plant_harv_dict_irr[h][0]) - 1) + count)) + '   >   ' + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - count)*3) + count)) + '\n'
            community_plant['agr_{}'.format(h)] = 'bana_comm'
            journal_crops[h] = 'bana_irr'
        elif item == 5.0:  # sugarcane
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][1][4])+'0000', 16) + (
                '   -   '*count + '   =   ') + ('   -   '*((len(plant_harv_dict_irr[h][0]) - 1) + count)) + '   >   ' + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - count)*3) + count)) + '\n'
            community_plant['agr_{}'.format(h)] = 'sugc_comm'
            journal_crops[h] = 'sugc_irr'
        # incase there is a pixel mismatch-and bare soil is gotten, use agrl but put the plant/harvest date as (plant 100, harvest 349 - common days for unknown)
        elif item == 6.0:  # agrl-bare
            line4 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('100.0000', 16) + ('   -   '*count +
                                                                                                                                                                 '   =   ') + ('   -   '*((len(plant_harv_dict_irr[h][0]) - 1) + count)) + '   >   ' + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - count)*3) + count)) + '\n'
            community_plant['agr_{}'.format(h)] = 'agrl_comm'
            journal_crops[h] = 'agrl_irr'
        g.write(line4)

    for count_h, item_h in enumerate(plant_harv_dict_irr[h][0], start=0):
        if item_h == 1.0:  # wheat
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][2][0])+'0000', 16) + ('   -   '*(len(
                plant_harv_dict_irr[h][0]) + 1)) + ('   -   '*(count_h*2) + '   <   ') + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_h)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_h)) + '\n'
        elif item_h == 2.0:  # maize
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][2][1])+'0000', 16) + ('   -   '*(len(
                plant_harv_dict_irr[h][0]) + 1)) + ('   -   '*(count_h*2) + '   <   ') + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_h)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_h)) + '\n'
        elif item_h == 3.0:  # soyabean
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][2][2])+'0000', 16) + ('   -   '*(len(
                plant_harv_dict_irr[h][0]) + 1)) + ('   -   '*(count_h*2) + '   <   ') + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_h)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_h)) + '\n'
        elif item_h == 4.0:
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][2][4])+'0000', 16) + ('   -   '*(len(
                plant_harv_dict_irr[h][0]) + 1)) + ('   -   '*(count_h*2) + '   <   ') + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_h)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_h)) + '\n'
        elif item_h == 5.0:  # sugarcane
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('{}'.format(plant_harv_dict_irr[h][2][4])+'0000', 16) + ('   -   '*(len(
                plant_harv_dict_irr[h][0]) + 1)) + ('   -   '*(count_h*2) + '   <   ') + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_h)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_h)) + '\n'
        # incase there is a pixel mismatch-and bare soil is gotten, use agrl but put the plant/harvest date as (plant 100, harvest 349 - common days for unknown)
        elif item_h == 6.0:  # agrl-bare
            line5 = l_padding('jday', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding('345.0000', 16) + ('   -   '*(len(plant_harv_dict_irr[h][0]) + 1)) + (
                '   -   '*(count_h*2) + '   <   ') + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_h)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_h)) + '\n'
        g.write(line5)

    for count_yr, item_yr in enumerate(plant_harv_dict_irr[h][0]):
        line6 = l_padding('year_rot', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + l_padding(str(1 + count_yr)+'.00000', 16) + ('   -   '*count_yr + '   =   ') + (
            '   -   '*(len(plant_harv_dict_irr[h][0])-1 + count_yr)) + (('   =   ')*2) + ('   -   '*(((len(plant_harv_dict_irr[h][0]) - 1)*2)-count_yr)) + '   =   ' + ('   -   '*(len(plant_harv_dict_irr[h][0]) - count_yr)) + '\n'
        g.write(line6)

    line7 = l_padding('year_rot', 27) + l_padding('hru', 12) + l_padding('0', 15) + l_padding('null', 21) + l_padding('-', 8) + \
        l_padding(str(len(plant_harv_dict_irr[h][0]))+'.00000', 16) + (
            '   -   '*(alternatives_irr-1) + '   >   ') + '\n'
    line8 = l_padding('act_typ', 27) + l_padding('obj', 6) + l_padding('obj_num', 21) + l_padding('name', 16) + l_padding(
        'option', 15) + l_padding('const', 10) + l_padding('const2', 14) + l_padding('fp', 9) + l_padding('outcome', 9) + '\n'

    g.write(line7)
    g.write(line8)

    for count_p, item_p in enumerate(plant_harv_dict_irr[h][0], start=0):
        if item_p == 1.0:  # wheat
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_swht', 22) + l_padding('swht', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives_irr - count_p - 1)) + '\n'
        elif item_p == 2.0:  # maize
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_corn', 22) + l_padding('corn', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives_irr - count_p - 1)) + '\n'
        elif item_p == 3.0:  # soyb
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_soyb', 22) + l_padding('soyb', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives_irr - count_p - 1)) + '\n'
         # remember that since i dont have dates for banana, i am using dates for sugarcane {rice dates(3) not being used[saved for asia]}
        elif item_p == 4.0:
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_bana', 22) + l_padding('bana', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives_irr - count_p - 1)) + '\n'
        elif item_p == 5.0:  # sugarcane
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_sugc', 22) + l_padding('sugc', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives_irr - count_p - 1)) + '\n'
        elif item_p == 6.0:  # agrl-bare
            line9 = l_padding('plant', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('plant_agrl', 22) + l_padding('agrl', 12) + l_padding(
                '0.00000', 12) + l_padding('1.00000', 11) + l_padding('null', 14) + ('   n   '*count_p + '   y   ') + ('   n   '*(alternatives_irr - count_p - 1)) + '\n'
        g.write(line9)

    for count_k, item_k in enumerate(plant_harv_dict_irr[h][0], start=0):
        if item_k == 1.0:  # wheat
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('swht', 12) + l_padding('0.00000', 12) + l_padding(
                '1.00000', 11) + l_padding('grain', 14) + ('   n   '*((len(plant_harv_dict_irr[h][0])*3) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict_irr[h][0]) - count_k)) + '\n'
        elif item_k == 2.0:  # maize
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('corn', 12) + l_padding('0.00000', 12) + l_padding(
                '1.00000', 11) + l_padding('grain', 14) + ('   n   '*((len(plant_harv_dict_irr[h][0])*3) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict_irr[h][0]) - count_k)) + '\n'
        elif item_k == 3.0:  # soyb
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('soyb', 12) + l_padding('0.00000', 12) + l_padding(
                '1.00000', 11) + l_padding('grain', 14) + ('   n   '*((len(plant_harv_dict_irr[h][0])*3) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict_irr[h][0]) - count_k)) + '\n'
        elif item_k == 4.0:
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('bana', 12) + l_padding('0.00000', 12) + l_padding(
                '1.00000', 11) + l_padding('grain', 14) + ('   n   '*((len(plant_harv_dict_irr[h][0])*3) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict_irr[h][0]) - count_k)) + '\n'
        elif item_k == 5.0:  # sugarcane
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('sugc', 12) + l_padding('0.00000', 12) + l_padding(
                '1.00000', 11) + l_padding('grain', 14) + ('   n   '*((len(plant_harv_dict_irr[h][0])*3) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict_irr[h][0]) - count_k)) + '\n'
        elif item_k == 6.0:  # agrl-bare
            line10 = l_padding('harvest_kill', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('grain_harv', 22) + l_padding('agrl', 12) + l_padding('0.00000', 12) + l_padding(
                '1.00000', 11) + l_padding('grain', 14) + ('   n   '*((len(plant_harv_dict_irr[h][0])*3) + count_k) + '   y   ') + ('   n   '*(len(plant_harv_dict_irr[h][0]) - count_k)) + '\n'
        g.write(line10)

    line_yy = l_padding('irr_demand', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('drip_high', 24) + l_padding('drip', 10) + l_padding('{}'.format(irr_water)+'.00000', 12) + l_padding(
        '0.00000', 11) + l_padding('unlim', 14) + ('   n   ' * len(plant_harv_dict_irr[h][0])) + ('   y   '*(len(plant_harv_dict_irr[h][0])*2)) + ('   n   '*(len(plant_harv_dict_irr[h][0])+1)) + '\n'
    g.write(line_yy)

    for count_z, item_z in enumerate(nfery_pfery_dict_irr[h][0], start=0):
        line_z = l_padding('fertilize', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('ElementalN', 22) + l_padding('elem_n', 12) + l_padding(str(round(item_z*10, 0)) +
                                                                                                                                                           '0000', 12) + l_padding('1.00000', 11) + l_padding('broadcast', 14) + ('   n   '*count_z + '   y   ') + ('   n   '*(alternatives_irr - count_z - 1)) + '\n'
        g.write(line_z)

    for count_w, item_w in enumerate(nfery_pfery_dict_irr[h][1], start=0):
        line_w = l_padding('fertilize', 27) + l_padding('hru', 12) + l_padding('0', 9) + l_padding('ElementalP', 22) + l_padding('elem_p', 12) + l_padding(str(round(item_w*10, 0)) +
                                                                                                                                                           '0000', 12) + l_padding('1.00000', 11) + l_padding('broadcast', 14) + ('   n   '*count_w + '   y   ') + ('   n   '*(alternatives_irr - count_w - 1)) + '\n'
        g.write(line_w)

    line11 = l_padding('rot_reset', 27) + l_padding('hru', 12) + l_padding('0', 12) + l_padding('reset_1', 21) + l_padding('null', 11) + \
        l_padding('1.00000', 11) + l_padding('0.00000', 11) + l_padding('null',
                                                                        14) + ('   n   '*(alternatives_irr-1) + '   y   ') + '\n'+'\n'
    g.write(line11)


# managament file
    mf3 = l_padding('agrl_{}'.format(h), 34) + \
        l_padding('0', 13) + l_padding('1', 21) + '\n'
    mf4 = r_padding('fh{0}'.format(h), 66) + '\n'
    mf.write(mf3), mf.write(mf4)
# landuse file
    with open(lum_file, 'a') as l_file:
        if agrl_line_1[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_1[1], 13) + l_padding(agrl_line_1[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_1[4], 17) + l_padding(agrl_line_1[5], 25) + l_padding(
                agrl_line_1[6], 18) + l_padding(agrl_line_1[7], 10) + l_padding(agrl_line_1[8], 26) + l_padding(agrl_line_1[9], 18) + l_padding(agrl_line_1[10], 18) + l_padding(agrl_line_1[11], 18) + l_padding(agrl_line_1[12], 18) + l_padding(agrl_line_1[13], 18) + '\n'
            l_file.write(lum)
            list_mgt_irr.append('agr_{0}'.format(h))
        elif agrl_line_2[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_2[1], 13) + l_padding(agrl_line_2[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_2[4], 17) + l_padding(agrl_line_2[5], 25) + l_padding(
                agrl_line_2[6], 18) + l_padding(agrl_line_2[7], 10) + l_padding(agrl_line_2[8], 26) + l_padding(agrl_line_2[9], 18) + l_padding(agrl_line_2[10], 18) + l_padding(agrl_line_2[11], 18) + l_padding(agrl_line_2[12], 18) + l_padding(agrl_line_2[13], 18) + '\n'
            l_file.write(lum)
            list_mgt_irr.append('agr_{0}'.format(h))
        elif agrl_line_3[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_3[1], 13) + l_padding(agrl_line_3[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_3[4], 17) + l_padding(agrl_line_3[5], 25) + l_padding(
                agrl_line_3[6], 18) + l_padding(agrl_line_3[7], 10) + l_padding(agrl_line_3[8], 26) + l_padding(agrl_line_3[9], 18) + l_padding(agrl_line_3[10], 18) + l_padding(agrl_line_3[11], 18) + l_padding(agrl_line_3[12], 18) + l_padding(agrl_line_3[13], 18) + '\n'
            l_file.write(lum)
            list_mgt_irr.append('agr_{0}'.format(h))
        elif agrl_line_4[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_4[1], 13) + l_padding(agrl_line_4[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_4[4], 17) + l_padding(agrl_line_4[5], 25) + l_padding(
                agrl_line_4[6], 18) + l_padding(agrl_line_4[7], 10) + l_padding(agrl_line_4[8], 26) + l_padding(agrl_line_4[9], 18) + l_padding(agrl_line_4[10], 18) + l_padding(agrl_line_4[11], 18) + l_padding(agrl_line_4[12], 18) + l_padding(agrl_line_4[13], 18) + '\n'
            l_file.write(lum)
            list_mgt_irr.append('agr_{0}'.format(h))
        elif agrl_line_5[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_5[1], 13) + l_padding(agrl_line_5[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_5[4], 17) + l_padding(agrl_line_5[5], 25) + l_padding(
                agrl_line_5[6], 18) + l_padding(agrl_line_5[7], 10) + l_padding(agrl_line_5[8], 26) + l_padding(agrl_line_5[9], 18) + l_padding(agrl_line_5[10], 18) + l_padding(agrl_line_5[11], 18) + l_padding(agrl_line_5[12], 18) + l_padding(agrl_line_5[13], 18) + '\n'
            l_file.write(lum)
            list_mgt_irr.append('agr_{0}'.format(h))
        elif agrl_line_6[0] == 'agrl_lum':
            lum = l_padding('agr_{0}'.format(h), 34) + l_padding(agrl_line_6[1], 13) + l_padding(agrl_line_6[2], 19) + l_padding('agrl_{}'.format(h), 16) + l_padding(agrl_line_6[4], 17) + l_padding(agrl_line_6[5], 25) + l_padding(
                agrl_line_6[6], 18) + l_padding(agrl_line_6[7], 10) + l_padding(agrl_line_6[8], 26) + l_padding(agrl_line_6[9], 18) + l_padding(agrl_line_6[10], 18) + l_padding(agrl_line_6[11], 18) + l_padding(agrl_line_6[12], 18) + l_padding(agrl_line_6[13], 18) + '\n'
            l_file.write(lum)
            list_mgt_irr.append('agr_{0}'.format(h))
g.close()
mf.close()


# creating the hru-data file
hru_dict = dict(zip(col_hru, list_mgt))
# subject to change
hru_data_file = '/{}/Scenarios/Default/TxtInOut/hru-data.hru'.format(datazone)
file_df = pd.read_csv(hru_data_file, delimiter='\s+', skiprows=1)
for k in hru_dict:
    file_df.loc[file_df['id'] == k, 'lu_mgt'] = hru_dict[k]
# irrigated hrus
hru_dict_irr = dict(zip(col_hru_irr, list_mgt_irr))
for x in hru_dict_irr:
    file_df.loc[file_df['id'] == x, 'lu_mgt'] = hru_dict_irr[x]
file_df = file_df.replace(np.nan, 'null')
# subject to change
draft_hru_data = '/{}/Scenarios/Default/management_tables/useless_data.hru'.format(
    datazone)
file_df.to_csv(draft_hru_data, sep='\t', index=False)

# subject to change
new_hru_data_file = '/{}/Scenarios/Default/management_tables/hru-data.hru'.format(
    datazone)
hru_txt = open(new_hru_data_file, 'w+')
with open(draft_hru_data, 'r') as f:
    for line in f:
        new_line = line.split()
        re_line = r_padding(new_line[0], 8) + r_padding(new_line[1], 18) + r_padding(new_line[2], 18) + r_padding(new_line[3], 18) + \
            r_padding(new_line[4], 18) + r_padding(new_line[5], 18) + \
            r_padding(new_line[6], 18) + r_padding(new_line[7], 18) + \
            r_padding(new_line[8], 18) + r_padding(new_line[9], 18) + '\n'
        hru_txt.write(re_line)
hru_txt.close()

with open(new_hru_data_file, 'r') as contents:
    save = contents.read()
with open(new_hru_data_file, 'w') as contents:
    contents.write(
        'hru-data.hru:-{0}- nkwasa 2021 {1}'.format(datazone,  now.strftime("%Y-%m-%d %H:%M:%S")) + '\n')
with open(new_hru_data_file, 'a') as contents:
    contents.write(save)


# creating the landuse.lum file

# subject to change
landuse_file = '/{}/Scenarios/Default/management_tables/landuse_copy.lum'.format(
    datazone)
landuse_df = pd.read_csv(landuse_file, delimiter='\s+', skiprows=1)
for k in community_plant:
    landuse_df.loc[landuse_df['name'] == k, 'plnt_com'] = community_plant[k]
landuse_df = landuse_df.replace(np.nan, 'null')
landuse_df.to_csv('landuse_draft.lum', sep='\t', index=False)
# subject to change
draft_landuse_file = '/{}/Scenarios/Default/management_tables/landuse_draft.lum'.format(
    datazone)
# subject to change
new_landuse_file = '/{}/Scenarios/Default/management_tables/landuse.lum'.format(
    datazone)
hru_txt = open(new_landuse_file, 'w+')
with open(draft_landuse_file, 'r') as f:
    for line in f:
        new_line = line.split()
        re_line = r_padding(new_line[0], 8) + r_padding(new_line[1], 18) + r_padding(new_line[2], 18) + r_padding(new_line[3], 18) + \
            r_padding(new_line[4], 18) + r_padding(new_line[5], 18) + \
            r_padding(new_line[6], 18) + r_padding(new_line[7], 18) + \
            r_padding(new_line[8], 18) + r_padding(new_line[9], 18) + r_padding(
                new_line[10], 18) + r_padding(new_line[11], 18) + r_padding(new_line[12], 18) + r_padding(new_line[13], 18) + '\n'
        hru_txt.write(re_line)
hru_txt.close()
with open(new_landuse_file, 'r') as contents:
    save = contents.read()
with open(new_landuse_file, 'w') as contents:
    contents.write(
        'landuse.lum:-{0}- nkwasa 2021 {1}'.format(datazone,  now.strftime("%Y-%m-%d %H:%M:%S")) + '\n')
with open(new_landuse_file, 'a') as contents:
    contents.write(save)
os.remove('landuse_copy.lum'), os.remove('landuse_draft.lum')


# creating the plant.ini file

# subject to change
plant_file_dir = '/{}/Scenarios/Default/management_tables/plant.ini'.format(
    datazone)
plnt_file = open(plant_file_dir, 'w+')
line100 = 'plant.ini:-{}- written by SWAT+ editor version nkwasa 2021 {}'.format(
    datazone,  now.strftime("%Y-%m-%d %H:%M:%S")) + '\n'
line200 = l_padding('pcom_name', 27) + l_padding('plt_cnt', 18) + l_padding('rot_yr_ini', 18) + l_padding('plt_name', 18) + l_padding('lc_status', 18) + l_padding(
    'lai_init', 18) + l_padding('bm_init', 18) + l_padding('phu_init', 18) + l_padding('plnt_pop', 18) + l_padding('yrs_init', 18) + l_padding('rsd_init', 18) + '\n'
line300 = l_padding('agrl_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line400 = r_padding('agrl', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line500 = l_padding('corn_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line600 = r_padding('corn', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line700 = l_padding('soyb_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line800 = r_padding('soyb', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line900 = l_padding('rice_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line1000 = r_padding('rice', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line1100 = l_padding('swht_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line1200 = r_padding('swht', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line1300 = l_padding('sugc_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line1400 = r_padding('sugc', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line1500 = l_padding('bana_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line1600 = r_padding('bana', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
line1700 = l_padding('frst_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line1800 = r_padding('frst', 62) + r_padding('y', 12) + r_padding('2.00000', 17) + r_padding('50000.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('1.00000', 19) + r_padding('10000.00000', 19) + '\n'
line1900 = l_padding('past_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line2000 = r_padding('past', 62) + r_padding('y', 12) + r_padding('0.00000', 17) + r_padding('20000.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('1.00000', 19) + r_padding('10000.00000', 19) + '\n'
line2100 = l_padding('rnge_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line2200 = r_padding('rnge', 62) + r_padding('y', 12) + r_padding('2.00000', 17) + r_padding('20000.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('1.00000', 19) + r_padding('10000.00000', 19) + '\n'
line2300 = l_padding('agrr_comm', 33) + l_padding('1', 21) + \
    l_padding('1', 27) + '\n'
line2400 = r_padding('agrr', 62) + r_padding('n', 12) + r_padding('0.00000', 17) + r_padding('0.00000', 17) + \
    r_padding('0.00000', 17) + r_padding('0.00000', 19) + \
    r_padding('0.00000', 19) + r_padding('10000.00000', 19) + '\n'
plnt_file.write(line100), plnt_file.write(line200), plnt_file.write(line300), plnt_file.write(line400), plnt_file.write(line500), plnt_file.write(line600), plnt_file.write(line700), plnt_file.write(line800), plnt_file.write(
    line900), plnt_file.write(line1000), plnt_file.write(line1100), plnt_file.write(line1200), plnt_file.write(line1300), plnt_file.write(line1400), plnt_file.write(line1500), plnt_file.write(line1600), plnt_file.write(line1700), plnt_file.write(line1800), plnt_file.write(line1900), plnt_file.write(line2000), plnt_file.write(line2100), plnt_file.write(line2200), plnt_file.write(line2300), plnt_file.write(line2400)
plnt_file.close()
os.remove('useless_data.hru')

# saving the crop journal just for reference and checking
df_crop_journal = pd.DataFrame(journal_crops.items(), columns=['hru', 'crop'])
df_crop_journal.to_csv('crop_journal_{}.csv'.format(
    datazone), sep=',', index=False)

with open('crop_journal_{}.con'.format(datazone), 'w') as f:
    f.write('crop journal: written by nkwasa 2021 on {} \n'.format(
        now.strftime("%Y-%m-%d %H:%M:%S")))
    df_crop_journal.to_string(f, col_space=5, index=False)


time_range = []
for k in os.listdir(crop_path):
    time_range.append(k.split()[0].split('.')[0][-4:])

print('default crops generated for {0} model for land use maps {1} --> {2}'.format(
    datazone, time_range[0], time_range[-1]))

print('\t >finished')
