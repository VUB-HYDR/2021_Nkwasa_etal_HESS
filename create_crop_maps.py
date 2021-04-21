#!/usr/bin/env python
'''this script creates a reclassified agrl land use map for crops - ISIMIP based on the crop
that occupies the largest percentage of a pixel

Author  : albert nkwasa
Contact : albert.nkwasa@vub.be
Date    : 2021.04.13

'''

import os
import xarray as xr
import rioxarray as rio
import rasterio
import numpy as np
import pandas as pd
import shutil


# specify the year (s) for which the crop map(s) are required
year_1 = 2000  # subject to change
year_2 = 2015  # subject to change

# path for storage of the generated map
os.chdir('/crop_maps')  # subject to change
wkdir = os.getcwd()


def del_path(file_name):
    'delete any existing folder'
    try:
        shutil.rmtree(file_name)
    except:
        pass


del_path(wkdir + '\\reclassified_landuse')


def new_path(file_name):
    try:
        os.mkdir(file_name)
    except:
        pass


new_path(wkdir + '\\reclassified_landuse')


# extraction of the years
# using a dummy of years and decoded_times
landuse_yr_path = '/landuse_years.csv'  # subject to change
df_landuse_yrs = pd.read_csv(landuse_yr_path)
df_slice = (df_landuse_yrs[(df_landuse_yrs['time_normal'] >= year_1) & (
    df_landuse_yrs['time_normal'] <= year_2)])
dec_time = df_slice['decode_time'].to_list()
nor_time = df_slice['time_normal'].to_list()
yr_interest = dict(zip(dec_time, nor_time))


# creation of the landuse tiffs
# path to states file from LUH2 - ISIMIP
landuse_file = '/states.nc'  # subject to change
file = xr.open_dataset(landuse_file, decode_times=False)
variable = []
for i in file.variables:
    variable.append(i)

drop_variables = {'time', 'lat', 'lon', 'lat_bounds',
                  'lon_bounds', 'secmb', 'secma', 'primf', 'primn', 'secdf', 'secdn', 'urban', 'pastr', 'range'}
variable_list = [i for i in variable if i not in drop_variables]

for k in variable_list:
    for l in yr_interest:
        parameter = file[k].sel(time=l)
        parameter.rio.set_spatial_dims('lon', 'lat', inplace=True)
        parameter.rio.set_crs('epsg:4326')
        fl_out = wkdir + \
            '\\reclassified_landuse\\{0}_{1}.tif'.format(k, yr_interest[l])
        parameter.rio.to_raster(fl_out)


os.chdir(wkdir + '\\reclassified_landuse')
for i in yr_interest:
    for k in os.listdir():
        if k.split()[0].split('.')[0][6:10] == str(yr_interest[i]):
            if k.startswith('c3ann'):
                with rasterio.open(k) as src:
                    src_c3ann = src.read(1)
                    src_meta = src.profile
                    src_c3ann = np.where(np.isnan(src_c3ann), -9999, src_c3ann)
            if k.startswith('c3nfx'):
                with rasterio.open(k) as src:
                    src_c3fnx = src.read(1)
                    src_c3fnx = np.where(np.isnan(src_c3fnx), -9999, src_c3fnx)
            if k.startswith('c3per'):
                with rasterio.open(k) as src:
                    src_c3per = src.read(1)
                    src_c3per = np.where(np.isnan(src_c3per), -9999, src_c3per)
            if k.startswith('c4per'):
                with rasterio.open(k) as src:
                    src_c4per = src.read(1)
                    src_c4per = np.where(np.isnan(src_c4per), -9999, src_c4per)
            if k.startswith('c4ann'):
                with rasterio.open(k) as src:
                    src_c4ann = src.read(1)
                    src_c4ann = np.where(np.isnan(src_c4ann), -9999, src_c4ann)

            # reclassifying agrl land use
    c3ann_1 = np.where((src_c3ann > src_c3fnx) & (src_c3ann > src_c3per) & (
        src_c3ann > src_c4ann) & (src_c3ann > src_c4per), 1, src_c3ann)
    c4ann_2 = np.where((src_c4ann > src_c3ann) & (src_c4ann > src_c3fnx) & (
        src_c4ann > src_c3per) & (src_c4ann > src_c4per), 2, src_c4ann)
    c3nfx_3 = np.where((src_c3fnx > src_c3ann) & (src_c3fnx > src_c3per) & (
        src_c3fnx > src_c4ann) & (src_c3fnx > src_c4per), 3, src_c3fnx)
    c3per_4 = np.where((src_c3per > src_c3ann) & (src_c3per > src_c3fnx) & (
        src_c3per > src_c4ann) & (src_c3per > src_c4per), 4, src_c3per)
    c4per_5 = np.where((src_c4per > src_c3ann) & (src_c4per > src_c3fnx) & (
        src_c4per > src_c3per) & (src_c4per > src_c4ann), 5, src_c4per)

    # joining/merging the reclassified maps
    step_1 = np.where((c4ann_2 != -9999) &
                      (c4ann_2 != 2), c3ann_1, c4ann_2)
    step_2 = np.where((c3nfx_3 != -9999) & (c3nfx_3 != 3), step_1, c3nfx_3)
    step_3 = np.where((c3per_4 != -9999) & (c3per_4 != 4), step_2, c3per_4)
    step_4 = np.where((c4per_5 != -9999) & (c4per_5 != 5), step_3, c4per_5)
    step_5 = np.where(step_4 == 0, 6, step_4)
    # Takes care of inconsistences
    step_6 = np.where((step_5 != 1) & (step_5 != 2) & (step_5 != 3) & (
        step_5 != 4) & (step_5 != 5) & (step_5 != 6) & (step_5 != -9999), 2, step_5)
    # reading the meta data
    src_meta['dtype'] = step_6.dtype
    src_meta['nodata'] = -9999

    agrl_raster = 'reclassified_landuse_{}.tif'.format(yr_interest[i])
    with rasterio.open(agrl_raster, 'w', **src_meta) as dst:
        dst.write(step_6, 1)


for k in os.listdir():
    if k.startswith('c3ann') or k.startswith('c4ann') or k.startswith('c3per') or k.startswith('c4per') or k.startswith('c3nfx'):
        os.remove(k)

print('\t > finished')
