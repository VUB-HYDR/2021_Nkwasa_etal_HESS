#!/usr/bin/env python3
'''this script is used to extract fertilizer maps from management.nc(LUH-ISIMIP) 
                          and P fertilizer from Lu Tian ascii files
                          
Author  : albert nkwasa
Contact : albert.nkwasa@vub.be
Date    : 2021.04.13

'''
import pandas as pd
import os
import xarray as xr
import rioxarray
import shutil
import rasterio
import glob
import numpy as np

# years for which fertilizer is required
year_1 = 2000
year_2 = 2010

# path to N fert storage
home = '/nitrogen_fertilizer'  # subject to change
os.chdir(home)

try:
    for k in os.listdir(home):
        if k.endswith('.tif'):
            os.remove(k)
except:
    pass


def del_path(file_name):
    'delete any existing folder'
    try:
        shutil.rmtree(file_name)
    except:
        pass


del_path('working_dir')


def new_path(file_name):
    try:
        os.mkdir(file_name)
    except:
        pass


new_path('working_dir')


wkdir = home + '\\working_dir'  # subject to change

os.chdir(wkdir)

# nitrogen fertilizer

# path to management file in LUH2-ISIMIP (Hurtt et al., 2020)
path_mgt = '/management.nc'
read_mgt = xr.open_dataset(path_mgt, decode_times=False)

# extracting variables from the netcdf file
variable_list = []
for f in read_mgt.variables:
    variable_list.append(f)
# delete the unwanted variables
unwanted_variables = {'lat', 'lon', 'time', 'lat_bounds', 'lon_bounds', 'combf', 'fulwd', 'rndwd', 'flood', 'fharv_c4per', 'fharv_c3per',
                      'crpbf_c3nfx', 'irrig_c3nfx', 'crpbf_c4per', 'irrig_c4per', 'crpbf_c3per', 'irrig_c3per', 'crpbf_c4ann', 'irrig_c4ann', 'crpbf_c3ann', 'irrig_c3ann'}
variable_list = [i for i in variable_list if i not in unwanted_variables]


# extraction of the years
# using a dummy of years and decoded_times
landuse_yr_path = '/landuse_years.csv'  # subject to change
df_landuse_yrs = pd.read_csv(landuse_yr_path)
df_slice = (df_landuse_yrs[(df_landuse_yrs['time_normal'] >= year_1) & (
    df_landuse_yrs['time_normal'] <= year_2)])
dec_time = df_slice['decode_time'].to_list()
nor_time = df_slice['time_normal'].to_list()
yr_interest = dict(zip(dec_time, nor_time))

# We iterate through the variable list and the time period to generate specific rasters
for k in variable_list:
    try:
        new_path(wkdir + '\\temp')  # subject to change
    except:
        pass
    path_files = wkdir + '\\temp'  # subject to change
    for l in yr_interest:
        parameter = read_mgt[k].sel(time=l)
        parameter.rio.set_spatial_dims('lon', 'lat', inplace=True)
        parameter.rio.set_crs('epsg:4326')
        fl_out = path_files + "\\{0}_{1}.tif".format(k, yr_interest[l])
        parameter.rio.to_raster(fl_out)

# getting the land use fraction occupied by the agric crops

# path to the landuse file in LUH2-ISIMIP (Hurtt et al., 2020)
path_static = '/states.nc'
file_read_static = xr.open_dataset(path_static, decode_times=False)
variables_static = []
for f in file_read_static.variables:
    variables_static.append(f)
# delete the unwanted variables
unwanted_variables_static = {'lat', 'lon', 'time',
                             'lat_bounds', 'lon_bounds', 'secma', 'secmb', 'pastr', 'primf', 'primn', 'urban', 'secdf', 'secdn', 'range'}
variables_static = [
    i for i in variables_static if i not in unwanted_variables_static]
for k in variables_static:
    static_files = wkdir + '\\temp'  # subject to change
    for h in yr_interest:
        parameter = file_read_static[k].sel(time=h)
        parameter.rio.set_spatial_dims('lon', 'lat', inplace=True)
        parameter.rio.set_crs('epsg:4326')  # subject to change
        fl_out = static_files + "\\{0}_{1}.tif".format(k, yr_interest[h])
        parameter.rio.to_raster(fl_out)


def read_tiff(file):
    with rasterio.open(file) as src:
        ras_data = src.read()
    return ras_data

# test_transition.tif is a presaved raster with the required meta data. # subject to change


def save_tiff(result, out_path):
    with rasterio.open('/test_transitions.tif') as src:
        ras_meta = src.profile
        ras_meta['nodata'] = -9999
    with rasterio.open(out_path, 'w', **ras_meta) as dst:
        dst.write(result)
    return out_path


def weight_lus_frst(dir_path, lu_c3ann, fert_c3ann, lu_c4ann, fert_c4ann, lu_c3per, fert_c3per, lu_c4per, fert_c4per, lu_c3nfx, fert_c3nfx, out_dir):
    os.chdir(dir_path)
    file_multiply = read_tiff(lu_c3ann)*read_tiff(fert_c3ann) + read_tiff(lu_c4ann)*read_tiff(fert_c4ann) + read_tiff(
        lu_c3per)*read_tiff(fert_c3per) + read_tiff(lu_c4per)*read_tiff(fert_c4per) + read_tiff(lu_c3nfx)*read_tiff(fert_c3nfx)
    file_sum = read_tiff(lu_c3ann) + read_tiff(lu_c4ann) + \
        read_tiff(lu_c3per) + read_tiff(lu_c4per) + read_tiff(lu_c3nfx)
    file_div = np.divide(file_multiply, file_sum, where=file_sum != 0)
    # converts kg/ha --> g/m2 (catered for in the mgt scripts)
    file_units = np.divide(file_div, 10, where=file_div != -9999)
    save_tiff(file_units, out_dir)


for i in nor_time:
    weight_lus_frst(wkdir + '\\temp', 'c3ann_{}.tif'.format(i), 'fertl_c3ann_{}.tif'.format(i), 'c4ann_{}.tif'.format(i), 'fertl_c4ann_{}.tif'.format(i), 'c3per_{}.tif'.format(i),
                    'fertl_c3per_{}.tif'.format(i), 'c4per_{}.tif'.format(i), 'fertl_c4per_{}.tif'.format(i), 'c3nfx_{}.tif'.format(i), 'fertl_c3nfx_{}.tif'.format(i), home + '\\nfery_agg_{}.tif'.format(i))

os.chdir(home)
del_path('working_dir')

# phosphorus fertilizer

phosp_dir = '/phosphorus_fertilizer'
os.chdir(phosp_dir)
try:
    for k in os.listdir('/phosphorus_fertilizer'):
        if k.endswith('.tif'):
            os.remove(k)
except:
    pass

# path to the phosphorus files (Zhang et al., 2017)
phop_path = '/zhang_phosphorus'
path_join = glob.glob(os.path.join(phop_path, '*.asc'))

for k in nor_time:
    for file in path_join:
        year = int(file.split('\\')[-1].split('.')[0][-4:])
        if year == k:
            output_file = phosp_dir + "\\pfery_agg_{}.tif".format(k)
            os.system(
                'gdal_translate -of GTiff -a_srs EPSG:4326 {} {}'.format(file, output_file))


print('fertilizer maps generated for model for a period {0} --> {1}'.format(
    nor_time[0], nor_time[-1]))
print('\t >finished')
