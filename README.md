# 2021_Nkwasa_et_al

Scripts used in the analysis of Nkwasa et al (submitted).

## To install
Python scripts should be run on python3.

## For users
This repository includes the processing  scripts used in Nkwasa et al, (submitted). Python is used to extract data from global datasets and generate decision tables for agricultural HRUs in the SWAT+ model 
for large scale hydrological applications.
For the python code:  
 (a) 'create_crop_maps.py' --> generates crop ids for respective crops that occupy the largest agricultural pixels from the LUH2-ISIMIP map(s) (Hurtt et al., 2020).
 (b) 'extract_fertilizer.py' --> creates fertilizer maps for the required simulation period using the ISIMIP N fertilizer files (Hurtt et al., 2020)  and the P fertilizer files from (Lu and Tian., 2017) 
 (c) 'global_managment.py' --> creates the management files using the global datasets and the model TxtInOut folder files. 


## Versions
Version 0.1.0 - April 2021 

## License
This project is licensed under the MIT License.

