#######################################################################
#
# NDR sensitivity analysis example
#
# 1/ Runs NDR a bunch of times, varying input parameters
# 2/ Compares #1 to observed data and creates output tables
#       with difference/error values
#
# We're modifying both biophysical coefficients and
#   single-value UI parameters - each coefficient/parameter
#   has its own section below. Only one parameter is changed at
#   a time, and individual coefficients for all LULC classes are 
#   changed at once 
#
# For re-using this script, look for the comments
#   # MODIFY THESE # or # MODIFY THIS # and do so to customize
#   the paths to inputs/outputs.
#
# Run this script with python -u, or you won't get the status
#   information printed out while running, it will all cache
#   until the very end. Haven't figured out how to fix this
#   programmatically and don't feel like spending time on it.
#
# All of this could surely be done more efficiently,
#   but still a useful example.
#
# Original script by Rich Sharp, modified by Stacie Wolny
#
#
#######################################################################

import os
import pygeoprocessing
import numpy
import natcap.invest.ndr.ndr
import shutil
import csv
# Pandas provides easier csv/table data processing - very useful here
import pandas
import time
# Needed for post-processing NDR results
import osgeo
from osgeo import gdal
from osgeo import ogr
import glob
import shutil
import sys


#                                                                  #
# MODIFY THESE - INPUT/OUTPUT FOLDERS, NODATA AND OBSERVED VALUE   #
#                                                                  #

# Folder where input data lives
DATA_BASE_DIR = (
    r'H:\Restoration\Data\Model_Inputs\NDR\NDR_inputs_baseline_sensitivity_expand')
# Folder where NDR outputs will be written
WORKSPACE_DIR = (
    r'H:\Restoration\Data\Model_Outputs\NDR\Sensitivity_test1_p')
# Folder where biophysical tables will be written that contain updated
#   calibration values for each run - just for a sanity/debugging check
ADJUSTED_TABLE_DIR = (
    r'H:\Restoration\Data\Model_Outputs\NDR\Sensitivity_test1_p\adjust_coeff_tables')
# Folder where post-processing output will be written
#   - Multiplying by effective export, reaggregating, etc
POST_PROC_DIR = (
    r'H:\Restoration\Data\Model_Outputs\NDR\Sensitivity_test1_p\post_process')\

# Output file name of the table where final difference/error results will be written (CSV)
OUTPUT_RESULT_TABLE_PATH = (os.path.join(POST_PROC_DIR, 'sensitivity_results_NDRp.csv'))

# Effective export raster output from Rich's dam retention script
##EFFECTIVE_EXPORT_RASTER_PATH = (
##    r'D:\Stacie\US\EPA_SNEP\data\Inputs_Moshassuck_clip\nutrient_eff_export_1fill.tif')

# Set NoData value for post-processing outputs - needed for pygeoprocessing calls
NODATA_VALUE = -1

# Observed nutrient export value that we're comparing our results to
#   for calculating error
# P value is annual average for the Catawba River station, 2008-2017
# Sites CatawbaR, LinvilleR, ElkCrk
SITE = "CatawbaR"
OBSERVED_VALUE_P = 4520.509


###############################################################################

def main():

    # Keep track of how long it takes to run the script
    start_time = time.time()
    

    #                                           #
    # MODIFY THESE - BASE DATA INPUTS TO NDR    #
    #                                           #
    args = {
        'biophysical_table_path': os.path.join(
            DATA_BASE_DIR, 'NDR_Biophysical.csv'),
        'calc_n': False,
        'calc_p': True,
        'dem_path': os.path.join(
            DATA_BASE_DIR, 'USGS_10m_snap.tif'),
        'flow_dir_algorithm': 'D8',
        'k_param': '2',
        'lulc_path': os.path.join(
            DATA_BASE_DIR, '2024_NLCD_expandclip_fixed.tif'),
        'runoff_proxy_path': os.path.join(
            DATA_BASE_DIR, 'prism_ppt_30yr_800m_ClipedExpandSite.tif'),
        'watersheds_path': os.path.join(
            DATA_BASE_DIR, f"MS_Drainage_{SITE}.shp"),
        'threshold_flow_accumulation': '500',
        'workspace_dir': WORKSPACE_DIR
    }
    
    # Create output folders if they don't exist
    if not os.path.exists(WORKSPACE_DIR):
        os.makedirs(WORKSPACE_DIR)
        print("Created output folder " + str(WORKSPACE_DIR))
    if not os.path.exists(POST_PROC_DIR):
        os.makedirs(POST_PROC_DIR)
        print("Created output folder " + str(POST_PROC_DIR))
    if not os.path.exists(ADJUSTED_TABLE_DIR):
        os.makedirs(ADJUSTED_TABLE_DIR)
        print("Created output folder " + str(ADJUSTED_TABLE_DIR))


    # Open final output sensitivity table
    #   - will be written to as each NDR run is done and post-processed
    results_file = open(OUTPUT_RESULT_TABLE_PATH, 'w')
    results_file.write("ws_id, SITE, param, param_value, modeled_p_export, observed_p, percent_error\n")

    #############################################################################
    
    ### Tweaking biophysical table parameters, one at a time
    
    #############################################################################
    
    
    # Read in original biophysical table
    orig_biophys_table = os.path.join(
        DATA_BASE_DIR, 'NDR_Biophysical.csv')
    # Copy the original table into the Workspace to work with subsequently, just in case
    shutil.copy(orig_biophys_table, WORKSPACE_DIR)
    # MODIFY THIS #
    orig_table_path = os.path.join(WORKSPACE_DIR, 'NDR_Biophysical.csv')   

    

#########################################################################
# run the changes to biophysical table for each variable individually
#########################################################################

    print("\n********************************************************")
    print("*** Processing load_p_perclass")
    print("********************************************************")

    # Multipliers applied to load_p, eff_p
    #   (.5x default value) -> (2x default value), in increments of .25x
    multiplier_sample_space = numpy.linspace(0.5, 0.5, num=1, endpoint=True)

    # Try pandas for table/data manipulation
    df_original = pandas.read_csv(orig_table_path)

    lucodes = df_original['lucode'].unique()

    # Multipliers are applied across all LULC classes individually
    for lucode in lucodes:
        if lucode != 82:
            continue
        for mult in multiplier_sample_space:
            df_biophys = df_original.copy()

            # Change biophysical table values for load_n
            df_biophys.loc[df_biophys['lucode'] == lucode, 'load_p'] *= mult

            # Write updated table to a new CSV, which will be used to run NDR
            mult_str = str(mult).replace('.', '-')
            new_biophys_table_path = (
                os.path.join(
                    ADJUSTED_TABLE_DIR, f'new_biophys_{SITE}_lucode_{lucode}_load_p_{mult_str}.csv'))
            df_biophys.to_csv(new_biophys_table_path, index=False)
           
            # Run NDR #
            run_workspace = os.path.join(WORKSPACE_DIR, f'workspace_{lucode}_load_p_{mult_str}')
            os.makedirs(run_workspace, exist_ok=True)
            local_args = args.copy()
            local_args['results_suffix'] = f"{SITE}_lucode_{lucode}_load_p_{mult_str}"
            local_args['biophysical_table_path'] = new_biophys_table_path
            local_args['workspace_dir'] = run_workspace

            print(f"\nRunning NDR with {lucode}_load_p * {mult}")
            natcap.invest.ndr.ndr.execute(local_args)
            print(f"Done running NDR with {lucode}_load_p * {mult}")


            print("Post-processing...")

            P_EXPORT_SHAPEFILE_PATH = os.path.join(run_workspace, 'watershed_results_ndr_' + str(local_args['results_suffix']) + '.gpkg')

            # Extract aggregated data from the aggregate output shapefile for P
            with ogr.Open(P_EXPORT_SHAPEFILE_PATH) as watershed_vector:
                watershed_layer = watershed_vector.GetLayer()
                watershed_feature = watershed_layer.GetNextFeature()
                ws_id = watershed_feature.GetField('Id')
                modeled_p = float(watershed_feature.GetField('p_surface_export'))
                error_p = ((modeled_p - OBSERVED_VALUE_P) / float(OBSERVED_VALUE_P)) * 100
                print(modeled_p)

            # Write data to the output sensitivity table
            results_file.write('%s,%s,%s,%f,%f,%f,%f\n' % (ws_id, SITE, f"{lucode}_load_p", mult, modeled_p, OBSERVED_VALUE_P, error_p))

        
            print("Done post-processing\n")
        
        print("Done post-processing biophysical\n")


    ###############################################################################

    #
    # Finish things up
    #

    # Close output sensivity table file
    results_file.close()
    
    # How long did it take to run all of this?
    elapsed_time = (time.time() - start_time)/60
    print("\nelapsed time: " + str(elapsed_time) + " minutes")
    
#####################################################################################


if __name__ == '__main__':
    main()

#####################################################################################
