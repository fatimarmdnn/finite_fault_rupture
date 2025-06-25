"""Input for AxiSEM3D
"""

import numpy as np
import argparse
from inparam_templates import *
import os
import shutil


def get_args():
    
    parser = argparse.ArgumentParser(description="Input for AxiSEM3D Simulations")
    
    ### AxiSEM3D Input Parameters ###
    parser.add_argument("--mesh_filename", type=str, help="exodus mesh filename", default = 'AxiSEMCartesian_sfba_m500_1s.e')
    parser.add_argument("--center_lat", type=float, help="latitude of the center of simulation cylinder",  default = 37.7)
    parser.add_argument("--center_lon", type=float, help="longitude of the center of simulation cylinder", default = -122.1)
    parser.add_argument("--vel_filenames", nargs='+', type=str, help="velocity filenames (.nc)", default= ['SFBA_layer1_min500.nc', 'SFBA_layer2.nc', 'SFBA_layer3.nc', 'SFBA_layer4.nc'])
    parser.add_argument("--vel_coords", type=str, choices=["DISTANCE_AZIMUTH", "XY_CARTESIAN", "LATITUDE_LONGITUDE"], default = 'DISTANCE_AZIMUTH')
    parser.add_argument("--eventfilename", type=str, help="event filename containing slip distribution data", required = True)
    parser.add_argument("--half_duration", type=float, help="half duration of the source (in seconds)", default = 2)
    parser.add_argument("--stations_name", type=str, help="example: SFBA_local", default = 'uniform')
    parser.add_argument("--stations_filename", type=str, help="stations filename (.txt)", default = 'STATIONS_OUTPUT.txt')
    parser.add_argument("--stations_coords", type=str, choices=["DISTANCE_AZIMUTH", "XY_CARTESIAN", "LATITUDE_LONGITUDE"], default = 'XY_CARTESIAN')
    parser.add_argument("--sampling_period", type=float, help="sampling period(in seconds)", default = 0.2)
    parser.add_argument("--nr_type", type=str, choices=["CONSTANT", "ANALYTICAL", "POINTWISE", "STRUCTURED"], default = 'POINTWISE')
    parser.add_argument("--nr", type=float, help = "nr value", default = 1000)
    parser.add_argument("--shared_input_dir", type=str, help="directory holding shared input files", default = 'input_files/')


    args = parser.parse_args()
    return args


def main():
    
    args = get_args() 
    
    sim_dir = "Hayward_Rupture"
    os.mkdir(sim_dir)
                                   
    ## Create Input Folder ##
    
    sim_input_dir = os.path.join(sim_dir, "input")

    # Create the input directory
    os.mkdir(sim_input_dir)
    
    with open(os.path.join(sim_input_dir, 'inparam.model.yaml'), 'w') as file:
        file.write(inparam_model_file(args.mesh_filename, args.center_lat, args.center_lon, args.vel_filenames, args.vel_coords))
        
    with open(os.path.join(sim_input_dir, 'inparam.source.yaml'), 'w') as file:
        file.write(inparam_source_file(args.eventfilename))
            
    with open(os.path.join(sim_input_dir, 'inparam.output.yaml'), 'w') as file:
        file.write(inparam_output_file(args.stations_filename, args.stations_coords, args.sampling_period))
        
    with open(os.path.join(sim_input_dir, 'inparam.nr.yaml'), 'w') as file:
        file.write(inparam_nr_file(args.nr_type, args.nr))
        
    with open(os.path.join(sim_input_dir, 'inparam.advanced.yaml'), 'w') as file:
        file.write(inparam_advanced_file())
    
    
    shutil.copy(os.path.join(args.shared_input_dir, 'AxiSEMCartesian_sfba_m500_1s.e'), sim_input_dir)
    shutil.copy(os.path.join(args.shared_input_dir, 'SFBA_finite_rupture_Nr0.nc'), sim_input_dir)
    shutil.copy(os.path.join(args.shared_input_dir, 'STATIONS_OUTPUT.txt'), sim_input_dir)
    # shutil.copy(os.path.join(args.shared_input_dir, 'SFBA_layer1_min500.nc'), sim_input_dir)
    # shutil.copy(os.path.join(args.shared_input_dir, 'SFBA_layer2.nc'), sim_input_dir)
    # shutil.copy(os.path.join(args.shared_input_dir, 'SFBA_layer3.nc'), sim_input_dir)
    # shutil.copy(os.path.join(args.shared_input_dir, 'SFBA_layer4.nc'), sim_input_dir)
    # shutil.copy(os.path.join(args.shared_input_dir, 'moho_smooth_24km.nc'), sim_input_dir)
    # shutil.copy(os.path.join(args.shared_input_dir, 'topo_smooth_100m_upsampled.nc'), sim_input_dir)


if __name__ == "__main__":
    main()

