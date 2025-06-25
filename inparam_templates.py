import numpy as np 
import pickle

def inparam_source_file(eventfile):

# takes in the event (dict) and returns the list of sources and their moment tensors 
# formatted in accordance with the AxiSEM3D inparam.source.yaml file 


    with open(eventfile, 'rb') as f:
        event = pickle.load(f)

    distances, azimuths, depths, areas, slips, strikes, dips, rakes, slip_times = (
        event['distances'],
        event['azimuths'],
        event['depths'],
        event['patch_areas'],
        event['slips'],
        event['strikes'],
        event['dips'],
        event['rakes'],
        event['slip_times']
    )

    source_input_strings = []
    c  = 0
    
    for i in range(slips.shape[0]):
        for j in range(slips.shape[1]):
            
            moment_tensor =  patch_moment_tensor(strikes[i, j], dips[i, j], rakes[i, j], slips[i, j], areas[i, j], azimuths[i,j])
            source_input_string = f"""
        - point{c}:
            location:
                distance_azimuth: [{distances[i,j]:.3f}, {azimuths[i,j]:.3f}]
                depth: {depths[i,j]*1e3:.3e}
                ellipticity: false
                depth_below_solid_surface: true
                undulated_geometry: true
            mechanism:
                type: MOMENT_TENSOR
                data: [{moment_tensor[0]:.3e}, {moment_tensor[1]:.3e}, {moment_tensor[2]:.3e}, {moment_tensor[3]:.3e}, {moment_tensor[4]:.3e}, {moment_tensor[5]:.3e}]
                unit: 1
            source_time_function:
                class_name: GaussianSTF
                half_duration: 1.
                decay_factor: 1.628
                time_shift: {slip_times[i,j]+1:.3f}
                use_derivative_integral: ERF"""
                    
            c = c + 1       
            source_input_strings.append(source_input_string)
            

    input_string_1 = f"""#
#  inparam.source.yaml
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 6/20/20.
#  Copyright © 2020 Kuangdai Leng. All rights reserved.
#

#  This is an AxiSEM3D input parameter file for
#  time-axis and sources


############################# time axis #############################
# parameters for the time axis of the simulation
time_axis:
    # what: record length (the end time in seismograms)
    # type: double
    # note: the start time depends on the source-time functions
    record_length: 60.

    # what: user-specified Δt
    # type: string / doublectober
    # only: NONE / value
    # note: use NONE to automatically determine Δt by mesh
    enforced_dt: NONE

    # what: the Courant number for determining Δt by mesh
    # type: double
    # note: 1) Δt increases with the Courant number; decrease it when
    #          numerical instability occurs
    #       2) [safe] 0.5 <===> 1.0 [aggressive]; 0.6~0.7 normally works
    #       3) if Courant_number < 0.3 but instability still occurs,
    #          it is likely to be an issue caused by an input 3D model
    #          (e.g., mislocation near a model boundary)
    Courant_number: 0.6

    # what: time integrator
    # type: string
    # only: NEWMARK / SYMPLECTIC
    # note: 1) NEWMARK is faster while SYMPLECTIC is less dispersive
    #       2) use SYMPLECTIC for ultra-long simulations
    #       3) Δt can be larger for SYMPLECTIC than for NEWMARK
    integrator: NEWMARK
    

######################### sources #########################
# what: list of sources
# type: array of objects
# note: 1) multiple sources are allowed
#       2) use [] if no source presents
list_of_sources:
    #==========================================================================#
    # this key can be arbitrary"""
         
    input_strings = input_string_1 + "\n".join(source_input_strings)   
    
    #with open('inparam.source.yaml', 'w') as file:
        #file.write(input_strings)    
        
    return input_strings 



def patch_moment_tensor(strike, dip, rake, slip, area, azimuth):

    # this function computes the moment tensor from strike, rake, and dip, slip and area
    # the output moment tesnor follows the same convention as that in AxiSEM3D source file
    # Fatima, October 24,2023


    # first, we define the normal vector n and the slip vector s are in the Aki and Richards system (North, East, Down)
    # change the strike, dip and rake from degrees to radians

    strike = np.deg2rad(strike)
    dip    = np.deg2rad(dip)
    rake   = np.deg2rad(rake)

    cst = 3.2e+10
    M0  = slip * area * cst

    Mxx = -M0 * (np.sin(dip) * np.cos(rake) * np.sin(2 * strike) + np.sin(2 * dip) * np.sin(rake) * (np.sin(strike))**2)
    Mxy =  M0 * (np.sin(dip) * np.cos(rake) * np.cos(2 * strike) + 0.5 * np.sin(2 * dip) * np.sin(rake) * np.sin(2 * strike))
    Mxz = -M0 * (np.cos(dip) * np.cos(rake) * np.cos(strike) + np.cos(2 * dip) * np.sin(rake) * np.sin(strike))
    Myy =  M0 * (np.sin(dip) * np.cos(rake) * np.sin(2 * strike) - np.sin(2 * dip) * np.sin(rake) * (np.cos(strike))**2)
    Myz = -M0 * (np.cos(dip) * np.cos(rake) * np.sin(strike) - np.cos(2 * dip) * np.sin(rake) * np.cos(strike))
    Mzz =  M0 * np.sin(2 * dip) * np.sin(rake)
    
    ## write the moment tensor components for AxiSEM3D

    moment_tensor = np.array([[Mzz, Mxz, Myz],
                              [Mxz, Mxx, Mxy],
                              [Myz, Mxy, Myy]])

    
    # rotate the moment tensor by the azimuth (because the source coordinates are defined using distance-azimuth)
    
    R_x = np.array([[1, 0, 0],
                    [0, np.cos((np.pi)/4 - azimuth), -np.sin((np.pi)/4 - azimuth)],
                    [0, np.sin((np.pi)/4 - azimuth),  np.cos((np.pi)/4-   azimuth)]])

    rot_tensor = R_x @ moment_tensor @ R_x.T
    
    reduced_tensor = np.array([rot_tensor[0,0], rot_tensor[1,1], rot_tensor[2,2], rot_tensor[0,1], rot_tensor[0,2], rot_tensor[1,2]]).ravel()

    return reduced_tensor 

def inparam_model_file(exodus_filename, center_lat, center_lon, nc_filename, vel_coords):
    """
    Generates the inparam.model.yaml file for AxiSEM3D
    Input:
        exodus_filename: mesh file
        center_lat     : latitude of center of cylinder(numerical domain)
        center_lon     : longitude of center of cylinder(numerical domain)
        nc_filename    : velocity nc file
        vel_coords     : coordinate system of velocity model: choose between DISTANCE_AZIMUTH / XY_CARTESIAN / LATITUDE_LONGITUDE
    Output: 
    Fatima, December 19, 2023
    """

    if vel_coords == 'XY_CARTESIAN':
        nc_variables = '[x, y, depth]'
    elif vel_coords == 'DISTANCE_AZIMUTH':
        nc_variables = '[r, phi, depth]'
        data_rank = '[2, 1, 0]'
    elif vel_coords == 'LATITUDE_LONGITUDE':
        nc_variables = '[latitude, longitude, depth]'
        data_rank = '[0, 1, 2]'

    input_string_1 = f"""#
#  inparam.model.yaml
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 6/20/20.
#  Copyright © 2019 Kuangdai Leng. All rights reserved.
#

#  This is an AxiSEM3D input parameter file for
#  1D and 3D models


############################# 1D model #############################
# parameters for 1D model (the mesh)
model1D:
    # what: Exodus mesh file created by salvus mesher
    # type: filename
    exodus_mesh: {exodus_filename}


############################# geodesy #############################
# parameters for geodesy
geodesy:
    # what: geographic location of the north pole in the mesh
    # type: array of double / SOURCE
    # note: 1) this reference location enables the usage of geographic
    #          coordinates for locating sources, receivers and 3D models,
    #          compatible with Cartesian meshes
    #       2) array of double: [latitude, longitude]
    #       3) SOURCE: determined this location by the FIRST source
    #                  presented in list_of_sources in inparam.source.yaml;
    #                  always use SOURCE for a single-source simulation
    lat_lon_north_pole_mesh: [{center_lat:.2f}, {center_lon:.2f}]

    # what: flattening on the surface
    # type: string / double
    # only: SPHERE / WGS84 / GRS80 / SPECFEM3D_GLOBE / value
    # note: 1) ellipticity is ignored for a Cartesian mesh
    #       2) 0 for a perfect sphere; ~0.0033 for the Earth
    #       3) ellipticity will be used in the transformation between
    #          the geographic and the geocentric co-latitudes;
    #          see eq. (14.32) in Dahlen & Tromp, 1998
    #       4) to actually deform the entire mesh, add 3D model
    #          "Ellipticity" to list_of_3D_models
    flattening_on_surface: WGS84


######################## absorbing boundary ########################
# parameters for absorbing boundary condition
absorbing_boundary:
    # what: model boundaries regarded as absorbing boundaries
    # type: array of string
    # only: a subset of [RIGHT, BOTTOM, TOP]
    # note: 1) an AxiSEM3D mesh may contain four outer boundaries:
    #          left (axial), right, bottom and top (surface); the right,
    #          bottom and top ones can be absorbing boundaries (the left
    #          or axial one is non-physical)
    #       2) use [] to disable absorbing boundary condition
    #          (so that all model boundaries will be stress-free)
    #       3) the most common case in seismology is [RIGHT, BOTTOM]
    boundaries: [RIGHT, BOTTOM]

    # what: enable the Clayton-Enquist approach
    # type: bool
    # note: the simplest linear approach by Clayton & Engquist (1977)
    enable_Clayton_Enquist: true

    # the sponge approach by Kosloff & Kosloff (1986)
    Kosloff_Kosloff:
        # what: enable the Kosloff-Kosloff approach
        # type: bool
        # note: Clayton-Enquist and Kosloff-Kosloff can be used together,
        #       but one of them has to be enabled at least
        enable: true

        # what: relative spans of the sponge layers
        # type: array of double
        # note: 1) must be presented in the same order as absorbing_boundaries
        #       2) to use Kosloff-Kosloff, the mesh should be a little larger
        #          than the required computational domain; for example, given
        #          a required domain spans from 0 to 100 km in depth, one can
        #          generate a mesh from 0 to 110 km and set the relative span
        #          to 0.05, so the thickness of the sponge layer at the mesh
        #          bottom will be determined as 110 * 0.05 = 5.5 km, leaving
        #          an unaffected depth range from 0 to 104.5 km for normal
        #          wave propagation and analysis
        #          allowed range: .01 ~ 0.25
        relative_spans: [.1, .1]

        # what: expression of γ-factor in solid domain
        # type: math expression
        # note: 1) γ-factor represents the absorbing strength at a point
        #       2) allowed arguments include (case sensitive):
        #          - VP, VS: P- and S- wave velocities at the point
        #          - RHO   : density at the point
        #          - SPAN  : span of the sponge layer
        #          - T0    : mesh period
        #          * VP, VS and RHO are the 1D values in the Exodus mesh
        #       3) this expression will be further multiplied by a pattern
        #          function that equals to 1 on the outermost edge of the
        #          sponge layer (i.e., on the mesh boundary) and gradually
        #          decreases house outline to 0 on the the innermost edge; such a decreasing
        #          pattern is automatically handled by the solver
        #       4) the default is an empirical expression from
        #          Haindl et al., 2020
        gamma_expr_solid: 4.4 / T0 * (VS / VP)^2 * exp(-0.08 * SPAN / (VP * T0))

        # what: expression of γ-factor in fluid domain
        # type: math expression
        # note: same as gamma_expr_solid but without VS dependency
        gamma_expr_fluid: 1.76 / T0 * exp(-0.08 * SPAN / (VP * T0))


######################## attenuation ########################
# what: attenuation mode
# type: string
# only: NONE / FULL / CG4
# note: 1) NONE: turn off attenuation
#       2) FULL: compute attenuation on all GLL points
#       3) CG4:  compute attenuation on 4 GLL points per element;
#                CG4 is mostly as accurate as FULL but more efficient
#                than FULL, see van Driel & Nissen-​Meyer, 2014;
#                CG4 requires set(NPOL 4) in CMakeLists.txt;
attenuation: CG4

############################# 3D models #############################
# what: list of 3D models
# type: array of objects
# note: 1) the order in this list can affect the final 3D model
#       2) use [] if no 3D model presents
list_of_3D_models:
    #==========================================================================#
    # this key can be arbitrary """

    input_string_2 = []

    for i in range(len(nc_filename)):
        input_string_2.append(f"""
    - vol_layer{i+1}:
        # what: activate this model
        # type: bool
        activated: true
        # what: class name
        # type: string
        # note: current built-in classes include
        #       - StructuredGridV3D: volumetric 3D model on a structured grid
        #       - StructuredGridG3D: geometric 3D model on a structured grid
        #       - StructuredGridO3D: ocean-load 3D model on a structured grid
        #       - Ellipticity: deform the mesh with global ellipticity
        class_name: StructuredGridV3D
        # -------------------------------- #
        # parameters for StructuredGridV3D #
        # -------------------------------- #
        # what: NetCDF data file
        # type: filename
        nc_data_file: {nc_filename[i]}
        # parameters for grid coordinates
        coordinates:
            # what: type of horizontal coordinates
            # type: string
            # only: DISTANCE_AZIMUTH / XY_CARTESIAN / LATITUDE_LONGITUDE
            horizontal: DISTANCE_AZIMUTH
            # what: type of vertical coordinate
            # type: string
            # only: RADIUS / DEPTH
            vertical: DEPTH
            # what: correct for ellipticity when locating the model
            # type: bool
            # note: used only when horizontal = LATITUDE_LONGITUDE
            ellipticity: false
            # what: use solid surface as depth origin
            # type: bool
            # note: used only when vertical = DEPTH
            depth_below_solid_surface: true
            # what: NetCDF variables for the coordinates
            # type: array of string
            nc_variables: {nc_variables}
            # what: rank of the coordinates in data
            # type: array of int
            data_rank: {data_rank}
            # what: length unit of the coordinates
            # type: string / value
            # only: km / m / number
            length_unit: m
            # what: angle unit of the coordinates
            # type: string
            # only: degree / radian
            angle_unit: degree
            # what: use undulated (otherwise reference) geometry to
            #       determine the vertical location
            # type: bool
            # note: compatible only with vertical = RADIUS
            undulated_geometry: false
            # what: use element center for model scope check
            # type: bool
            # note: this feature allows for in-plane discontinuities
            whole_element_inplane: false
        # parameters for properties
        properties:
            - VP:
                # what: NetCDF variable
                nc_var: VP
                # what: factor or unit
                factor: 1
                # what: reference kind
                # only: ABS / REF1D / REF3D / REF_PERTURB
                reference_kind: ABS
            - VS:
                # what: NetCDF variable
                nc_var: VS
                # what: factor or unit
                factor: 1
                # what: reference kind
                # only: ABS / REF1D / REF3D / REF_PERTURB
                reference_kind: ABS
            - RHO:
                # what: NetCDF variable
                nc_var: RHO
                # what: factor or unit
                factor: 1
                # what: reference kind
                # only: ABS / REF1D / REF3D / REF_PERTURB
                reference_kind: ABS
            - QMU:
                # what: NetCDF variable
                nc_var: QMU
                # what: factor or unit
                factor: 1
                # what: reference kind
                # only: ABS / REF1D / REF3D / REF_PERTURB
                reference_kind: ABS
            - QKAPPA:
                # what: NetCDF variable
                nc_var: QKAPPA
                # what: factor or unit
                factor: 1
                # what: reference kind
                # only: ABS / REF1D / REF3D / REF_PERTURB
                reference_kind: ABS
        # what: store grid data only on the leader processors
        # type: bool
        # note: turn this on if the model is large; set mpi:nproc_per_group
        #       in inparam.advanced.yaml to the number of processors per
        #       node to minimize memory usage
        store_grid_only_on_leaders: true """)

    input_string = input_string_1 + "\n".join(input_string_2)

    return input_string

   


def inparam_output_file(stations_filename, stations_coords, dt):
    
# Input:
#       -stations_coords: only: LATITUDE_LONGITUDE / DISTANCE_AZIMUTH / XY_CARTESIAN
#       - dt: sampling period
    
    
    input_string = f"""#
#
#  inparam.output.yaml
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 6/20/20.
#  Copyright © 2020 Kuangdai Leng. All rights reserved.
#

#  This is an AxiSEM3D input parameter file for
#  station-wise and element-wise output


############################# station-wise #############################
# what: list of station groups
# type: array of objects
# note: 1) different options such as channels and sampling rates can
#          be used for different station groups; for example, one may
#          have one group of real-world seismic stations to record the
#          displacement vector at a high sampling rate and another group
#          of animation stations to record only the vertical displacement
#          at a low sampling rate
#       2) use [] if no station group presents
list_of_station_groups: 
    - fixed:
        # station locations
        locations:
            # what: station location file
            # type: filename
            # note: 1) must be an ascii file with five or six columns:
            #          NAME NETWORK x1 x2 [useless] x3
            #          [useless] is for compatibility with SPECFEM and AxiSEM
            #       2) (x1, x2, x3) are specified by the next two options
            #       3) empty lines and comment lines (led by #) are allowed
            station_file: {stations_filename}
            # what: horizontal coordinates x1 and x2
            # type: string
            # only: LATITUDE_LONGITUDE / DISTANCE_AZIMUTH / XY_CARTESIAN
            # note: 1) the unit for LATITUDE_LONGITUDE is degree
            #       2) the unit for DISTANCE is either radian (for a
            #          spherical mesh) or meter (for a Cartesian mesh);
            #          the unit for AZIMUTH is radian
            #       3) the unit for XY_CARTESIAN is meter
            horizontal_x1_x2: {stations_coords}
            # what: vertical coordinate x3
            # type: string
            # only: RADIUS / DEPTH
            # note: the unit is meter
            vertical_x3: DEPTH
            # what: correct for ellipticity when locating the stations
            # type: bool
            # note: used only when horizontal_x1_x2 = LATITUDE_LONGITUDE
            ellipticity: true
            # what: use solid surface as depth origin
            # type: bool
            # note: used only when vertical_x3 = DEPTH
            depth_below_solid_surface: true
            # what: use undulated (otherwise reference) geometry to
            #       determine the vertical locations
            # type: bool
            # note: compatible with vertical_x3 = both RADIUS and DEPTH
            undulated_geometry: true
        # wavefields to be recorded
        wavefields:
            # what: coordinate frame of wavefields
            # type: string
            # only: spz / RTZ / ENZ / xyz
            # note: 1) spz: (s, phi, z) or AxiSEM3D-intrinsic
            #       2) RTZ: (radial, transpose, vertical)
            #       3) ENZ: (east, north, vertical)
            #       4) xyz: (x, y, z) in source-centered frame
            coordinate_frame: ENZ
            # what: type of medium
            # type: string
            # only: SOLID / FLUID
            # note: all stations in a group must be located in either
            #       the solid or the fluid domain
            medium: SOLID
            # what: list of channels
            # type: array of string
            # note: 1) allowed channels for medium = SOLID
            #          * displacement:
            #            U, U1, U2, U3, U_NORM (or |U|)
            #          * gradient of displacement:
            #            G, G11, G12, G13, G21, G22, G23, G31, G32, G33,
            #            Gii (or G_I1)
            #          * strain:
            #            E, E11, E12, E13, E21, E22, E23, E31, E32, E33,
            #            Eii (or E_I1), E_J2
            #          * stress:
            #            S, S11, S12, S13, S21, S22, S23, S31, S32, S33,
            #            Sii (or S_I1), S_J2
            #          * curl:
            #            R, R1, R2, R3, R_NORM (or |R|)
            #      2) allowed channels for medium = FLUID
            #          * displacement:
            #            U, U1, U2, U3, U_NORM (or |U|)
            #          * scalar potential of displacement (U = ∇X):
            #            X
            #          * pressure:
            #            P
            #      3) (1, 2, 3) are determined by coordinate_frame
            #      4) using U means [U1, U2, U3], and similarly for G, E, S
            #         and R; duplicated channels are automatically removed
            channels: [U]
        # temporal sampling
        temporal:
            # what: sampling period
            # type: string / double
            # only: DT / DTx2 / DTx3 / ... / value
            # note: DT stands for Δt of the simulation; DTx3 means
            #       sampling period = Δt * 3
            sampling_period: {dt:.1f}
            # what: time window
            # type: string / array of double
            # only: FULL / [t0, t1]
            # note: use FULL to record the whole simulation
            time_window: FULL
        # file options
        file_options:
            # what: output file format
            # type: string
            # only: ASCII_STATION / ASCII_CHANNEL / NETCDF
            # note: 1) ASCII_STATION: one ascii file contains all channels at
            #                         one station, available only for a small
            #                         number of stations
            #       2) ASCII_CHANNEL: one ascii file contains one channel at
            #                         all stations, available for many stations
            #       3) NETCDF: much more efficient than ascii, available for
            #                  many stations; parallel NetCDF can be activated
            #                  in CMakeLists.txt
            format: NETCDF
            # what: number of sampled time steps to be buffered
            # type: int
            # note: 1) the solver buffers wave data during the time loop for
            #          efficient writing; increase this buffer size to save
            #          writing time and decrease it to save memroy (useful for
            #          a large number of stations)
            #       2) this parameter does not affect the final results
            buffer_size: 1000
            # what: flush file after writing a buffer to it
            # type: bool
            # note: 1) pro: minimizes data loss in case of abnormal termination
            #          con: hits output performace if buffer_size is small
            #       2) this parameter does not affect the final results
            flush: true
            
############################# element-wise #############################
# what: list of element groups
# type: array of objects
# note: 1) different options such as channels and sampling rates can
#          be used for different element groups
#       2) use [] if no element group presents
list_of_element_groups: []
            """
    return input_string

    
def inparam_advanced_file():

    input_string = f"""#
#
#  inparam.advanced.yaml
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 6/20/20.
#  Copyright © 2019 Kuangdai Leng. All rights reserved.
#

#  This is an AxiSEM3D input parameter file for
#  advanced settings


##################### verbosity #####################
# parameters for verbosity
verbose:
    # what: verbose to stdout or a file
    # type: string
    # only: STDOUT / filename
    channel: STDOUT

    # what: verbosity level
    # type: string
    # only: NONE / ESSENTIAL / DETAILED
    level: ESSENTIAL

    # what: show/hide runtime warnings
    # type: bool
    warnings: true

    # what: time step interval to display time loop info
    # type: int
    loop_info_interval: 1000

    # what: time step interval for stability check
    # type: int
    # note: use 1 to precisely locate the origin of instability
    stability_interval: 1


##################### mpi #####################
# parameters for mpi
mpi:
    # what: number of processors in a mpi group
    # type: int
    # note: 1) AxiSEM3D uses a two-level MPI structure where the
    #          processors are divided into groups to avoid broadcasting
    #          a large input dataset (e.g., the Exodus mesh or a large
    #          3D model) on every processor; instead, a large dataset can
    #          be stored only on a "leader" processor in each group, which
    #          handles data requests from its members
    #       2) increase this number (from 1 to the number of processors per
    #          per node) to save memory
    nproc_per_group: 1
    
    # what: weight for load balancing
    # type: string
    # only: ELEMENT / ELEMENT_POINT
    # note: 1) ELEMENT:       use cost measurement on elements
    #       2) ELEMENT_POINT: use cost measurement on both elements and points
    weight_for_load_balancing: ELEMENT_POINT

    # what: plot domain decomposition
    # type: bool
    # note: the output netcdf file contains three variables:
    #       * coords,   double, (X, 2), (s,z) of the element centers
    #       * mpi_rank, int,    (X, ),  mpi rank of the elements
    #       * weights,  double, (X, ),  element weights for decomposition
    #       where X is the number of elements
    plot_domain_decomposition: false


##################### developers #####################
# parameters for developers
develop:
    # what: enable/disable preloop diagnosis
    # type: bool
    # note: 1) output/develop/preloop_diagnosis.log for runtime and memory
    #       2) output/develop/cost_measurements.log for cost measurements
    #          on elements and GLL points
    diagnose_preloop: true

    # what: maximum time steps for running
    # type: int
    # note: use 0 to free this limit
    max_num_time_steps: 0

    # what: wall-clock time limit (sec) for FFTW planning
    # type: double
    time_limit_for_fftw_planning: 60.

    # what: enforce FFTW lucky numbers
    # type: bool
    # note: FFTW is good at handling logical sizes of the form:
    #       n = (2^a)*(3^b)*(5^c)*(7^d)*(11^e)*(13^f), where e + f < 2,
    #       as called the lucky numbers; users should use true.
    fftw_lucky_numbers: true"""
    
    return input_string

    
def inparam_nr_file(nr_type, nr):
    
#Input: nr_type: CONSTANT / ANALYTICAL / POINTWISE / STRUCTURED
#            nr: resolution in the azimuthal direction

    input_string = f"""#
#
#  inparam.nr.yaml
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 6/20/20.
#  Copyright © 2019 Kuangdai Leng. All rights reserved.
#

#  This is an AxiSEM3D input parameter file for
#  the Nr field: Nr = Nr(s,z)
#  s:  horizontal coordinate in the 2D mesh
#  z:  vertical coordinate in the 2D mesh
#  Nr: the azimuthal resolution of solution, i.e., the number of
#      "cardinal points" along the "ring" generated by rotating the
#      point (s,z) around the axis by 2π; Nu = Nr/2 is the Fourier
#      expansion order of solution at (s,z)


##################### type #####################
# what: type of Nr(s,z)
# type: string
# only: CONSTANT / ANALYTICAL / POINTWISE / STRUCTURED
# note: 1) CONSTANT:   Nr(s,z) = const
#       2) ANALYTICAL: analytical Nr(s,z) defined in NrFieldAnalytical.cpp
#       3) POINTWISE:  Nr provided at discrete control points
#       4) STRUCTURED: Nr provided on a structured grid
type_Nr: {nr_type}

# what: bound Nr(s,z) from above by inplane resolution
# type: bool
# note: there is no reason to use an azimuthal resolution higher than
#       the inplane (or mesh) resolution; users should use true.
bound_Nr_by_inplane: true


##################### constant #####################
# what: the constant value for type_Nr = CONSTANT
# type: int
# note: for an axisymmetric model with a single axial source, use
#       1) 5 for a moment tensor (earthquake)
#       2) 3 for a force vector (impact)
#       3) 1 for a pressure (explosion) in either solid or fluid
constant: {nr}


##################### analytical #####################
# parameters for type_Nr = ANALYTICAL
analytical:
    # what: code ID to match NrFieldAnalytical::sCodeID
    # type: string
    # note: to ensure that AxiSEM3D has been compiled with the wanted
    #       NrFieldAnalytical.cpp, repeat here the code ID defined by
    #       NrFieldAnalytical::sCodeID in NrFieldAnalytical.cpp (line 18)
    code_ID: depth-dependent (AxiSEM3D default)

    # what: parameters used by the default NrFieldAnalytical.cpp
    # note: 1) the default NrFieldAnalytical.cpp implements a
    #          depth-dependent Nr(s,z), i.e., Nr(s,z) = Nr(depth),
    #          with code ID = "depth-dependent (AxiSEM3D default)"
    #       2) linear interpolation is applied between two control depths
    depth_dependent_AxiSEM3D_default:
        # what: the control depths
        # type: array of double
        control_depths: [0., 50e3, 100e3, 6371e3]

        # what: Nr at the control depths
        # type: array of double
        Nr_at_control_depths: [100, 100, 50, 50]

    # what: any user-defined parameters for NrFieldAnalytical.cpp
    # type: any
    # note: these parameters can have arbitrary names and types, depending
    #       on how they are read and used in NrFieldAnalytical.cpp
    any_user_defined_parameters:
        example__bool: true
        example__string: Hello world!
        example__array_of_double: [1., 2., 3.]
        example__array_of_string: [path, file1, file2]


##################### pointwise #####################
# parameters for type_Nr = POINTWISE
pointwise:
    # what: netcdf data file
    # type: filename
    # note: 1) this file must contain the following two variables:
    #          * pointwise_sz, double, (X, 2), (s,z) of X control points
    #          * pointwise_Nr, int,    (X, ),  Nr at the X control points
    #       2) the unit is meter for s and z
    #       3) interpolation is based on inverse distance weighting
    #       4) another variable starting_Nr_for_scanning will exist if
    #          this file has been created by wavefield scanning
    nc_data_file: SFBA_finite_rupture_Nr0.nc

    # what: factor multiplied to Nr(s,z)
    # type: double
    # note: useful if nc_data_file was created by wavefield scanning;
    #       for example, Nr(s,z) obtained by scanning s20rts may be
    #       applied to s40rts by using a factor of 2.0
    multip_factor: 1.0


##################### structured #####################
# parameters for type_Nr = STRUCTURED
structured:
    # what: netcdf data file
    # type: filename
    # note: 1) for a Cartesian mesh, this file must contain three variables:
    #          * structured_s,  double, (M, ),  s of M grid points
    #          * structured_z,  double, (N, ),  z of N grid points
    #          * structured_Nr, int,    (M, N), Nr at the M*N grid points
    #       2) for a sphercial mesh, replace _s with _r and _z with _t (theta)
    #       3) the unit is meter for s, z and r and radian for theta
    nc_data_file: structured.nc

    # what: value of Nr at any location out of the grid range
    # type: int
    value_out_of_range: 5


##################### wavefield scanning #####################
# ~~~~~~~~~~~~~~~
# Q: What is wavefield scanning?
# A: AxiSEM3D "learns" a sub-optimal Nr(s,z) during a simulation and
#    dump the resultant Nr(s,z) into a file, which can be re-used in
#    subsequent simulations with similar input parameters, such as
#    similar (or simpler) 3D models, source depth and record length.
# ~~~~~~~~~~~~~~~
# Q: How wavefield scanning works?
# A: Starting from the current Nr(s,z), AxiSEM3D analyzes the Fourier
#    series of the wavefield at any mesh point upon the arrival of an
#    energy peak (measured in H2-norm) and determines if any "small"
#    higher-order terms can be truncated away; at the end of the
#    simulation, the required Nr (maximum over time) is dumped to file.
# ~~~~~~~~~~~~~~~
# Q: How to enable wavefield scanning?
# A: Setting enable_scanning = true. The starting (current) Nr(s,z) can be
#    any of CONSTANT, ANALYTICAL, POINTWISE, STRUCTURED. The resultant
#    Nr(s,z) will be no greater than the starting Nr(s,z).
# ~~~~~~~~~~~~~~~
# Q: How to re-use Nr(s,z) obtained by scanning?
# A: Use type_Nr = POINTWISE, setting nc_data_file to the output file
#    of wavefield scanning (after moving it from output/ to input/).
# ~~~~~~~~~~~~~~~

# parameters for wavefield scanning
wavefield_scanning:
    # what: enable/disable wavefield scanning
    # type: bool
    # note: enabling wavefield scanning barely slows a simulation but
    #       will increase memory usage
    enable_scanning: false

    # what: output file
    # type: filename
    output_file: NZ_point_source_Nr_learned.nc

    # what: relative threshold for the convergence of Fourier series
    # type: double
    # note: 1) this parameter represents the accuracy loss by truncating
    #          the Fourier series of the wavefield
    #       2) allowed range: [1e-4, 1e-1]
    threshold_Fourier_convergence: 1e-2
    
    # what: relative amplitude skipped for scanning
    # type: double
    # note: 1) an energy peak with an amplitude smaller than
    #          "this relative amplitude * the largest energy peak"
    #          will be skipped for scanning
    #       2) using 1. means that the resultant Nr accounts only for
    #          the largest energy peak across the record length
    #       3) using 0. means that the resultant Nr accounts for all
    #          the energy peaks across the record length
    #       4) allowed range: [0., 1.]
    relative_amplitude_skipped: 0.

    # advanced scanning parameters (users are unlikely to change)
    advanced:
        # what: absolute amplitude skipped for scanning
        # type: double
        # note: 1) tiny values must be skipped to avoid numerical errors
        #       2) allowed range: [1e-14, 1e-10]
        absolute_amplitude_skipped: 1e-14

        # what: maximum number of energy peaks to be recorded
        # type: int
        # note: use a small one to consider only a few largest peaks
        max_num_peaks: 10

        # what: perform scanning only on vertex GLL points
        # type: bool
        # note: vertex-only scanninng can significantly decrease both
        #       runtime memory and output file size
        vertex_only: false

        # what: how many time steps per mesh period to detect energy peaks
        # type: int
        # note: must be no less than 4; recommended range: [8, 16]
        num_steps_per_mesh_period: 12"""
    
    return input_string



