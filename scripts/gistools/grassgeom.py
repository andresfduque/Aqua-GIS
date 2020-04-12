#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created by AndresD at 8/06/19.

Module that integrates best tools for geomorphometric analysis available in different GIS open source software

Features:
    + Raster Hydrological Correction
    + Watersheds Analysis
    + Watershed Discharges
    + Watershed Outlets
    + R Basin Preprocess
    + R Basin Simplified
    + R Basin
    + HAND Pointer
    + HAND Stream Profile
    + HAND Synthetic Inundation
    + HAND Specific Inundation
    + HAND Hydraulic Parameters

Pre-requisites:
    + SAGA GIS (compiled with python enabled)
    + WhiteBox Geospatial Data Analysis (pip installed)
    + GRASS python scripting library

@author:    Andres Felipe Duque Perez
Email:      aduquep@ingenevo.com.co

Credits:

Lindsay, J. B. (2016). Whitebox GAT: A case study in geomorphometric analysis. Computers & Geosciences, 95, 75-84.
http://dx.doi.org/10.1016/j.cageo.2016.07.003

Neteler, M., Bowman, M.H., Landa, M., Metz, M., 2012. GRASS GIS: A multi-purpose open source GIS. Environ Model Soft 31,
124–130.
https://doi.org/10.1016/j.envsoft.2011.11.014

Conrad, O., Bechtel, B., Bock, M., Dietrich, H., Fischer, E., Gerlitz, L., Wehberg, J., Wichmann, V., and Böhner, J.
(2015): System for Automated Geoscientific Analyses (SAGA) v. 2.1.4, Geosci. Model Dev., 8, 1991-2007,
doi:10.5194/gmd-8-1991-2015.

TODO: F2PY for extreme time consumption routines

"""
# modules import
import os
import sys
import math
import shutil
import sqlite3
from time import time
from subprocess import call

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import whitebox

import grass.script as grass
import grass.script.array as garray
from grass.exceptions import CalledModuleError


def check_r_bsn_requisites():
    """Check if all complements to run r_basin are installed."""
    found_missing = False
    for program in ('r.stream.basins', 'r.stream.distance', 'r.stream.extract', 'r.stream.order', 'r.stream.snap',
                    'r.stream.stats'):
        if not grass.find_program(program, '--help'):
            found_missing = True
            print("'%s' required. Please install '%s' first using 'g.extension %s'" % program, program, program)
    if found_missing:
        grass.fatal("An ERROR occurred running r.basin")


def f_geodesic_dist_wgs84(lat1, lat2, lon1, lon2):
    """
    Vicenty geodesic distance.

    Parameters
    ----------
    lat1: float
        latitude of point 1 [in degrees]
    lon1: float
        longitude of point 1 [in degrees]
    lat2: float
        latitude of point 2 [in degrees]
    lon2: float
        longitude of point 2 [in degrees]

    Returns
    -------
    Geodesic distance [meters]

    """
    # WGS84 ellipsoid parameters
    major_axis = 6378137.0  # Semi-mayor axis in meters [WGS84]
    minor_axis = 6356752.314245  # Semi-minor axis in meters [WGS84]
    flattening = (major_axis - minor_axis) / major_axis  # Flattening

    fill_float = 9.9692099683868690e+36

    # calculate distance
    delta_l = (lon2 - lon1) * math.pi / 180.  # Difference in longitude

    lambda1 = delta_l             # First approximation
    lambda_p = 2 * math.pi  # Convergence condition [lambda1-lambdaP]

    u_1 = math.atan((1 - flattening) * math.tan(lat1 * math.pi / 180))
    u_2 = math.atan((1 - flattening) * math.tan(lat2 * math.pi / 180))
    sin_u1 = math.sin(u_1)
    sin_u2 = math.sin(u_2)
    cos_u1 = math.cos(u_1)
    cos_u2 = math.cos(u_2)

    iter_limit = 100
    cos_sq_alpha = None

    while (abs(lambda1 - lambda_p) > 1e-12) & (iter_limit > 0):
        sin_lambda = math.sin(lambda1)
        cos_lambda = math.cos(lambda1)

        sin_sigma = math.sqrt((cos_u2 * sin_lambda) * (cos_u2 * sin_lambda) +
                              (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda) *
                              (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda))

        if sin_sigma == 0:
            print('co-incident points')
            sys.exit()

        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda
        sigma = math.atan2(sin_sigma, cos_sigma)

        sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma
        cos_sq_alpha = 1 - sin_alpha * sin_alpha

        cos2_sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alpha

        if np.isnan(cos2_sigma_m):  # equatorial line: cosSqAlpha=0
            cos2_sigma_m = 0

        c_1 = flattening / 16 * cos_sq_alpha * (4 + flattening * (4 - 3 * cos_sq_alpha))
        lambda_p = lambda1
        lambda1 = delta_l + (1 - c_1) * flattening * sin_alpha * \
            (sigma + c_1 * sin_sigma * (cos2_sigma_m + c_1 * cos_sigma * (-1 + 2 * cos2_sigma_m * cos2_sigma_m)))

        iter_limit = iter_limit - 1

        if iter_limit == 0:
            print('Formula failed to converge')
            distance = fill_float
            sys.exit()

    u_sq = cos_sq_alpha * (major_axis * major_axis - minor_axis * minor_axis) / (minor_axis * minor_axis)
    major_axis_delta = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    minor_axis_delta = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

    delta_sigma = minor_axis_delta * sin_sigma * \
        (cos2_sigma_m + minor_axis_delta / 4 *
         (cos_sigma * (-1 + 2 * cos2_sigma_m * cos2_sigma_m) - minor_axis_delta / 6 * cos2_sigma_m *
          (-3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos2_sigma_m * cos2_sigma_m)))

    distance = minor_axis * major_axis_delta * (sigma - delta_sigma)
    distance = np.int(distance * 1000.0)  # round to 1mm precision
    distance = distance / 1e3

    return distance


def f_geodesic_bearing(lat1, lat2, lon1, lon2):
    """
    Haversine geodesic bearing for geographic coordinates.

    source: http://www.movable-type.co.uk/scripts/latlong.html

    Parameters
    ----------
    lat1: float
        latitude of point 1 [in degrees]
    lon1: float
        longitude of point 1 [in degrees]
    lat2: float
        latitude of point 2 [in degrees]
    lon2: float
        longitude of point 2 [in degrees]

    Returns
    -------
    Geodesic bearing [degrees]

    """
    rho1 = lat1 * math.pi / 180
    rho2 = lat2 * math.pi / 180

    sigma1 = lon1 * math.pi / 180
    sigma2 = lon2 * math.pi / 180
    delta_sigma = sigma2 - sigma1

    delta_y = math.sin(sigma2 - sigma1) * math.cos(rho2)
    delta_x = math.cos(rho1) * math.sin(rho2) - math.sin(rho1) * math.cos(rho2) * math.cos(delta_sigma)

    brn = math.atan2(delta_y, delta_x)
    brn = brn * 180 / math.pi

    return brn


def f_bsn_morphometry(raster_names_df, vector_names_df, fields_names_df, out_dir, analysis_type='Full',
                      deleted_maps=True):
    """Morphometric analysis of the hydrologic basin, defined by closing point.

    Parameters
    ----------
    points_v : [type]
        [description]
    point_id_field : [type]
        [description]
    bsn_name_field : [type]
        [description]
    sub_bsn_name_field : [type]
        [description]
    raster_names_df : [type]
        [description]
    vector_names_df : [type]
        [description]
    out_dir : [type]
        [description]
    analysis_type : str, optional
        [description], by default 'Full'
    deleted_maps : bool, optional
        [description], by default True

    Returns
    -------
    morphometric_df: pandas.DataFrame
        Morphometric characterization for each basin.

    """
    # ==================================================================================================================
    # -- Main imports and preparations
    # ==================================================================================================================
    print('START BASIN MORPHOMETRIC CHARACTERIZATION')
    # check dependencies
    check_r_bsn_requisites()

    # define location type
    if grass.locn_is_latlong():
        coords_system = 'geo'
    else:
        coords_system = 'proj'

    # check if destination folder exists, if not, create it
    if not os.path.isdir(out_dir):
        folder = out_dir.split('/')
        parent_folder = "/".join(folder[:-2])
        basin_folder = "/".join(folder[:-1])
        if os.path.isdir(parent_folder):
            if not os.path.isdir(basin_folder):
                os.mkdir(basin_folder)
                os.mkdir(out_dir)
        else:
            print('!Can not create folder in the specified location!')

    # global raster variables
    r_slope = raster_names_df['slope_r']
    r_strahler = raster_names_df['strahler_network_r']
    r_elevation_global = raster_names_df['dem_original_r']
    r_drainage_global = raster_names_df['flow_direction_r']
    r_stream_network = raster_names_df['stream_network_thinned_r']
    r_accumulation_global = raster_names_df['flow_accumulation_r']
    r_downstream_distance_global = raster_names_df['downstream_distance_r']

    # global vector variables
    v_outlet_snap = 'outlet_snap'
    v_stream_network = vector_names_df['stream_network_v']
    v_hydrologic_points = vector_names_df['analysis_points_v']

    # global vector fields
    point_id_field = fields_names_df['bsn_id']
    bsn_name_field = fields_names_df['bsn_name']
    sub_bsn_name_field = fields_names_df['sub_bsn_name']

    # save current region
    grass.read_command('g.region', flags='p', save='original', overwrite=True)

    # ==================================================================================================================
    # -- Input points preparation
    # ==================================================================================================================
    # snap outlet to stream network [hardcoded to four times raster resolution]
    pixel_res = grass.raster_info(r_elevation_global)['ewres']
    try:
        grass.read_command('g.findfile', element='vector', file=v_outlet_snap)
        grass.run_command('v.db.connect', map=v_outlet_snap, table=v_hydrologic_points, flags='d')
        grass.run_command('r.stream.snap', input=v_hydrologic_points, output=v_outlet_snap,
                          stream_rast=r_stream_network, radius=2, overwrite=True, quiet=True)
        grass.run_command('v.db.connect', map=v_outlet_snap, table=v_hydrologic_points, flags='o')
    except grass.CalledModuleError:
        grass.run_command('r.stream.snap', input=v_hydrologic_points, output=v_outlet_snap,
                          stream_rast=r_stream_network, radius=2, overwrite=True, quiet=True)
        grass.run_command('v.db.connect', map=v_outlet_snap, table=v_hydrologic_points, flags='o')

    # set conditions for analysis type
    sql_path = grass.read_command('db.databases', driver='sqlite').replace('\n', '')
    con = sqlite3.connect(sql_path)
    sql_stat = 'SELECT * FROM %s ' % v_hydrologic_points
    analysis_points_df = pd.read_sql_query(sql_stat, con)
    bsn_id_list = analysis_points_df[point_id_field]
    bsn_name_list = analysis_points_df[bsn_name_field]
    sub_bsn_name_list = analysis_points_df[sub_bsn_name_field]

    # ==================================================================================================================
    # -- Create database to store morphometric parameters estimations
    # ==================================================================================================================
    parameter_index = ['1. Easting Centroid of basin',
                       '2. Northing Centroid of basin',
                       '3. Rectangle containing basin N',
                       '4. Rectangle containing basin W',
                       '5. Rectangle containing basin S',
                       '6. Rectangle containing basin E',
                       '7. Basin Area [km^2]',
                       '8. Basin Perimeter [km]',
                       '9. Basin Max Length [km]',
                       '10. Basin Width [km]',
                       '11. Basin Distance to Center of Gravity [km]',
                       '12. Basin Max Elevation [masl]',
                       '13. Basin Min Elevation [masl]',
                       '14. Basin Elevation Difference [m]',
                       '15. Basin Mean Elevation [masl]',
                       '16. Basin Mean Slope [m/m]',
                       '17. Basin Length of Directing Vector [km]',
                       '18. Basin Prevalent Orientation [degree from north, counterclockwise]',
                       '19. Mainchannel Length [km]',
                       '20. Mainchannel Mean Slope [m/m]',
                       '21. Mainchannel Mean Slope USGS [m/m]',
                       '22. Mainchannel Max Elevation [masl]',
                       '23. Stream Network Mean Hillslope Length [m]',
                       '24. Stream Network Magnitude O [und]',
                       '25. Stream Network Max Order [Strahler]',
                       '26. Stream Network Number of Streams [und]',
                       '27. Stream Network Total Stream Length [km]',
                       '28. Stream Network First Order Stream Frequency',
                       '29. Stream Network Drainage Density [km/km^2]',
                       '30. Stream Network Stream Frequency [num/km^2]',
                       '31. Stream Network Bifurcation Ratio [Horton]',
                       '32. Stream Network Length Ratio [Horton]',
                       '33. Stream Network Area ratio [Horton]',
                       '34. Stream Network Slope ratio [Horton]',
                       '35. Index Compactness Coefficient',
                       '36. Index Circularity Ratio',
                       '37. Index Topological Diameter',
                       '38. Index Elongation Ratio',
                       '39. Index Shape Factor',
                       '40. Index Massiveness coefficient',
                       '41. Index Orographic coefficient',
                       '42. Index Stability coefficient',
                       '43. Index Asymmetry']

    parameters_df = pd.DataFrame(data=None, index=parameter_index)

    # ==================================================================================================================
    # -- Morphometric characterization of each point [loop]
    # ==================================================================================================================
    for i in enumerate(analysis_points_df):
        basin_processing_time = time()

        # id, name
        bsn_id = bsn_id_list[i]
        bsn_name = bsn_name_list[i]
        sub_bsn_name = sub_bsn_name_list[i]

        # variable names for each basin analysis
        r_mask = 'r_mask'
        r_bsn = sub_bsn_name + '_basin'
        r_basin_slope = sub_bsn_name + '_slope'
        r_basin_stream = sub_bsn_name + '_stream'
        r_basin_strahler = sub_bsn_name + '_strahler'
        r_basin_drainage = sub_bsn_name + '_drainage'
        r_basin_distance = sub_bsn_name + '_dist2out'
        r_basin_elevation = sub_bsn_name + '_elevation'
        r_basin_half_basin = sub_bsn_name + '_half_basin'
        r_basin_accumulation = sub_bsn_name + '_accumulation'
        r_basin_slope_average = sub_bsn_name + '_slope_average'
        r_basin_height_average = sub_bsn_name + '_height_average'
        r_downstream_distance = sub_bsn_name + '_downstream_distance'
        r_basin_hillslope_distance = sub_bsn_name + '_hillslope_distance'

        v_basin = sub_bsn_name + '_basin'
        v_half_basin = sub_bsn_name + '_half_basin'
        v_mainchannel = sub_bsn_name + '_mainchannel'
        v_centroid = sub_bsn_name + '_centroid'
        v_mainchannel_split = sub_bsn_name + '_mainchannel_split'
        v_basin_stream_network = sub_bsn_name + '_stream_network'

        # ==============================================================================================================
        # -- Basin Delineation
        # ==============================================================================================================
        try:
            # extract analysis point
            point_id = bsn_id[i]
            grass.run_command('v.extract', input=v_hydrologic_points, output='Pour_Point%s' % point_id,
                              type='point', where='cat=%s' % point_id, overwrite=True, quiet=True)
            grass.run_command('v.to.rast', input='Pour_Point%s' % point_id, output='Pour_Point%s' % point_id,
                              use='cat', type='point', overwrite=True, quiet=True)
            pour_point_coordinates = grass.read_command('v.info', map='Pour_Point%s' % point_id, flags='g')
            pour_point_coordinates = dict(x.split('=', 1) for x in pour_point_coordinates.split('\n') if '=' in x)

            # delineate basin
            grass.run_command('r.stream.basins', direction=r_drainage_global, basins=r_bsn,
                              points='Pour_Point%s' % point_id, overwrite=True, quiet=True)

            print('!Delineation of %s basin DONE!' % sub_bsn_name)
        except TypeError:
            print('!Delineation of %s basin FAILED!' % sub_bsn_name)

        # ==============================================================================================================
        # -- Basin Vectorization
        # ==============================================================================================================
        try:
            grass.run_command('r.to.vect', input=r_bsn, output=v_basin, type='area', flags='sv', overwrite=True,
                              quiet=True)

            # add two columns to the table: area and perimeter
            grass.run_command('v.db.addcolumn', map=v_basin, columns='area double precision', quiet=True)
            grass.run_command('v.db.addcolumn', map=v_basin, columns='perimeter double precision', quiet=True)

            # populate perimeter column
            grass.run_command('v.to.db', map=v_basin, option='perimeter', units='kilometers', columns='perimeter',
                              quiet=True)
            # read perimeter
            tmp = grass.read_command('v.to.db', map=v_basin, option='perimeter', units='kilometers',
                                     columns='perimeter', flags='p')
            perimeter_basin = float(tmp.split('\n')[1].split('|')[1])

            # populate area column
            grass.run_command('v.to.db', map=v_basin, option='area', columns='area', units='kilometers', quiet=True)

            # read area
            tmp = grass.read_command('v.to.db', map=v_basin, option='area', units='kilometers', columns='area',
                                     flags='p')
            area_basin = float(tmp.split('\n')[1].split('|')[1])

            print('!Vectorization of %s basin DONE!' % sub_bsn_name)
        except TypeError:
            print('!Vectorization of %s basin FAILED!' % sub_bsn_name)

        # ==============================================================================================================
        # -- Mask and Cropping
        # ==============================================================================================================
        if analysis_type == 'Full':
            try:
                grass.run_command('r.mask', raster=r_bsn)

                # add mask
                expr = '%s = %s / %s' % (r_mask, r_bsn, r_bsn)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop accumulation map
                expr = '%s = %s' % (r_basin_accumulation, r_accumulation_global)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop elevation map
                expr = '%s = %s' % (r_basin_elevation, r_elevation_global)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop flow directions map
                expr = '%s = %s' % (r_basin_drainage, r_drainage_global)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop stream network map
                expr = '%s = %s' % (r_basin_stream, r_stream_network)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop slope map
                expr = '%s = %s' % (r_basin_slope, r_slope)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop stream network strahler map
                expr = '%s = %s' % (r_basin_strahler, r_strahler)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # crop stream network strahler vector
                grass.run_command('v.overlay', ainput=v_stream_network, binput=v_basin, operator='and', atype='line',
                                  output=v_basin_stream_network, olayer='0,1,0', overwrite=True, quiet=True)

                # crop downstream distance and adjust
                expr = '%s = %s' % (r_downstream_distance, r_downstream_distance_global)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)
                global_downstream_distance = grass.raster_info(r_downstream_distance)['min']

                expr = '%s = %s - %0.4f' % (r_downstream_distance, r_downstream_distance, global_downstream_distance)
                grass.run_command('r.mapcalc', expression=expr, overwrite=True)

                # half basin
                max_accumulation = grass.raster_info(r_basin_accumulation)['max']
                grass.run_command('r.watershed', elevation=r_basin_elevation, threshold=int(max_accumulation * 0.9),
                                  half_basin=r_basin_half_basin, flags='s', overwrite=True, quiet=True)
                grass.run_command('r.to.vect', input=r_basin_half_basin, output=v_half_basin, type='area', flags='sv',
                                  overwrite=True, quiet=True)
                grass.run_command('v.db.addcolumn', map=v_half_basin, columns='area double precision', quiet=True)
                grass.run_command('v.to.db', map=v_half_basin, option='area', columns='area', units='kilometers',
                                  quiet=True)

                sql_stat2 = 'SELECT * FROM %s ' % v_half_basin
                half_basin_df = pd.read_sql_query(sql_stat2, con)
                asymmetry_coefficient = half_basin_df['area'].max() / half_basin_df['area'].min()

                print('!Making and Cropping of %s basin DONE!' % sub_bsn_name)
            except TypeError:
                print('!Making and Cropping of %s basin FAILED!' % sub_bsn_name)
            # remove mask
            grass.run_command('r.mask', flags='r')

            # ==========================================================================================================
            # -- Morphometric characterization
            # ==========================================================================================================
            try:
                # distance to outlet
                start_time = time()
                grass.run_command('r.stream.distance', stream_rast='Pour_Point%s' % point_id, overwrite=True,
                                  direction=r_basin_drainage, flags='o', distance=r_basin_distance, quiet=True)
                print('!Downstream Distance to Outlet of %s basin CALCULATED!' % sub_bsn_name)
                print('Processing time = %.2f minutes' % ((time() - start_time)/60))

                start_time = time()
                # hill-slope distance to river network
                grass.run_command("r.stream.distance", stream_rast=r_basin_stream, direction=r_basin_drainage,
                                  elevation=r_basin_elevation, distance=r_basin_hillslope_distance, overwrite=True,
                                  quiet=True)
                print('!Hillslope Distance to River Network of %s basin CALCULATED!' % sub_bsn_name)
                print('Processing time = %.2f minutes' % ((time() - start_time)/60))

                # mean elevation
                grass.run_command("r.stats.zonal", base=r_bsn, cover=r_basin_elevation, method="average",
                                  output=r_basin_height_average, quiet=True)
                mean_elev = grass.raster_info(r_basin_height_average)['min']
                grass.run_command('g.remove', type='raster', name=r_basin_height_average, flags='f')
                print('!Mean Elevation of %s basin CALCULATED [%s msnm]!' % (sub_bsn_name, mean_elev))

                # mean slope
                grass.run_command("r.stats.zonal", base=r_bsn, cover=r_basin_slope, method="average",
                                  output=r_basin_slope_average)
                mean_slope = grass.raster_info(r_basin_slope_average)['min']
                mean_slope = mean_slope / 100   # m/m
                grass.run_command('g.remove', type='raster', name=r_basin_slope_average, flags='f')
                print('!Mean Slope of %s basin CALCULATED [%.2f m/m]!' % (sub_bsn_name, mean_slope))

                # centroid and mean slope
                baricenter_slope = grass.read_command("r.volume", input=r_slope, clump=r_bsn)
                baricenter_slope = baricenter_slope.split()
                # mean_slope = float(baricenter_slope[30].decode())
                basin_centroid_east = float(baricenter_slope[33])
                basin_centroid_north = float(baricenter_slope[34])
                print('!Centroid of %s basin CALCULATED [%.2fE, %.2fN]!' % (sub_bsn_name, basin_centroid_east,
                                                                            basin_centroid_north))
                # rectangle coordinates
                info_region_basin = grass.read_command("g.region", vect=v_basin, flags='m')
                dict_region_basin = dict(x.split('=', 1) for x in info_region_basin.split('\n') if '=' in x)
                print('!Rectangle Containing %s basin CALCULATED!' % sub_bsn_name)

                # directing vector
                if coords_system == 'proj':
                    delta_x = abs(float(basin_centroid_east) - float(pour_point_coordinates['east']))
                    delta_y = abs(float(basin_centroid_north) - float(pour_point_coordinates['north']))
                    length_orienting_vector = math.sqrt((delta_x ** 2) + (delta_y ** 2)) / 1000
                    print('!Directing Vector of %s basin CALCULATED!' % sub_bsn_name)
                else:
                    length_orienting_vector = f_geodesic_dist_wgs84(float(basin_centroid_north),
                                                                    float(pour_point_coordinates['north']),
                                                                    float(basin_centroid_east),
                                                                    float(pour_point_coordinates['east'])) / 1000
                    print('!Directing Vector of %s basin CALCULATED!' % sub_bsn_name)

                # prevalent orientation
                if coords_system == 'proj':
                    if delta_y != 0:
                        prevalent_orientation = math.atan(abs(delta_x) / abs(delta_y))
                        prevalent_orientation = prevalent_orientation * 180 / math.pi
                        if delta_x >= 0 > delta_y:
                            prevalent_orientation = 180 - prevalent_orientation
                        elif delta_x < 0 and delta_y < 0:
                            prevalent_orientation = 180 + prevalent_orientation
                        elif delta_x < 0 <= delta_y:
                            prevalent_orientation = 360 - prevalent_orientation
                    else:
                        if delta_x >= 0:
                            prevalent_orientation = 90
                        else:
                            prevalent_orientation = 270
                    print('!Prevalent Orientation of %s basin CALCULATED [%.2f degrees]!' % (sub_bsn_name,
                                                                                             prevalent_orientation))
                else:
                    prevalent_orientation = f_geodesic_bearing(float(basin_centroid_north),
                                                               float(pour_point_coordinates['north']),
                                                               float(basin_centroid_east),
                                                               float(pour_point_coordinates['east']))
                    print('!Prevalent Orientation of %s basin CALCULATED [%.2f]!' % (sub_bsn_name,
                                                                                     prevalent_orientation))

                # compactness coefficient [AFDP]
                compactness_coefficient = 0.28 * perimeter_basin / math.sqrt(area_basin)
                print('!Compactness Coefficient of %s basin CALCULATED [%.2f]!' % (sub_bsn_name,
                                                                                   compactness_coefficient))

                # circularity ratio
                circularity_ratio = (4 * math.pi * area_basin) / (perimeter_basin ** 2)
                print('!Circularity Ratio of %s basin CALCULATED [%.2f]!' % (sub_bsn_name, circularity_ratio))

                # main channel length and slope
                sql_stat_hack = 'SELECT * FROM %s ' % v_basin_stream_network
                hack_df = pd.read_sql_query(sql_stat_hack, con)
                hack_order = hack_df['hack'].min()
                expr = 'hack=%s' % hack_order
                grass.run_command('v.extract', input=v_basin_stream_network, output=v_mainchannel_split, type='line',
                                  where=expr, overwrite=True, quiet=True)
                grass.run_command("v.build.polylines", input=v_mainchannel_split, output=v_mainchannel, type='line',
                                  cats='first', overwrite=True, quiet=True)
                grass.run_command('v.to.db', map=v_mainchannel, option='length', units='kilometers', columns='length',
                                  quiet=True)

                sql_stat1 = 'SELECT * FROM %s ' % v_mainchannel_split
                stream_network_df = pd.read_sql_query(sql_stat1, con)
                stream_length = stream_network_df['length']
                stream_slope = stream_network_df['gradient']
                mainchannel_out_elevation = stream_network_df['outlet_elev']
                mainchannel_source_elevation = stream_network_df['source_elev']
                weighted_slope = stream_length * stream_slope / stream_length.sum()
                mainchannel_length = stream_length.sum() / 1000     # kilometers
                mainchannel_slope = weighted_slope.sum()            # m/m
                mainchannel_max_elevation = mainchannel_source_elevation.max()
                mainchannel_outlet_elevation = mainchannel_out_elevation.min()
                print('!Mainchannel Length of %s basin CALCULATED [%.2f km]!' % (sub_bsn_name, mainchannel_length))
                print('!Mainchannel Slope of %s basin CALCULATED [%.2f percent]!' % (sub_bsn_name,
                                                                                     mainchannel_slope * 100))

                # main channel mean slope (USGS: https://pubs.usgs.gov/sir/2006/5312/pdf/sir2006-5312.pdf)
                # n_85 = 0
                # n_10 = 0
                # while stream_network_df['cum_length'][n_85] > 0.15 * mainchannel_length * 1000:
                #     e_85 = stream_network_df['source_elev'][n_85]
                #     l_85 = stream_network_df['cum_length'][n_85]
                #     n_85 += 1
                # while stream_network_df['cum_length'][n_10] > 0.90 * mainchannel_length * 1000:
                #     e_10 = stream_network_df['source_elev'][n_10]
                #     l_10 = stream_network_df['cum_length'][n_10]
                #     n_10 += 1
                # mainchannel_slope_usgs = (e_85 - e_10) / (l_10 - l_85)
                mainchannel_slope_usgs = mainchannel_slope

                print('!Mainchannel Length of %s basin CALCULATED [%.2fkm]!' % (sub_bsn_name, mainchannel_length))
                print('!Mainchannel Slope of %s basin CALCULATED [%.2f m/m]!' % (sub_bsn_name, mainchannel_slope))
                print('!Mainchannel Max Elevation of %s basin CALCULATED [%.2f msnm]!' % (sub_bsn_name,
                                                                                          mainchannel_max_elevation))
                print('!Mainchannel Outlet Elevation of %s basin CALCULATED [%.2f msnm]!' %
                      (sub_bsn_name, mainchannel_outlet_elevation))

                # main channel profile
                mainchannel_profile = stream_network_df.plot(x='cum_length', y='outlet_elev', color="blue",
                                                             label='Talweg')
                mainchannel_profile.set_xlabel('Abscisa [m]')
                mainchannel_profile.set_ylabel('Elevación [msnm]')
                mainchannel_profile_plot = mainchannel_profile.get_figure()
                mainchannel_profile_plot.savefig(out_dir + '/' + bsn_name + '/' + sub_bsn_name + '/' + 'profile.png')

                # cumulative areas over basin
                if coords_system == 'proj':
                    accumulation = stream_network_df['flow_accum'] * pixel_res ** 2 / 1000000   # square kilometers
                elif coords_system == 'geo':    # [area equals pixel resolution (degrees) by earth radius power of 2]
                    accumulation = stream_network_df['flow_accum'] * (pixel_res * math.pi / 180 * 6371) ** 2

                stream_network_df['flow_accum_km2'] = accumulation
                stream_network_df['relative_flow_accum'] = stream_network_df['flow_accum_km2'] / accumulation.max()
                stream_network_df['relative_elev'] = (stream_network_df['source_elev'] - mainchannel_outlet_elevation) \
                    / (mainchannel_max_elevation - mainchannel_outlet_elevation)

                # hypsometric curve
                mainchannel_hypsometric = stream_network_df.plot(x='relative_flow_accum', y='outlet_elev',
                                                                 color='brown', label='Curva hispométrica')
                mainchannel_hypsometric.set_xlabel('Area relativa [-]')
                mainchannel_hypsometric.set_ylabel('Elevación [msnm]')
                mainchannel_hypsometric_plot = mainchannel_hypsometric.get_figure()
                mainchannel_hypsometric_plot.savefig(out_dir + '/' + bsn_name + '/' + sub_bsn_name + '/' +
                                                     'hypsometric.png')
                plt.close(mainchannel_hypsometric_plot)

                # hypsographic curve
                mainchannel_hypsographic = stream_network_df.plot(x='relative_flow_accum', y='relative_elev',
                                                                  color='brown', label='Curva hipsográfica')
                mainchannel_hypsographic.set_xlabel('Area relativa [-]')
                mainchannel_hypsographic.set_ylabel('Elevación relativa')
                mainchannel_hypsographic_plot = mainchannel_hypsographic.get_figure()
                mainchannel_hypsographic_plot.savefig(out_dir + '/' + bsn_name + '/' + sub_bsn_name + '/' +
                                                      'hypsographic.png')
                plt.close(mainchannel_hypsographic_plot)

                # export elevation, hypsometric and hypsographic profile
                stream_network_df.to_csv(out_dir + '/' + bsn_name + '/' + sub_bsn_name + '/' + 'mainchannel_info.csv')
                print('!Mainchannel Profile of %s basin EXPORTED!' % sub_bsn_name)
                print('!Hypsometric and Hypsographic Curves of %s basin EXPORTED!' % sub_bsn_name)

                # distance to centroid [AFDP]
                grass.run_command('v.buffer', input=v_basin, output='buffered_basin', distance=0.001, overwrite=True,
                                  quiet=True)

                grass.run_command('v.extract', input='buffered_basin', type='centroid', output='centroid',
                                  overwrite=True, quiet=True)

                grass.run_command('v.type', input='centroid', output=v_centroid, from_type='centroid', to_type='point',
                                  overwrite=True, quiet=True)

                grass.run_command('v.db.addtable', map=v_centroid, table=v_centroid, overwrite=True, quiet=True)

                grass.run_command('v.db.addcolumn', map=v_centroid, columns='dist_mc int', overwrite=True, quiet=True)

                grass.run_command('v.distance', from_=v_centroid, from_type='point', to=v_mainchannel, to_type='line',
                                  column='dist_mc', upload='dist', output='distance', overwrite=True, quiet=True)

                grass.run_command('v.to.points', input='distance', output='distance_p', use='node', overwrite=True,
                                  quiet=True)

                grass.run_command('v.to.db', map='distance_p', layer='2', option='cat', columns='cat')

                grass.run_command('v.what.rast', map='distance_p', raster=r_basin_distance, layer='2', column='along',
                                  overwrite=True, quiet=True)

                # distance to outlet from nearest point to center of gravity
                dist_out2 = grass.read_command('v.db.select', flags='c', map='distance_p', layer='2', col='along',
                                               overwrite=True, quiet=True)
                dist_out2 = (dist_out2.split("\n"))
                dist_out_centroid = float(dist_out2[1]) / 1000
                print('!Distance to Centroid of %s basin CALCULATED [%.2fkm]!' % (sub_bsn_name, dist_out_centroid))

                # topological diameter
                # [https: // www.iihr.uiowa.edu / rmantilla / files / 2013 / 01 / GeomTopol_CALI.pdf]
                topological_diameter = len(stream_network_df)
                print('!Topology Diameter of %s basin CALCULATED [%.2f]!' % (sub_bsn_name, topological_diameter))

                # basin characteristic altitudes
                basin_elevation_info = grass.raster_info(r_basin_elevation)
                basin_min_elevation = basin_elevation_info['min']
                basin_max_elevation = basin_elevation_info['max']
                basin_elevation_delta = basin_max_elevation - basin_min_elevation
                print('!Characteristic Altitudes of %s basin CALCULATED!' % sub_bsn_name)

                # basin max length [AFDP]
                distance_info = grass.raster_info(r_basin_distance)
                basin_max_length = distance_info['max'] / 1000      # km

                # shape factor [AFDP]
                shape_factor = area_basin / basin_max_length ** 2

                # elongation Ratio [AFDP]
                elongation_ratio = (2 * math.sqrt(area_basin / math.pi)) / basin_max_length

                # mean hillslope length [AFDP]
                hillslope_info = grass.read_command('r.univar', map=r_basin_hillslope_distance)
                mean_hillslope_length = float(hillslope_info.split('\n')[9].split('=')[0].split(':')[1])

                # stream network statistics
                txt_horton = out_dir + '/' + bsn_name + '/' + sub_bsn_name + '/' + 'horton_stats.txt'

                grass.run_command('g.region', raster=r_basin_strahler)

                stream_stats = grass.read_command('r.stream.stats', stream_rast=r_basin_strahler,
                                                  direction=r_basin_drainage, elevation=r_basin_elevation)
                text_horton_stats = open(txt_horton, "w")
                text_horton_stats.write("%s" % stream_stats)
                text_horton_stats.close()
                stream_stats_summary = stream_stats.split('\n')[4].split('|')
                stream_stats_mom = stream_stats.split('\n')[8].split('|')
                max_order = stream_stats_summary[0]
                num_streams = stream_stats_summary[1]
                len_streams = stream_stats_summary[2]
                stream_freq = stream_stats_summary[5]
                bif_ratio = stream_stats_mom[0]
                len_ratio = stream_stats_mom[1]
                area_ratio = stream_stats_mom[2]
                slope_ratio = stream_stats_mom[3]

                # drainage density
                drainage_density = float(len_streams) / float(area_basin)

                # magnitude [Corrected]
                line_id = 22 + 2 * int(max_order)     # line identifier
                stream_stats_mag = stream_stats.split('\n')[line_id].split('|')
                magnitude_o = stream_stats_mag[1]

                # first order stream frequency [Corrected]
                first_order_stream_frequency = int(magnitude_o) / area_basin

                # Other index [AFDP]
                mass_coefficient = float(mean_elev) / float(area_basin)
                oro_coefficient = float((basin_elevation_delta / 1000) ** 2) / float(area_basin)
                stab_coefficient = float(1) / float(drainage_density)
                basin_width = float(area_basin) / float(basin_max_length)

                # ======================================================================================================
                # -- Populate morphometric parameters data frame
                # ======================================================================================================
                parameters_df[sub_bsn_name] = [basin_centroid_east,
                                               basin_centroid_north,
                                               float(dict_region_basin['n']),
                                               float(dict_region_basin['w']),
                                               float(dict_region_basin['s']),
                                               float(dict_region_basin['e']),
                                               area_basin,
                                               perimeter_basin,
                                               basin_max_length,
                                               basin_width,
                                               dist_out_centroid,
                                               basin_max_elevation,
                                               basin_min_elevation,
                                               basin_elevation_delta,
                                               mean_elev,
                                               mean_slope,
                                               length_orienting_vector,
                                               prevalent_orientation,
                                               mainchannel_length,
                                               mainchannel_slope,
                                               mainchannel_slope_usgs,
                                               mainchannel_max_elevation,
                                               mean_hillslope_length,
                                               magnitude_o,
                                               max_order,
                                               num_streams,
                                               len_streams,
                                               first_order_stream_frequency,
                                               drainage_density,
                                               stream_freq,
                                               bif_ratio,
                                               len_ratio,
                                               area_ratio,
                                               slope_ratio,
                                               compactness_coefficient,
                                               circularity_ratio,
                                               topological_diameter,
                                               elongation_ratio,
                                               shape_factor,
                                               mass_coefficient,
                                               oro_coefficient,
                                               stab_coefficient,
                                               asymmetry_coefficient]
                print('!MORPHOMETRIC CHARACTERIZATION of %s basin DONE!' % sub_bsn_name)

                # ======================================================================================================
                # -- Delete maps if flagged
                # ======================================================================================================
                if deleted_maps:
                    deleted_maps = [r_mask, r_basin_slope, r_basin_stream, r_basin_strahler, r_basin_drainage,
                                    r_basin_distance, r_basin_elevation, r_basin_half_basin, r_basin_accumulation,
                                    r_basin_slope_average, r_basin_height_average, r_basin_hillslope_distance]
                    grass.run_command('g.remove', type='raster', name=deleted_maps, flags='f')

                    deleted_vectors = [v_half_basin, v_mainchannel, v_centroid, v_mainchannel_split,
                                       v_basin_stream_network]
                    grass.run_command('g.remove', type='vector', name=deleted_vectors, flags='f')

            except TypeError:
                print('!Morphometric Characterization of %s basin FAILED!' % sub_bsn_name)

        # ==============================================================================================================
        # -- Restore original region
        # ==============================================================================================================
        grass.read_command('g.region', flags='p', region='original')

    # ==================================================================================================================
    # -- Export morphometric parameters data frame
    # ==================================================================================================================
    parameters_df.to_csv(out_dir + '/' + bsn_name + '/' + sub_bsn_name + '/' + 'Morphometric Characterization.csv')
    print('BASIN %s PROCESSING TIME = %.2f minutes' % (sub_bsn_name, (time() - basin_processing_time)/60))
    grass.run_command('g.remove', flags='f', type='region', name='original')
    print('BASIN MORPHOMETRIC CHARACTERIZATION ENDED')


def f_bsn_simple_morphometry(input_dem, analysis_points, analysis_points_field, slope_r, flow_direction_r,
                             flow_accumulation_r, stream_network_r, stream_network_v, downstream_distance_r,
                             max_hydrology_file):
    """
    Preprocess if r.basin is going to be applied to the hole study area [recommended].

    Parameters
    ----------
    input_dem: str
        DEM used to create flow accumulation, flow direction, and streams network maps
    basin_directory: str
        Folder path where morphometric characterization reports will be stored
    analysis_points: str
        Name of the vector containing the analysis points
    analysis_points_field: str
        Analysis point field containing datavalues for basins codification
    raster_names_df: object
        Dataframe [pandas] containing predefined raster map names
    vector_names_df: object
        Dataframe [pandas] containing predefined raster map names
    slope_r: str
        Slope raster map [%]
    flow_accumulation: str
        Flow accumulation map
    flow_direction: str
        Flow directions map
    stream_network_r: str
        Stream network map [thinned]
    stream_network_v: str
        Stream network vector [generated from thinned stream network map]
    stream_strahler_r: str
        Stream network raster classified with Strahler order
    downstream_distance_r: str
        Downstream distance map

    Returns
    -------
    TODO: what does this function returns

    """
    # ==================================================================================================================
    # -- Main imports and preparations
    # ==================================================================================================================
    print('START BASIN MORPHOMETRIC CHARACTERIZATION')
    # check dependencies
    check_r_bsn_requisites()

    # define location type
    if grass.locn_is_latlong():
        coords_system = 'geo'

    else:
        coords_system = 'proj'

    # global variables
    r_slope = slope_r
    r_elevation_global = input_dem
    r_drainage_global = flow_direction_r
    r_stream_network = stream_network_r
    r_accumulation_global = flow_accumulation_r
    r_downstream_distance_global = downstream_distance_r

    v_outlet_snap = 'outlet_snap'
    v_stream_network = stream_network_v

    # save current region
    grass.read_command('g.region', flags='p', save='original', overwrite=True)

    # ==================================================================================================================
    # -- Input points preparation
    # ==================================================================================================================
    # snap outlet to stream network [hardcoded to four times raster resolution]
    pixel_res = grass.raster_info(r_elevation_global)['ewres']
    try:
        grass.read_command('g.findfile', element='vector', file=v_outlet_snap)
        grass.run_command('v.db.connect', map=v_outlet_snap, table=analysis_points, flags='d')
        grass.run_command('r.stream.snap', input=analysis_points, output=v_outlet_snap,
                          stream_rast=stream_network_r, radius=pixel_res * 4, overwrite=True, quiet=True)
        grass.run_command('v.db.connect', map=v_outlet_snap, table=analysis_points, flags='o')
    except CalledModuleError:
        grass.run_command('r.stream.snap', input=analysis_points, output=v_outlet_snap,
                          stream_rast=stream_network_r, radius=pixel_res * 4, overwrite=True, quiet=True)
        grass.run_command('v.db.connect', map=v_outlet_snap, table=analysis_points, flags='o')

    # set conditions for analysis type
    sql_path = grass.read_command('db.databases', driver='sqlite').replace('\n', '')
    con = sqlite3.connect(sql_path)
    sql_stat = 'SELECT * FROM %s ' % analysis_points
    analysis_points_df = pd.read_sql_query(sql_stat, con)
    prefix_list = analysis_points_df[analysis_points_field]

    # ==================================================================================================================
    # -- Create database to store morphometric parameters estimations
    # ==================================================================================================================
    parameter_index = ['1. Easting Centroid of basin',
                       '2. Northing Centroid of basin',
                       '3. Rectangle containing basin N',
                       '4. Rectangle containing basin W',
                       '5. Rectangle containing basin S',
                       '6. Rectangle containing basin E',
                       '7. Basin Area [km^2]',
                       '8. Basin Perimeter [km]',
                       '9. Basin Max Length [km]',
                       '10. Basin Width [km]',
                       '11. Basin Distance to Center of Gravity [km]',
                       '12. Basin Max Elevation [masl]',
                       '13. Basin Min Elevation [masl]',
                       '14. Basin Elevation Difference [m]',
                       '15. Basin Mean Elevation [masl]',
                       '16. Basin Mean Slope [m/m]',
                       '17. Basin Length of Directing Vector [km]',
                       '18. Basin Prevalent Orientation [degree from north, counterclockwise]',
                       '19. Mainchannel Length [km]',
                       '20. Mainchannel Mean Slope [m/m]',
                       '21. Mainchannel Max Elevation [masl]',
                       '22. Index Compactness Coefficient',
                       '23. Index Circularity Ratio',
                       '24. Index Topological Diameter',
                       '25. Index Elongation Ratio',
                       '26. Index Shape Factor']
    parameters_df = pd.DataFrame(data=None, index=parameter_index)

    # ==================================================================================================================
    # -- Morphometric characterization of each point [loop]
    # ==================================================================================================================
    for i in range(0, len(analysis_points_df)):
        basin_processing_time = time()

        # set conditions on analysis type
        prefix = 'stream_' + str(int(prefix_list[i]))

        # variable names for each basin analysis
        r_mask = 'r_mask'
        r_bsn = prefix + '_basin'
        r_basin_slope = prefix + '_slope'
        r_basin_stream = prefix + '_stream'
        r_basin_drainage = prefix + '_drainage'
        r_basin_distance = prefix + '_dist2out'
        r_basin_elevation = prefix + '_elevation'
        r_basin_half_basin = prefix + '_half_basin'
        r_basin_accumulation = prefix + '_accumulation'
        r_basin_slope_average = prefix + '_slope_average'
        r_basin_height_average = prefix + '_height_average'
        r_basin_hillslope_distance = prefix + '_hillslope_distance'

        v_basin = prefix + '_basin'
        v_half_basin = prefix + '_half_basin'
        v_mainchannel = prefix + '_mainchannel'
        v_centroid = prefix + '_centroid'
        v_mainchannel_split = prefix + '_mainchannel_split'
        v_basin_stream_network = prefix + '_stream_network'

        # ==============================================================================================================
        # -- Basin Delineation
        # ==============================================================================================================
        try:
            # extract analysis point
            point_id = i + 1
            grass.run_command('v.extract', input=analysis_points, output='Pour_Point%s' % point_id,
                              type='point', where='cat=%s' % point_id, overwrite=True, quiet=True)
            grass.run_command('v.to.rast', input='Pour_Point%s' % point_id, output='Pour_Point%s' % point_id,
                              use='cat', type='point', overwrite=True, quiet=True)
            pour_point_coordinates = grass.read_command('v.info', map='Pour_Point%s' % point_id, flags='g')
            pour_point_coordinates = dict(x.split('=', 1) for x in pour_point_coordinates.split('\n') if '=' in x)

            # delineate basin
            grass.run_command('r.water.outlet', input=r_drainage_global, output=r_bsn,
                              coordinates=[pour_point_coordinates['east'], pour_point_coordinates['north']],
                              overwrite=True, quiet=True)

            print('!Delineation of %s basin DONE!' % prefix)
        except TypeError:
            print('!Delineation of %s basin FAILED!' % prefix)

        # ==============================================================================================================
        # -- Basin Vectorization
        # ==============================================================================================================
        try:
            grass.run_command('r.to.vect', input=r_bsn, output=v_basin, type='area', flags='sv', overwrite=True,
                              quiet=True)

            # add two columns to the table: area and perimeter
            grass.run_command('v.db.addcolumn', map=v_basin, columns='area double precision', quiet=True)
            grass.run_command('v.db.addcolumn', map=v_basin, columns='perimeter double precision', quiet=True)

            # populate perimeter column
            grass.run_command('v.to.db', map=v_basin, option='perimeter', units='kilometers', columns='perimeter',
                              quiet=True)
            # read perimeter
            tmp = grass.read_command('v.to.db', map=v_basin, option='perimeter', units='kilometers',
                                     columns='perimeter', flags='p')
            perimeter_basin = float(tmp.split('\n')[1].split('|')[1])

            # populate area column
            grass.run_command('v.to.db', map=v_basin, option='area', columns='area', units='kilometers', quiet=True)

            # read area
            tmp = grass.read_command('v.to.db', map=v_basin, option='area', units='kilometers', columns='area',
                                     flags='p')
            area_basin = float(tmp.split('\n')[1].split('|')[1])

            print('!Vectorization of %s basin DONE!' % prefix)
        except TypeError:
            print('!Vectorization of %s basin FAILED!' % prefix)

        # ==============================================================================================================
        # -- Mask and Cropping
        # ==============================================================================================================
        try:
            grass.run_command('r.mask', raster=r_bsn)

            # add mask
            expr = '%s = %s / %s' % (r_mask, r_bsn, r_bsn)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            # crop accumulation map
            expr = '%s = %s' % (r_basin_accumulation, r_accumulation_global)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            # crop elevation map
            expr = '%s = %s' % (r_basin_elevation, r_elevation_global)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            # crop flow directions map
            expr = '%s = %s' % (r_basin_drainage, r_drainage_global)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            # crop stream network map
            expr = '%s = %s' % (r_basin_stream, r_stream_network)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            # crop slope map
            expr = '%s = %s' % (r_basin_slope, slope_r)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            # crop stream network strahler vector
            grass.run_command('v.overlay', ainput=v_stream_network, binput=v_basin, operator='and', atype='line',
                              output=v_basin_stream_network, olayer='0,1,0', overwrite=True, quiet=True)

            # crop downstream distance [dist2out] and adjust
            expr = '%s = %s' % (r_basin_distance, r_downstream_distance_global)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)
            global_downstream_distance = grass.raster_info(r_basin_distance)['min']

            expr = '%s = %s - %0.4f' % (r_basin_distance, r_basin_distance, global_downstream_distance)
            grass.run_command('r.mapcalc', expression=expr, overwrite=True)

            print('!Making and Cropping of %s basin DONE!' % prefix)
        except TypeError:
            print('!Making and Cropping of %s basin FAILED!' % prefix)
        # remove mask
        grass.run_command('r.mask', flags='r')

        # ==========================================================================================================
        # -- Morphometric characterization
        # ==========================================================================================================
        try:
            # mean elevation
            grass.run_command("r.stats.zonal", base=r_bsn, cover=r_basin_elevation, method="average",
                              output=r_basin_height_average, quiet=True)
            mean_elev = grass.raster_info(r_basin_height_average)['min']
            grass.run_command('g.remove', type='raster', name=r_basin_height_average, flags='f')
            print('!Mean Elevation of %s basin CALCULATED [%s]!' % (prefix, mean_elev))

            # mean slope
            grass.run_command("r.stats.zonal", base=r_bsn, cover=r_basin_slope, method="average",
                              output=r_basin_slope_average)
            mean_slope = grass.raster_info(r_basin_slope_average)['min']
            mean_slope = mean_slope / 100   # m/m
            grass.run_command('g.remove', type='raster', name=r_basin_slope_average, flags='f')
            print('!Mean Slope of %s basin CALCULATED [%.2f m/m]!' % (prefix, mean_slope))

            # centroid and mean slope
            baricenter_slope = grass.read_command("r.volume", input=r_slope, clump=r_bsn)
            baricenter_slope = baricenter_slope.split()
            basin_centroid_east = float(baricenter_slope[33])
            basin_centroid_north = float(baricenter_slope[34])
            print('!Centroid of %s basin CALCULATED [%.2fE, %.2fN]!' % (prefix, basin_centroid_east,
                                                                        basin_centroid_north))
            # rectangle coordinates
            info_region_basin = grass.read_command("g.region", vect=v_basin, flags='m')
            dict_region_basin = dict(x.split('=', 1) for x in info_region_basin.split('\n') if '=' in x)
            print('!Rectangle Containing %s basin CALCULATED!' % prefix)

            # directing vector
            if coords_system == 'proj':
                delta_x = abs(float(basin_centroid_east) - float(pour_point_coordinates['east']))
                delta_y = abs(float(basin_centroid_north) - float(pour_point_coordinates['north']))
                length_orienting_vector = math.sqrt((delta_x ** 2) + (delta_y ** 2)) / 1000
                print('!Directing Vector of %s basin CALCULATED!' % prefix)
            else:
                length_orienting_vector = f_geodesic_dist_wgs84(float(basin_centroid_north),
                                                                float(pour_point_coordinates['north']),
                                                                float(basin_centroid_east),
                                                                float(pour_point_coordinates['east'])) / 1000
                print('!Directing Vector of %s basin CALCULATED!' % prefix)

            # prevalent orientation
            if coords_system == 'proj':
                if delta_y != 0:
                    prevalent_orientation = math.atan(abs(delta_x) / abs(delta_y))
                    prevalent_orientation = prevalent_orientation * 180 / math.pi
                    if delta_x >= 0 > delta_y:
                        prevalent_orientation = 180 - prevalent_orientation
                    elif delta_x < 0 and delta_y < 0:
                        prevalent_orientation = 180 + prevalent_orientation
                    elif delta_x < 0 <= delta_y:
                        prevalent_orientation = 360 - prevalent_orientation
                else:
                    if delta_x >= 0:
                        prevalent_orientation = 90
                    else:
                        prevalent_orientation = 270
                print('!Prevalent Orientation of %s basin CALCULATED [%.2f degrees]!' % (prefix,
                                                                                         prevalent_orientation))
            else:
                prevalent_orientation = f_geodesic_bearing(float(basin_centroid_north),
                                                           float(pour_point_coordinates['north']),
                                                           float(basin_centroid_east),
                                                           float(pour_point_coordinates['east']))
                print('!Prevalent Orientation of %s basin CALCULATED [%.2f]!' % (prefix, prevalent_orientation))

            # compactness coefficient [AFDP]
            compactness_coefficient = 0.28 * perimeter_basin / math.sqrt(area_basin)
            print('!Compactness Coefficient of %s basin CALCULATED [%.2f]!' % (prefix, compactness_coefficient))

            # circularity ratio
            circularity_ratio = (4 * math.pi * area_basin) / (perimeter_basin ** 2)
            print('!Circularity Ratio of %s basin CALCULATED [%.2f]!' % (prefix, circularity_ratio))

            # main channel length and slope
            sql_stat_hack = 'SELECT * FROM %s ' % v_basin_stream_network
            hack_df = pd.read_sql_query(sql_stat_hack, con)
            hack_order = hack_df['hack'].min()
            expr = 'hack=%s' % hack_order
            grass.run_command('v.extract', input=v_basin_stream_network, output=v_mainchannel_split, type='line',
                              where=expr, overwrite=True, quiet=True)
            grass.run_command("v.build.polylines", input=v_mainchannel_split, output=v_mainchannel, type='line',
                              cats='first', overwrite=True, quiet=True)
            grass.run_command('v.to.db', map=v_mainchannel, option='length', units='kilometers', columns='length',
                              quiet=True)

            sql_stat1 = 'SELECT * FROM %s ' % v_mainchannel_split
            stream_network_df = pd.read_sql_query(sql_stat1, con)
            stream_length = stream_network_df['length']
            stream_slope = stream_network_df['gradient']
            mainchannel_out_elevation = stream_network_df['outlet_elev']
            mainchannel_source_elevation = stream_network_df['source_elev']
            weighted_slope = stream_length * stream_slope / stream_length.sum()
            mainchannel_length = stream_length.sum() / 1000     # kilometers
            mainchannel_slope = weighted_slope.sum()            # m/m
            mainchannel_max_elevation = mainchannel_source_elevation.max()
            mainchannel_outlet_elevation = mainchannel_out_elevation.min()
            print('!Mainchannel Length of %s basin CALCULATED [%.2f]!' % (prefix, mainchannel_length))
            print('!Mainchannel Slope of %s basin CALCULATED [%.2f] percent!' % (prefix, mainchannel_slope * 100))
            print('!Mainchannel Max Elevation of %s basin CALCULATED [%.2fm]!' % (prefix,
                                                                                  mainchannel_max_elevation))
            print('!Mainchannel Outlet Elevation of %s basin CALCULATED [%.2fm]!' % (prefix,
                                                                                     mainchannel_outlet_elevation))

            # distance to centroid [AFDP]
            grass.run_command('v.buffer', input=v_basin, output='buffered_basin', distance=0.001, overwrite=True,
                              quiet=True)

            grass.run_command('v.extract', input='buffered_basin', type='centroid', output='centroid',
                              overwrite=True, quiet=True)

            grass.run_command('v.type', input='centroid', output=v_centroid, from_type='centroid', to_type='point',
                              overwrite=True, quiet=True)

            grass.run_command('v.db.addtable', map=v_centroid, table=v_centroid, overwrite=True, quiet=True)

            grass.run_command('v.db.addcolumn', map=v_centroid, columns='dist_mc int', overwrite=True, quiet=True)

            grass.run_command('v.distance', from_=v_centroid, from_type='point', to=v_mainchannel, to_type='line',
                              column='dist_mc', upload='dist', output='distance', overwrite=True, quiet=True)

            grass.run_command('v.to.points', input='distance', output='distance_p', use='node', overwrite=True,
                              quiet=True)

            grass.run_command('v.to.db', map='distance_p', layer='2', option='cat', columns='cat')

            grass.run_command('v.what.rast', map='distance_p', raster=r_basin_distance, layer='2', column='along',
                              overwrite=True, quiet=True)

            # distance to outlet from nearest point to center of gravity
            dist_out2 = grass.read_command('v.db.select', flags='c', map='distance_p', layer='2', col='along',
                                           overwrite=True, quiet=True)
            dist_out2 = (dist_out2.split("\n"))
            dist_out_centroid = float(dist_out2[1]) / 1000
            print('!Distance to Centroid of %s basin CALCULATED [%.2fkm]!' % (prefix, dist_out_centroid))

            # topological diameter
            # [https: // www.iihr.uiowa.edu / rmantilla / files / 2013 / 01 / GeomTopol_CALI.pdf]
            topological_diameter = len(stream_network_df)
            print('!Topology Diameter of %s basin CALCULATED [%.2f]!' % (prefix, topological_diameter))

            # basin characteristic altitudes
            basin_elevation_info = grass.raster_info(r_basin_elevation)
            basin_min_elevation = basin_elevation_info['min']
            basin_max_elevation = basin_elevation_info['max']
            basin_elevation_delta = basin_max_elevation - basin_min_elevation
            print('!Characteristic Altitudes of %s basin CALCULATED!' % prefix)

            # basin max length [AFDP]
            distance_info = grass.raster_info(r_basin_distance)
            basin_max_length = distance_info['max'] / 1000      # km

            # shape factor [AFDP]
            shape_factor = area_basin / basin_max_length ** 2

            # elongation Ratio [AFDP]
            elongation_ratio = (2 * math.sqrt(area_basin / math.pi)) / basin_max_length

            # Other index [AFDP]
            basin_width = float(area_basin) / float(basin_max_length)

            # ======================================================================================================
            # -- Populate morphometric parameters data frame
            # ======================================================================================================
            parameters_df[prefix] = [basin_centroid_east, basin_centroid_north, float(dict_region_basin['n']),
                                     float(dict_region_basin['w']), float(dict_region_basin['s']),
                                     float(dict_region_basin['e']), area_basin, perimeter_basin, basin_max_length,
                                     basin_width, dist_out_centroid, basin_max_elevation, basin_min_elevation,
                                     basin_elevation_delta, mean_elev, mean_slope, length_orienting_vector,
                                     prevalent_orientation, mainchannel_length, mainchannel_slope,
                                     mainchannel_max_elevation, compactness_coefficient, circularity_ratio,
                                     topological_diameter, elongation_ratio, shape_factor]

            print('!MORPHOMETRIC CHARACTERIZATION of %s basin DONE!' % prefix)

            # ==========================================================================================================
            # -- Delete maps if flagged
            # ==========================================================================================================
            deleted_maps = [r_mask, r_basin_slope, r_basin_stream, r_basin_drainage, r_basin_elevation,
                            r_basin_half_basin, r_basin_accumulation, r_basin_slope_average, r_basin_height_average,
                            r_basin_hillslope_distance, 'Pour_Point%s' % point_id]
            grass.run_command('g.remove', type='raster', name=deleted_maps, flags='f')

            deleted_vectors = [v_half_basin, v_mainchannel, v_mainchannel_split, v_basin_stream_network,
                               'Pour_Point%s' % point_id]
            grass.run_command('g.remove', type='vector', name=deleted_vectors, flags='f')
        except TypeError:
            print('!Morphometric Characterization of %s basin FAILED!' % prefix)
        print('%s OUT OF %s BASINS DONE' % (i + 1, len(analysis_points_df)))

        # ==============================================================================================================
        # -- Restore original region
        # ==============================================================================================================
        grass.read_command('g.region', flags='p', region='original')

    # ==================================================================================================================
    # -- Delete global vectors
    # ==================================================================================================================
    deleted_vectors = ['centroid', 'distance', 'distance_p', 'buffered_basin']
    grass.run_command('g.remove', type='vector', name=deleted_vectors, flags='f')

    # ==================================================================================================================
    # -- Export morphometric parameters data frame
    # ==================================================================================================================
    parameters_df.to_excel(max_hydrology_file, sheet_name='1. Morphometric Parameters')
    print('BASIN %s PROCESSING TIME = %.2f minutes' % (prefix, (time() - basin_processing_time)/60))
    grass.run_command('g.remove', flags='f', type='region', name='original')
    print('BASIN MORPHOMETRIC CHARACTERIZATION ENDED')


def f_pointer(flow_direction):
    """Tells the cell where to flow.

    Parameters
    ----------
    flow_direction : str
        Flow direction map

    Returns
    -------
    int
        Row and column pointing the cell where te actual cell is going to flow

    """
    # %% Flow direction pointer based on grass format
    if flow_direction == 1:
        i = -1
        j = 1
    elif flow_direction == 2:
        i = -1
        j = 0
    elif flow_direction == 3:
        i = -1
        j = -1
    elif flow_direction == 4:
        i = 0
        j = -1
    elif flow_direction == 5:
        i = 1
        j = -1
    elif flow_direction == 6:
        i = 1
        j = 0
    elif flow_direction == 7:
        i = 1
        j = 1
    elif flow_direction == 8:
        i = 0
        j = 1

    return i, j


def v_discharge_points(stream_network_r, flow_direction_r, discharges_r, discharges_v):
    """
    Extract discharges points [points upstream a convergence].

    Parameters
    ----------
    stream_network_r: str
        Stream network map name
    flow_direction_r: str
        Flow direction map name
    discharges_r: str
        Stream network discharges map
    discharges_v: str
        Stream network discharges vector [points]


    Returns
    -------
    Each discharge point in the defined stream network [inside GRASS mapset]

    """
    # read maps
    stream_rast = garray.array()
    stream_rast.read(stream_network_r)

    direction_rast = garray.array()
    direction_rast.read(flow_direction_r)

    n_rows = np.shape(stream_rast)[0]
    n_cols = np.shape(stream_rast)[1]

    # identify discharges points
    discharges_rast = garray.array()
    discharges_rast[discharges_rast == 0] = np.nan
    for row in range(0, n_rows):
        for col in range(0, n_cols):
            stream_code = stream_rast[row, col]
            if stream_code != 0:
                flow_dir = direction_rast[row, col]
                if flow_dir > 0:
                    i, j = f_pointer(flow_dir)
                    downstream_stream = stream_rast[row + i, col + j]
                    if downstream_stream != stream_code:
                        discharges_rast[row, col] = stream_code

    discharges_rast.write(discharges_r, overwrite=True)
    grass.run_command('r.to.vect', input=discharges_r, output=discharges_v, type='point', overwrite=True)

    # add coordinates to discharges points
    grass.run_command('v.db.addcolumn', map=discharges_v, column='east double')
    grass.run_command('v.db.addcolumn', map=discharges_v, column='north double')
    grass.run_command('v.to.db', map=discharges_v, option='coor', columns=['east', 'north'])

    print('!Discharges points identification DONE [AFDP]!')


def v_outlet_points(stream_network_r, flow_direction_r, outlets_r, outlets_v):
    """
    Extract outlet points [points that flow outside study region].

    Parameters
    ----------
    grass: object
        GRASS environment (grass created with grass_environment module)
    garray: object
        GRASS garray module
    stream_network_r: str
        Stream network map name
    flow_direction_r: str
        Flow direction map name
    outlets_r: str
        Stream network outlets map
    outlets_v: str
        Stream network outlets vector [points]

    Returns
    -------
    Each point flowing outside study region [inside GRASS mapset]

    """
    # TODO more efficient
    # read maps
    stream_rast = garray.array()
    stream_rast.read(stream_network_r)

    direction_rast = garray.array()
    direction_rast.read(flow_direction_r)

    n_rows = np.shape(stream_rast)[0]
    n_cols = np.shape(stream_rast)[1]

    # identify discharges points
    outlets_rast = garray.array()
    outlets_rast[outlets_rast == 0] = np.nan
    for row in range(0, n_rows):
        for col in range(0, n_cols):
            stream_code = stream_rast[row, col]
            if stream_code != 0:
                flow_dir = direction_rast[row, col]
                if flow_dir < 0:
                    outlets_rast[row, col] = stream_code

    outlets_rast.write(outlets_r, overwrite=True)
    grass.run_command('r.to.vect', input=outlets_r, output=outlets_v, type='point', overwrite=True)
    expr = '%s = int(%s)' % (outlets_r, outlets_r)
    grass.run_command('r.mapcalc', expression=expr, overwrite=True)
    print('!Discharges points identification DONE [AFDP]!')


def r_terrain_correction(dem_raster, hydro_corrected_dem, stream_network=None, correction_method=2,
                         stream_burn_method='1', epsilon=2.0):
    """
    Hydrological correction of DEM using SAGA GIS and WhiteBox GAT.

    Parameters
    ----------
    dem_raster: str
        Digital Elevation Model map name
    stream_network_raster: str
        Streams Network map name
    hydro_corrected_dem: str
        Hydrological Corrected map name
    correction_method: int
        Select hydrological correction method [1: Fill - Planchon & Darboux (SAGA), 2: Fill - Wang & Liu (SAGA),
        3: Fill - r.fill.dir (GRASS), 4: Breach depressions - Lindsay (WhiteBox),
        5: Stream Burn - Saunders (WitheBox, not recommended), 6: Stream Burn - O. Conrad (SAGA)]
    stream_burn_method: int
        Stream burn methods [1: simply decrease, 2: lower neighbors, 3: trace stream]
    epsilon: float
        Depression depth in stream burning functions

    Returns
    -------
    Hydrological corrected DEM [inside GRASS mapset]

    """
    # preparations
    wbt = whitebox.WhiteboxTools()
    # call(['saga_cmd', 'ta_preprocessor', '3'])
    # call(['saga_cmd', '-h'])
    temporal_stream_network = None
    metadata = None

    # create temporal folder
    temporal_folder = os.getcwd() + '/Saga_Temp'
    if not os.path.exists(temporal_folder):
        os.mkdir(temporal_folder)

    # export raster
    temporal_raster = temporal_folder + '/' + dem_raster + '.tif'
    grass.run_command('r.out.gdal', input=dem_raster, output=temporal_raster, format='GTiff', overwrite=True)
    if stream_network:
        temporal_stream_network = temporal_folder + '/' + stream_network + '.tif'
        grass.run_command('r.out.gdal', input=stream_network, output=temporal_stream_network, format='GTiff',
                          overwrite=True)

    # DTM hydrological correction algorithms
    if correction_method == 1:     # Fill with Planchon & Darboux (2001) algorithm
        call(['saga_cmd', 'ta_preprocessor', '3', '-DEM', temporal_raster, '-RESULT',
              temporal_folder + '/FILL_DEM', '-MINSLOPE', '0.01'])
        metadata = 'FILLED DEM with Planchon & Darboux (2001) algorithm, minimum slope = 0.01'
        temporal_raster = temporal_folder + '/FILL_DEM.sdat'
    elif correction_method == 2:   # Fill with Wang & Liu (2007) algorithm
        call(['saga_cmd', 'ta_preprocessor', '4', '-ELEV', temporal_raster, '-FILLED',
              temporal_folder + '/FILL_DEM', '-FDIR', temporal_folder + '/DIR', '-WSHED', temporal_folder + '/WSH',
              '-MINSLOPE', '0.1'])
        temporal_raster = temporal_folder + '/FILL_DEM.sdat'
        metadata = 'FILLED DEM with Wang & Liu (2007) algorithm, minimum slope = 0.1'
    elif correction_method == 3:   # Fill with GRASS r.fill.dir algorithm
        grass.run_command('r.fill.dir', input=dem_raster, output=hydro_corrected_dem, direction='DIR', overwrite=True)
        metadata = 'FILLED DEM with GRASS r.fill.dir algorithm'
        grass.run_command('g.remove', type='raster', name='DIR', flags='f')
    elif correction_method == 4:   # Breach with Lindsay(2016) algorithm
        wbt.set_working_dir(temporal_folder)
        wbt.verbose = False
        wbt.breach_depressions(temporal_raster, temporal_folder + '/BREACHED_DEM.tif')
        temporal_raster = temporal_folder + '/BREACHED_DEM.tif'
        metadata = 'BREACHED DEM with Lindsay (2016) algorithm, minimum slope = 0.1'
    elif correction_method == 5:   # Stream Burn - Saunders (WitheBox)
        wbt.set_working_dir(temporal_folder)
        wbt.verbose = False
        wbt.fill_burn(temporal_raster, temporal_stream_network, temporal_folder + '/BURNED_DEM.tif')
        temporal_raster = temporal_folder + '/BURNED_DEM.tif'
        metadata = 'STREAM BURNED DEM with Saunders (1999) algorithm'
    elif correction_method == 6:   # Stream Burn - O. Conrad (SAGA)
        call(['saga_cmd', 'ta_preprocessor', '6', '-DEM', temporal_raster, '-STREAM', temporal_stream_network,
              '-BURN', temporal_folder + '/BURNED_DEM', '-METHOD', str(stream_burn_method), '-EPSILON', str(epsilon)])
        temporal_raster = temporal_folder + '/BURNED_DEM.sdat'
        metadata = 'STREAM BURNED DEM with O. Conrad (2011) algorithm'

    # import results to GRASS working environment
    if correction_method != 3:     # if not a GRASS method
        grass.run_command('r.in.gdal', input=temporal_raster, output=hydro_corrected_dem, overwrite=True,
                          title=metadata)

    # delete temporal folder
    shutil.rmtree(temporal_folder, ignore_errors=True)


def r_watershed(dem_raster, threshold=500, flow_direction=None, flow_accumulation=None, stream_network=None,
                basins=None, half_basins=None, topographic_wetness_index=None, stream_power_index=None,
                slope_length=None, slope_steepness=None, flags='s'):
    """
    Hydrological processing using WATERSHEDS module.

    Parameters
    ----------
    dem_raster: str
        DEM name
    coordinates: list
        Set region coordinates where analysis is going to be performed [north, east, south, west]
    threshold: int
        Define flow accumulation threshold for stream initiation
    flow_direction: str
        Flow drainage direction map name
    flow_accumulation: str
        Flow accumulation map name
    stream_network: str
        Stream network map name
    basins: str
        Basins map name
    half_basins: str
        Half basins map name
    topographic_wetness_index: str
        Topographic wetness index map name
    stream_power_index: str
        Stream power index map name
    slope_length: str
        Slope length map name [used for U.S.L.E methods]
    slope_steepness: str
        Slope steepness map name [used for U.S.L.E methods]
    flags: str
        Set WATERSHEDS processing flags [default = s, single flow drainage directions]

    Returns
    -------
    Flow Drainage Direction map [FDIR]      >> Optional
    Flow Accumulation map [FACC]            >> Optional
    Stream Network map [STR]                >> Optional
    Basins map [BSN]                        >> Optional
    Half Basins map [HLF]                   >> Optional
    Topographic Wetness Index map [TCI]     >> Optional
    Stream Power Index map [POW]            >> Optional
    Slope Steepness map [STP]               >> Optional
    Slope Length map [SLN]                  >> Optional

    """
    grass.run_command('r.watershed', elevation=dem_raster, threshold=threshold, accumulation=flow_accumulation,
                      tci=topographic_wetness_index, spi=stream_power_index, drainage=flow_direction,
                      basin=basins, half_basin=half_basins, stream=stream_network,
                      length_slope=slope_length, slope_steepness=slope_steepness, flags=flags,
                      overwrite=True, quiet=True)
    print('!Terrain Hydrological Processing DONE [r.watershed]!')


def r_bsn_preprocess(dem_original_r, hydro_corrected_dem_r, flow_direction_r, flow_accumulation_r, stream_network_r,
                     out_raster_names_df, out_vector_names_df):
    """
    Preprocess if r.basin is going to be applied to the hole study area [recommended].

    Parameters
    ----------
    dem_original_r: str
        Hydrological corrected DEM used to create flow accumulation, flow direction, and streams network maps
    hydro_corrected_dem_r: str
        Hydrological corrected DEM used to create flow accumulation, flow direction, and streams network maps
    raster_names_df: object
        Dataframe [pandas] containing predefined raster map names
    vector_names_df: object
        Dataframe [pandas] containing predefined raster map names
    flow_direction_r: str
        Flow directions map
    flow_accumulation_r: str
        Flow accumulation map
    stream_network_r: str
        Stream network map
    Returns
    -------
    Slope map
    Aspect map
    Stream network thinned map
    Stream network vector
    Downstream length map

    """
    # check grass dependencies
    check_r_bsn_requisites()

    # ==================================================================================================================
    # -- Set raster maps and vector names
    # ==================================================================================================================
    # global variables
    r_drainage_global = flow_direction_r
    r_elevation_global = hydro_corrected_dem_r
    r_stream_network_global = stream_network_r
    r_accumulation_global = flow_accumulation_r
    r_original_elevation_global = dem_original_r

    # result maps
    r_hack_global = out_raster_names_df['Code']['hack_network_r']
    r_slope_global = out_raster_names_df['Code']['slope_original_r']
    r_shreve_global = out_raster_names_df['Code']['shreve_network_r']
    r_horton_global = out_raster_names_df['Code']['horton_network_r']
    r_aspect_global = out_raster_names_df['Code']['aspect_original_r']
    r_strahler_global = out_raster_names_df['Code']['strahler_network_r']
    r_hillslope_distance = out_raster_names_df['Code']['hillslope_dist_r']
    r_stream_network_outlets_global = out_raster_names_df['Code']['outlets_r']
    r_downstream_distance_global = out_raster_names_df['Code']['downstream_distance_r']
    r_stream_network_thin_global = out_raster_names_df['Code']['stream_network_thinned_r']

    # result vectors
    v_stream_network_outlets_global = out_vector_names_df['Code']['outlet_points']
    v_stream_network_global = out_vector_names_df['Code']['stream_network_polyline']

    # ==================================================================================================================
    # -- Terrain analysis
    # ==================================================================================================================
    overall_start_time = time()
    try:
        # ==============================================================================================================
        # -- Slope and aspect maps
        # ==============================================================================================================
        start_time = time()
        grass.run_command('r.slope.aspect', elevation=r_original_elevation_global, slope=r_slope_global,
                          aspect=r_aspect_global, format='percent', overwrite=True, quiet=True)
        print('!Slope and Aspect maps CALCULATED! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        # correct aspect map
        # start_time = time()
        # In Grass, aspect categories represent the number degrees of east and they increase
        # counterclockwise: 90deg is North, 180 is West, 270 is South 360 is East.
        # The aspect value 0 is used to indicate undefined aspect in flat areas with slope=0.
        # We calculate the number of degree from north, increasing counterclockwise.
        # expr = "%s = if(%s == 0, 0, if(%s > 90, 450 - %s, 90 - %s))" % (r_aspect_mod_global, r_aspect_global,
        #                                                                 r_aspect_global, r_aspect_global,
        #                                                                 r_aspect_global)
        # grass.run_command('r.mapcalc', expression=expr, overwrite=True)
        # print('!Aspect map CORRECTED! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        # ==============================================================================================================
        # -- Stream extract
        # ==============================================================================================================
        start_time = time()
        grass.run_command('r.thin', input=r_stream_network_global, output=r_stream_network_thin_global, overwrite=True,
                          quiet=True)
        print('!Stream raster thinned! - Process time = %.2f minutes' % ((time() - start_time) / 60))
        start_time = time()
        grass.run_command('r.to.vect', input=r_stream_network_thin_global, output=v_stream_network_global, type='line',
                          overwrite=True, quiet=True)
        print('!Streams vectors CALCULATED! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        # ==============================================================================================================
        # -- Stream order
        # ==============================================================================================================
        start_time = time()
        grass.run_command('r.stream.order', elevation=r_elevation_global, stream_rast=r_stream_network_global,
                          direction=r_drainage_global, accumulation=r_accumulation_global, hack=r_hack_global,
                          strahler=r_strahler_global, shreve=r_shreve_global, horton=r_horton_global,
                          stream_vect=v_stream_network_global, topo=None, overwrite=True, quiet=True)
        print('!Stream order maps CALCULATED! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        # ==============================================================================================================
        # -- Hill-slope distance to river network
        # ==============================================================================================================
        start_time = time()
        grass.run_command("r.stream.distance", stream_rast=r_stream_network_thin_global, direction=r_drainage_global,
                          elevation=r_elevation_global, distance=r_hillslope_distance, overwrite=True, quiet=True)
        print('!Distance to network maps calculated! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        # ==============================================================================================================
        # -- Identify outlet points
        # ==============================================================================================================
        start_time = time()
        # read maps
        stream_rast = garray.array()
        stream_rast.read(stream_network_r)

        direction_rast = garray.array()
        direction_rast.read(flow_direction_r)

        n_rows = np.shape(stream_rast)[0]
        n_cols = np.shape(stream_rast)[1]

        outlets_rast = garray.array()
        outlets_rast[outlets_rast == 0] = np.nan

        # identify outlet points
        for row in range(0, n_rows):
            for col in range(0, n_cols):
                stream_code = stream_rast[row, col]
                if stream_code != 0:
                    flow_dir = direction_rast[row, col]
                    if flow_dir < 0:
                        outlets_rast[row, col] = stream_code

        outlets_rast.write(r_stream_network_outlets_global, overwrite=True)
        grass.run_command('r.to.vect', input=r_stream_network_outlets_global, output=v_stream_network_outlets_global,
                          type='point', overwrite=True)
        expr = '%s = int(%s)' % (r_stream_network_outlets_global, r_stream_network_outlets_global)
        grass.run_command('r.mapcalc', expression=expr, overwrite=True)
        print('!Outlet points IDENTIFIED! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        # ==============================================================================================================
        # -- Downstream distance
        # ==============================================================================================================
        start_time = time()
        grass.run_command('r.stream.distance', stream_rast=r_stream_network_outlets_global, direction=flow_direction_r,
                          method='downstream', distance=r_downstream_distance_global, overwrite=True)
        print('!Downstream distance CALCULATED! - Process time = %.2f minutes' % ((time() - start_time) / 60))

        print('!MORPHOMETRIC PRE-PROCESSING TIME = %.2f minutes' % ((time() - overall_start_time) / 60))
    except TypeError:
        print('r.basin failed')


def r_hand_pointer(flow_direction_r, stream_network_r, hand_pointer_row_index, hand_pointer_column_index):
    """
    Construct HAND pointer.

    Parameters
    ----------
    flow_direction_r: str
        Flow directions map
    stream_network_r: str
        Stream network map [thinned]
    hand_pointer_row_index: str
        Output raster name
    hand_pointer_column_index: str
        Output raster name

    Returns
    -------
    HAND_pointer X Position
    HAND_pointer Y Position

    """
    # ==================================================================================================================
    # -- # Get region parameters
    # ==================================================================================================================
    region = grass.read_command('g.region', 'ug').splitlines()
    region_parmeters = []
    for i in enumerate(region):
        region_parmeters.append(region[i].split('='))
        region_parmeters[i][1] = float(region_parmeters[i][1])

    n_rows = region_parmeters[8][1]
    n_cols = region_parmeters[9][1]

    # ==================================================================================================================
    # -- # Construct HAND pointer rasters
    # ==================================================================================================================
    # raster to store the pointer for HAND algorithm
    hand_pointer_rows = garray.array()
    hand_pointer_columns = garray.array()

    hand_pointer_rows[hand_pointer_rows == 0] = np.nan
    hand_pointer_columns[hand_pointer_columns == 0] = np.nan

    # stream raster
    stream_rast = garray.array()
    stream_rast.read(stream_network_r)

    # directions raster
    flow_dir_rast = garray.array()
    flow_dir_rast.read(flow_direction_r)

    start_time = time()
    for row in range(0, int(n_rows)):
        for col in range(0, int(n_cols)):
            flow_direction_r = flow_dir_rast[row, col]
            is_stream = stream_rast[row, col]
            row_aux = row
            col_aux = col
            while flow_direction_r > 0 and is_stream == 0:
                i, j = f_pointer(flow_direction_r)
                row_aux += i
                col_aux += j
                flow_direction_r = flow_dir_rast[row_aux, col_aux]
                is_stream = stream_rast[row_aux, col_aux]

            if is_stream == 0:
                hand_pointer_rows[row, col] = np.nan
                hand_pointer_columns[row, col] = np.nan
            else:
                hand_pointer_rows[row, col] = row_aux
                hand_pointer_columns[row, col] = col_aux
        print(row, col)
    hand_pointer_rows.write(hand_pointer_row_index, overwrite=True)
    hand_pointer_columns.write(hand_pointer_column_index, overwrite=True)

    end_time = time()
    print(end_time - start_time, ' seconds')


def r_hand_stream_profiles(input_dem, stream_network_r, stream_network_elevation_r):
    """
    Get stream network elevations raster.

    Parameters
    ----------
    grass: object
        GRASS environment (created with grass_environment module)
    input_dem: str
        DEM used to create flow accumulation, flow direction, and streams network maps
    stream_network_r: str
        Stream network raster
    stream_network_elevation_r: str
        Stream network elevations raster [thalweg]

    Returns
    -------
    Stream network thalweg

    """
    expr = '%s = if(not(isnull(%s)), %s, null())' % (stream_network_elevation_r, stream_network_r, input_dem)
    grass.run_command('r.mapcalc', expression=expr, overwrite=True)


def r_hand_synthetic_inundation(input_dem, stream_network_r, stream_network_elevation_r, hand_pointer_row_index,
                                hand_pointer_column_index, water_depth, hand_inundation_wse, hand_inundation_depth,
                                hand_inundation_coding, stream_network_elevation_offset_r):
    """
    Inundation based on assumed water surface elevation.

    Parameters
    ----------
    input_dem: str
        DEM used to create flow accumulation, flow direction, and streams network maps
    stream_network_r: str
        Stream network raster
    stream_network_elevation_r: str
        Stream network elevations raster [thalweg]
    hand_pointer_row_index: str
        HAND row pointer index
    hand_pointer_column_index: str
        HAND column pointer index
    water_depth: float
        Assumed water depth [map set units]
    hand_inundation_wse: str
        Raster with inundation water surface elevation
    hand_inundation_depth: str
        Raster with inundation depth
    hand_inundation_coding: str
        Raster with inundation coded by stream
    stream_network_elevation_offset_r: str
        Raster with streams elevation plus water depth offset

    Returns
    -------
    Stream network thalweg

    """
    # ==================================================================================================================
    # -- Calculate water surface elevation
    # ==================================================================================================================
    expr = '%s = %s + %s' % (stream_network_elevation_offset_r, stream_network_elevation_r, water_depth)
    grass.run_command('r.mapcalc', expression=expr, overwrite=True)

    # ==================================================================================================================
    # -- Get region parameters
    # ==================================================================================================================
    region = grass.read_command('g.region', 'ug').splitlines()
    region_parameters = []
    for i in enumerate(region):
        region_parameters.append(region[i].split('='))
        region_parameters[i][1] = float(region_parameters[i][1])

    n_rows = region_parameters[8][1]
    n_cols = region_parameters[9][1]

    # ==================================================================================================================
    # -- Inundation with HAND
    # ==================================================================================================================
    # raster to store the pointer for HAND algorithm
    hand_pointer_rows = garray.array()
    hand_pointer_cols = garray.array()

    hand_pointer_rows.read(hand_pointer_row_index)
    hand_pointer_cols.read(hand_pointer_column_index)

    # read maps
    dem = garray.array()
    dem.read(input_dem)

    stream_network = garray.array()
    stream_network.read(stream_network_r)

    stream_network_energy = garray.array()
    stream_network_energy.read(stream_network_elevation_offset_r)

    # flood map
    flood = garray.array()
    flood[flood == 0] = np.nan

    wse = garray.array()
    wse[wse == 0] = np.nan

    # HAND inundation
    for row in range(0, int(n_rows)):
        for col in range(0, int(n_cols)):
            if not np.isnan(hand_pointer_rows[row, col]) and not np.isnan(hand_pointer_cols[row, col]):
                elevation = dem[row, col]
                energy = stream_network_energy[np.int(hand_pointer_rows[row, col]),
                                               np.int(hand_pointer_cols[row, col])]
                if energy >= elevation:
                    flood[row, col] = stream_network[np.int(hand_pointer_rows[row, col]),
                                                     np.int(hand_pointer_cols[row, col])]
                    wse[row, col] = energy

    flood.write(hand_inundation_coding, overwrite=True)
    wse.write(hand_inundation_wse, overwrite=True)

    expr = '%s = if(isnull(%s), null(), %s - %s)' % (hand_inundation_depth, hand_inundation_coding, hand_inundation_wse,
                                                     input_dem)
    grass.run_command('r.mapcalc', expression=expr, overwrite=True)

    expr = '%s = int(if(isnull(%s), null(), %s))' % (hand_inundation_coding, hand_inundation_coding,
                                                     hand_inundation_coding)
    grass.run_command('r.mapcalc', expression=expr, overwrite=True)


def r_hand_specific_inundation(input_dem, inundation_wse_r, hand_pointer_row_index, hand_pointer_column_index,
                               hand_inundation_map):
    """
    Inundation based on assumed water surface elevation.

    Parameters
    input_dem: str
        DEM used to create flow accumulation, flow direction, and streams network maps
    stream_network_r: str
        Stream network raster
    stream_network_elevation_r: str
        Stream network elevations raster [thalweg]
    hand_pointer_row_index: str
        HAND row pointer index
    hand_pointer_column_index: str
        HAND column pointer index
    streams_water_depth: object
        Assumed water depth [map set units]
    hand_inundation_wse: str
        Raster with inundation water surface elevation
    hand_inundation_depth: str
        Raster with inundation depth
    hand_inundation_coding: str
        Raster with inundation coded by stream

    Returns
    -------
    Stream network thalweg

    TODO: what does this function returns, check arguments, parameters etc

    """
    # ==================================================================================================================
    # -- Get region parameters
    # ==================================================================================================================
    region = grass.read_command('g.region', 'ug').splitlines()
    region_parameters = []
    for i in enumerate(region):
        region_parameters.append(region[i].split('='))
        region_parameters[i][1] = float(region_parameters[i][1])

    n_rows = region_parameters[8][1]
    n_cols = region_parameters[9][1]

    # ==================================================================================================================
    # -- Inundation with HAND
    # ==================================================================================================================
    hand_pointer_rows = garray.array()
    hand_pointer_cols = garray.array()

    hand_pointer_rows.read(hand_pointer_row_index)
    hand_pointer_cols.read(hand_pointer_column_index)

    # read maps
    dem = garray.array()
    dem.read(input_dem)

    wse = garray.array()
    wse.read(inundation_wse_r)

    # flood map
    flood = garray.array()
    flood[flood == 0] = np.nan

    # HAND inundation
    for row in range(0, int(n_rows)):
        for col in range(0, int(n_cols)):
            if not np.isnan(hand_pointer_rows[row, col]) and not np.isnan(hand_pointer_cols[row, col]):
                elevation = dem[row, col]
                energy = wse[np.int(hand_pointer_rows[row, col]), np.int(hand_pointer_cols[row, col])]
                if energy >= elevation:
                    flood[row, col] = 1

    flood.write(hand_inundation_map, overwrite=True)


def hand_hydraulic_parameters(slope_map, hand_inundation_coding, hand_inundation_depth, stream_network_v, zonal_map,
                              hand_depth, manning_map=None, manning=0.05):
    """
    Inundation based on assumed water surface elevation.

    Parameters
    ----------
    slope_map: str
        Slope map generated from original DEM
    stream_network_v: str
        Stream network vector
    hand_inundation_depth: str
        Raster with inundation depth
    hand_inundation_coding: str
        Raster with inundation coded by stream


    manning_map: str
        Raster containing manning values
    manning: float
        Generalized manning value [default = 0.05]

    Returns
    -------
    pandas.Dataframe
        HAND streams hydraulic parameters [Zheng et al, 2018]
    grass map
        Zonal statistics map with flooded stream codes and slope as label

    """
    # ==================================================================================================================
    # -- Calculate HAND stream parameters
    # ==================================================================================================================
    # load required maps
    slopes = garray.array()
    slopes.read(slope_map)

    hand_coding = garray.array()
    hand_coding.read(hand_inundation_coding)

    # calculate zonal statistics map with flooded stream codes and slope as label
    grass.run_command('r.stats.zonal', base=hand_inundation_coding, cover=slope_map, method='average', output=zonal_map,
                      flags='cr', overwrite=True)

    # get streams statistics
    stream_stats = grass.read_command('r.stats', input=zonal_map, flags='acl').splitlines()
    stream_stats = [i.split(' ') for i in stream_stats]
    stream_stats_codes = ['stream_' + i[0] for i in stream_stats[:-1]]
    stream_stats_slopes = [float(i[1]) for i in stream_stats[:-1]]      # percent slope
    stream_stats_areas = [float(i[2]) for i in stream_stats[:-1]]       # square meters area

    # get stream network attribute table in dataframe format
    sql_path = grass.read_command('db.databases', driver='sqlite').replace('\n', '')
    con = sqlite3.connect(sql_path)
    sql_stat = 'SELECT * FROM %s ' % stream_network_v
    stream_network_df = pd.read_sql_query(sql_stat, con)
    stream_codes = stream_network_df['a_cat']
    stream_lengths = stream_network_df['a_length']
    stream_slopes = stream_network_df['a_gradient']

    # get region resolution
    region = grass.read_command('g.region', 'ug').splitlines()
    region_parameters = []
    for i in enumerate(region):
        region_parameters.append(region[i].split('='))
        region_parameters[i][1] = float(region_parameters[i][1])

    # HAND hydraulic parameters for each stream
    complete_stream_codes = ['stream_' + str(i) for i in stream_codes]
    hydraulic_parameters = pd.DataFrame(np.zeros([len(stream_codes), 13]),
                                        index=complete_stream_codes,
                                        columns=['Stream_Length [m]',
                                                 'Stream_Slope [m/m]',
                                                 'Flooded_Area [m2]',
                                                 'Channel_Slope [m/m]',
                                                 'Channel_Bed_Area [m2]',
                                                 'Flooded_Volume [m3]',
                                                 'Average_Width [m]',
                                                 'Average_Cross_Section_Area [m2]',
                                                 'Average_Cross_Section_Wetted_Perimeter [m]',
                                                 'Average_Cross_Section_Hydraulic_Radius [m]',
                                                 'Stream_Discharge_Capacity [m3/s]',
                                                 'Stream_Average_Velocity [m/s]',
                                                 'Manning'])

    # fill known morphometric and hydraulic parameters [from stream network]
    hydraulic_parameters['Stream_Length [m]'] = stream_lengths.values
    hydraulic_parameters['Stream_Slope [m/m]'] = stream_slopes.values

    # delete duplicate rows in case stream network have this error [preferable to correct in vector map]
    hydraulic_parameters = hydraulic_parameters.loc[~hydraulic_parameters.index.duplicated(keep='first')]

    # fill known morphometric and hydraulic parameters [from zonal analysis]
    for i in enumerate(stream_stats_codes):
        hydraulic_parameters.loc[stream_stats_codes[i], 'Flooded_Area [m2]'] = stream_stats_areas[i]
        hydraulic_parameters.loc[stream_stats_codes[i], 'Channel_Slope [m/m]'] = stream_stats_slopes[i] / 100

    # calculate hydraulic parameters
    counter = 1
    for i in range(0, len(stream_stats_codes)):
        stream_code = hydraulic_parameters.index.values[i]
        stream_slope = hydraulic_parameters.loc[stream_code, 'Stream_Slope [m/m]']
        stream_length = hydraulic_parameters.loc[stream_code, 'Stream_Length [m]']
        flooded_area = hydraulic_parameters.loc[stream_code, 'Flooded_Area [m2]']
        channel_bedslope = hydraulic_parameters.loc[stream_code, 'Channel_Slope [m/m]']

        # correct zero slope values and short streams
        if stream_length <= 0.1:
            stream_length = 0.1
        if stream_slope <= 0.01:
            stream_slope = 0.01

        # calculate channel bedslope area
        channel_bedslope_factor = np.sqrt(1 + channel_bedslope ** 2)
        channel_bed_area = flooded_area * channel_bedslope_factor
        hydraulic_parameters.loc[stream_code, 'Channel_Bed_Area [m2]'] = channel_bed_area

        # calculate flood volume
        flood_volume = flooded_area * hand_depth
        hydraulic_parameters.loc[stream_code, 'Flooded_Volume [m3]'] = flood_volume

        # average flood width
        hydraulic_parameters.loc[stream_code, 'Average_Width [m]'] = flooded_area / stream_length

        # average cross section area
        flood_cross_section_average_area = flood_volume / stream_length
        hydraulic_parameters.loc[stream_code, 'Average_Cross_Section_Area [m2]'] = flood_cross_section_average_area

        # average wetted perimeter
        flood_cross_section_wetted_perimeter = channel_bed_area / stream_length
        hydraulic_parameters.loc[stream_code, 'Average_Cross_Section_Wetted_Perimeter [m]'] = \
            flood_cross_section_wetted_perimeter

        # average hydraulic radius
        flood_cross_section_hydraulic_radius = flood_cross_section_average_area / flood_cross_section_wetted_perimeter
        hydraulic_parameters.loc[stream_code, 'Average_Cross_Section_Hydraulic_Radius [m]'] =\
            flood_cross_section_hydraulic_radius

        # discharge capacity
        discharge_capacity = ((1 / manning) * flood_cross_section_average_area *
                              flood_cross_section_hydraulic_radius ** (2 / 3) * stream_slope ** 0.5)    # m3/s
        hydraulic_parameters.loc[stream_code, 'Stream_Discharge_Capacity [m3/s]'] = discharge_capacity

        # stream velocity
        stream_average_velocity = discharge_capacity / flood_cross_section_average_area
        hydraulic_parameters.loc[stream_code, 'Stream_Average_Velocity [m/s]'] = stream_average_velocity

        # Manning
        hydraulic_parameters.loc[stream_code, 'Manning'] = manning

        print(str(counter) + ' of ' + str(len(stream_stats_codes)) + ' streams')
        counter += 1

    return hydraulic_parameters
