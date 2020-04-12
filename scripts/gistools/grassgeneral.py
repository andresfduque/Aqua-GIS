#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""# Created by AndresD at 8/06/19.

Module that integrates best tools for geomorphometric analysis available in different GIS open source software

Features:
    + Create GRASS working environment

Pre-requisites:
    + GRASS GIS

@author:    Andres Felipe Duque Perez
Email:      aduquep@ingenevo.com.co

Credits:

GRASS Development Team, 2017. Geographic Resources Analysis Support System (GRASS) Software, Version 7.2.
Open Source Geospatial Foundation. Electronic document:. http://grass.osgeo.org
"""

# imports
import os
import sys
import subprocess
import numpy as np
from osgeo import gdal
import grass.script as grass


def grass_create_location(grass_bin, database_path, location, epsg):
    """Create GRASS location if it does not exist.

    Parameters
    ----------
    grass_bin: str
        CMD used to execute GRASS in your OS
    database_path: str
        Folder path defined to contain GRASS database
    location: str
        Location name
    epsg: int
        Coordinates system EPSG code

    """
    # create location
    if not os.path.isdir(database_path + '/' + location):
        start_cmd = grass_bin + ' -c epsg:' + str(epsg) + ' -e "' + database_path + '/' + location + '"'
        print(start_cmd)
        process = subprocess.Popen(start_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, err = process.communicate()
        if process.returncode != 0:
            print(sys.stderr, 'ERROR: %s' % err)
            print(sys.stderr, 'ERROR: Cannot generate location (%s)' % start_cmd)
            # sys.exit(-1)
        else:
            print('Created location %s' % database_path + '/' + location)
    else:
        print('Location already exists')


def append_grass_path(grass_path):
    """Append GRASS path to OS environmental variables.

    Parameters
    ----------
    grass_path : str
        GRASS library location

    """
    # define GRASS libraries path
    gis_base = os.environ['GISBASE'] = grass_path
    # Add GIS BASE to the python path
    sys.path.append(gis_base + '/etc/python')

    # elif sys.platform == 'win32' or sys.platform == 'win64':
    #     # Define GRASS libraries path
    #     gis_base = os.environ['GISBASE'] = grass_path
    #     # Add GISBASE to the python path
    #     sys.path.append(os.path.join(r'C:\OSGeo4W64\apps\grass\grass-7.0.5\etc\python'))
    #     import grass.script as grass


def grass_to_arcgis_fdir(flow_direction_grass, flow_direction_arcgis, output_folder):
    """Convert flow direction map from GRASS into ArcGIS.

        !---!---!---!		!----!----!-----!
        ! 3 ! 2 ! 1 !		! 32 ! 64 ! 128 !
        !---!---!---!		!----!----!-----!
        ! 4 !   ! 8 !		! 16 !    !  1  !
        !---!---!---!		!----!----!-----!
        ! 5 ! 6 ! 7 !		!  8 !  4 !  2  !
        !---!---!---!		!----!----!-----!
          !GRASS			!ArcGIS

    Parameters
    ----------
    flow_direction_grass: str
        Flow directions map in GRASS format
    flow_direction_arcgis: str
        Flow directions map in ARCGIS format
    output_folder: str
        Folder to store output maps

    """
    # constants
    fill_short = -32767

    # ==================================================================================================================
    # -- Export GRASS data to GeoTiff
    # ==================================================================================================================
    # files path
    grass_fdir_path = output_folder + '/' + flow_direction_grass + '.tif'
    arcgis_fdir_path = output_folder + '/' + flow_direction_arcgis + '.tif'

    # Export raster to GeoTiff format
    grass.run_command('r.out.gdal', input=flow_direction_grass, output=grass_fdir_path, nodata=fill_short,
                      type='Int16', format='GTiff', overwrite=True, quiet=True, flags="c")

    # ==================================================================================================================
    # -- Convert flow direction map
    # ==================================================================================================================
    # open the image
    data_source = gdal.Open(grass_fdir_path, gdal.GA_ReadOnly)
    if data_source is None:
        sys.exit('Could not open image dataset')

    # get image size
    n_rows = data_source.RasterYSize
    n_cols = data_source.RasterXSize

    # get geo-reference info
    trans = data_source.GetGeoTransform()
    projection = data_source.GetProjection()

    # get dir raster
    grass_fdir = data_source.ReadAsArray()

    # create empty
    arcgis_fdir = np.empty((n_rows, n_cols), dtype='uint16')
    arcgis_fdir.fill(fill_short)

    # convert flow direction map from ArcGIS into GRASS format
    arcgis_fdir[grass_fdir == 1] = 128
    arcgis_fdir[grass_fdir == 2] = 64
    arcgis_fdir[grass_fdir == 3] = 32
    arcgis_fdir[grass_fdir == 4] = 16
    arcgis_fdir[grass_fdir == 5] = 8
    arcgis_fdir[grass_fdir == 6] = 4
    arcgis_fdir[grass_fdir == 7] = 2
    arcgis_fdir[grass_fdir == 8] = 1
    arcgis_fdir[grass_fdir == 0] = 0
    arcgis_fdir[grass_fdir == -1] = 128
    arcgis_fdir[grass_fdir == -2] = 64
    arcgis_fdir[grass_fdir == -3] = 32
    arcgis_fdir[grass_fdir == -4] = 16
    arcgis_fdir[grass_fdir == -5] = 8
    arcgis_fdir[grass_fdir == -6] = 4
    arcgis_fdir[grass_fdir == -7] = 2
    arcgis_fdir[grass_fdir == -8] = 1

    # free memory
    data_source = None

    # create GeoTiff to store flow direction map in GRASS format
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(arcgis_fdir_path, n_cols, n_rows, 1, gdal.GDT_Byte,)

    # write the band
    dst_ds.GetRasterBand(1).WriteArray(arcgis_fdir)

    # set a no data value if required
    dst_ds.GetRasterBand(1).SetNoDataValue(fill_short)

    # geo-reference the image
    dst_ds.SetGeoTransform(trans)
    dst_ds.SetProjection(projection)

    # Close GeoTiff
    dst_ds.FlushCache()
