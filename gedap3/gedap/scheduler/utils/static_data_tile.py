"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
'''
module static_data_tile
-----------------------

:author: Alix Damman
:date: Created on 06 Sep 2016

This module contains the class StaticDataTile.
'''

import netCDF4
import os

class StaticDataTile:
    '''
    This class loads the content of a STATIC netCDF file (see C++ documentation of struct StaticData (Common)).
    
    :param file: path of the STATIC netCDF file
    :type file: string
    
    .. todo:: 
    
       [static_data_tile.py] If you update the C++ structure StaticData in Data.h and modify the method save_static_data() from \
       the C++ TileNetCDFManager class, you must update this class too.
    '''
    
    def __init__(self, file):
        # if not os.path.isfile(file):
        #     logger.error('The file {} does not exist'.format(file))
        #     raise IOError('The file {} does not exist'.format(file))

        nc = netCDF4.Dataset(file)
        
        # ##################################################################################### #
        #                                       DIMENSIONS                                      # 
        # ##################################################################################### #
        
        self.nb_lines = len(nc.dimensions['line'])
        self.nb_columns = len(nc.dimensions['column'])
        
        # ##################################################################################### #
        #                                       VARIABLES                                       # 
        # ##################################################################################### # 
        group = nc.groups['Static_Data']
        self.longitude  = group.variables['Longitude'][:,:]
        self.latitude   = group.variables['Latitude'][:,:]
        self.land_cover = group.variables['LandCover'][:,:]
        self.elevation  = group.variables['Elevation'][:,:]
        
        # ##################################################################################### #
        #                                         END                                           # 
        # ##################################################################################### # 
        
        nc.close()
