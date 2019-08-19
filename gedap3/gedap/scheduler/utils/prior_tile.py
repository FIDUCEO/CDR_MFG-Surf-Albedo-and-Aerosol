"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
"""
module prior_tile
-----------------

:author: Alix Damman
:date: Created on 06 Sep 2016

This module contains the class PriorTile.
"""

from .common_tile import Aerosol_Class
import numpy as np
import netCDF4

np.set_printoptions(linewidth=256, formatter={'float': '{: 0.5f}'.format})


class PriorTile(object):
    """
    This class loads the content of a PRIOR netCDF file (see C++ documentation of struct PriorData (Common)).

    :param file: path of the PRIOR netCDF file
    :type file: string

    .. todo::

       [prior_tile.py] If you update the C++ structure PriorData in Data.h and modify the method save_prior_data() from \
       the C++ TileNetCDFManager class, you must update this class too.
    """
    
    def __init__(self, file):

        nc = netCDF4.Dataset(file)
        
        # ##################################################################################### #
        #                                       DIMENSIONS                                      # 
        # ##################################################################################### #
        
        self.nb_lines     = len(nc.dimensions['line'])
        self.nb_columns   = len(nc.dimensions['column'])
        self.nb_aer_class = 0
        self.nb_bands     = 0
        
        # ##################################################################################### #
        #                             BANDS + AER. CLASSES LIST                                 # 
        # ##################################################################################### #

        self.v_aerosol_class = []
        for aer_class_as_string in nc.list_aerosol_class.split(', '):
            st_name, st_type = aer_class_as_string.split(' ')
            self.v_aerosol_class.append( Aerosol_Class(st_name, st_type[1:-1]) )
            
        self.nb_aer_class = len(self.v_aerosol_class)     

        self.v_bands = nc.list_bands.split(', ')
        self.nb_bands = len(self.v_bands)

        # print 'bands: ', self.v_bands
        # print 'aer. classes: ', self.v_aerosol_class
        # print ''

        # ##################################################################################### #
        #                                       VARIABLES                                       # 
        # ##################################################################################### #        
        
        self.last_update_time = nc.variables['last_update_time'][:,:]




        # ##################################################################################### #
        #                                      RPV PARAMS                                       # 
        # ##################################################################################### #
        
        self.rho_0       = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_rho_0 = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.k           = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_k     = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.theta       = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_theta = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.rho_c       = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_rho_c = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )

        i_band = 0
        for band_group in list(nc.groups.values()):
            self.rho_0[i_band, :, :]       = band_group.variables['rho_0'][:,:]
            self.sigma_rho_0[i_band, :, :] = band_group.variables['sigma_rho_0'][:,:]
            self.k[i_band, :, :]           = band_group.variables['k'][:,:]
            self.sigma_k[i_band, :, :]     = band_group.variables['sigma_k'][:,:]
            self.theta[i_band, :, :]       = band_group.variables['theta'][:,:]
            self.sigma_theta[i_band, :, :] = band_group.variables['sigma_theta'][:,:]
            self.rho_c[i_band, :, :]       = band_group.variables['rho_c'][:,:]
            self.sigma_rho_c[i_band, :, :] = band_group.variables['sigma_rho_c'][:,:]
            i_band += 1

        # ##################################################################################### #
        #                                      AOT PARAMS                                       #
        # ##################################################################################### #

        self.aot_coarse = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.aot_fine_abs = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.aot_fine_non_abs = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_aot = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )

        i_band = 0
        for band_group in list(nc.groups.values()):
            self.aot_coarse[i_band, :, :] = band_group.variables['aot_coarse'][:,:]
            self.aot_fine_abs[i_band, :, :] = band_group.variables['aot_fine_abs'][:,:]
            self.aot_fine_non_abs[i_band, :, :] = band_group.variables['aot_fine_non_abs'][:,:]
            self.sigma_aot[i_band, :, :] = band_group.variables['sigma_aot'][:,:]
            i_band += 1

        # ##################################################################################### #
        #                                         END                                           # 
        # ##################################################################################### # 
        
        # print 'finished to read Prior Tile ', file
        
        nc.close()
