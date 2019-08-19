"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
"""
module input_tile
-----------------

:author: Alix Damman
:date: Created on 06 Sep 2016

This module contains the class InputTile.
"""

import numpy as np
import netCDF4


class InputTile(object):
    """
    This class loads the content of an INPUT netCDF file (see C++ documentation of structures Observations, ModelParams, Mask, InputData (Common)).

    :param file: path of the INPUT netCDF file
    :type file: string

    .. todo::

       [input_tile.py] If you update on of the C++ structures Observations, ModelParams, Mask, InputData in Data.h and modify the method save_input_data() from \
       the C++ TileNetCDFManager class, you must update this class too.
    """
    
    def __init__(self, file):
        # print file
        nc = netCDF4.Dataset(file)
        
        # ##################################################################################### #
        #                                       DIMENSIONS                                      # 
        # ##################################################################################### #
        
        self.nb_lines     = len(nc.dimensions['line'])
        self.nb_columns   = len(nc.dimensions['column'])
        self.nb_bands     = len(nc.dimensions['band'])
        
        # ##################################################################################### #
        #                                         BANDS                                         # 
        # ##################################################################################### #
        
        self.v_bands = []   
        for band_group in nc.groups['Observations'].groups:
            self.v_bands.append( band_group )
            
        self.nb_bands = len(self.v_bands)

        # ##################################################################################### #
        #                                       ATTRIBUTS                                       # 
        # ##################################################################################### #         
        
        self.b_under_illuminated = nc.__dict__['is_Under_illuminated']
         
        # ##################################################################################### #
        #                                       VARIABLES                                       # 
        # ##################################################################################### # 
        
        group_obs = nc.groups['Observations']
        self.Sun_Zenith_Angles   = group_obs.variables['Sun_Zenith_Angles'][:,:] 
        self.Sun_Azimuth_Angles  = group_obs.variables['Sun_Azimuth_Angles'][:,:] 
        self.View_Zenith_Angles  = group_obs.variables['View_Zenith_Angles'][:,:] 
        self.View_Azimuth_Angles = group_obs.variables['View_Azimuth_Angles'][:,:] 
                
        self.TOA_BRF           = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.Radiometric_Noise = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.Acquisition_Time  = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        
        i_band = 0
        for band_group in list(group_obs.groups.values()):
            self.TOA_BRF[i_band,:,:]             = band_group.variables['TOA_BRF'][:,:] 
            self.Radiometric_Noise[i_band,:,:]   = band_group.variables['Radiometric_Noise'][:,:] 
            self.Acquisition_Time[i_band,:,:]    = band_group.variables['Acquisition_Time'][:,:] 
            i_band += 1
        
        group_mod_param = nc.groups['ModelParams']
        self.TCO3             = group_mod_param.variables['TCO3'][:,:] 
        self.TCWV             = group_mod_param.variables['TCWV'][:,:] 
        self.Surface_Pressure = group_mod_param.variables['Surface_Pressure'][:,:] 
        self.Wind_Speed       = group_mod_param.variables['Wind_Speed'][:,:] 
        self.Wind_Direction   = group_mod_param.variables['Wind_Direction'][:,:] 
        
        group_mask = nc.groups['Mask']
        self.Mask = group_mask.variables['Mask'][:,:] 
        
        # ##################################################################################### #
        #                                         END                                           # 
        # ##################################################################################### # 
        
        nc.close()
