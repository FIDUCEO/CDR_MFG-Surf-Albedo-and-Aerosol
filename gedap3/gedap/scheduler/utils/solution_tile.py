"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
"""
module solution_tile
--------------------

:author: Alix Damman
:date: Created on 06 Sep 2016

This module contains the class SolutionTile.
"""

from .common_tile import Aerosol_Class
import numpy as np
import netCDF4


np.set_printoptions(linewidth=256, formatter={'float': '{: 0.5f}'.format})


class SolutionTile(object):
    """
    This class loads the content of a SOLUTION netCDF file (see C++ documentation of class TileNetCDFManager).

    Because of the of datasets required in product files for the AEROSOL_CCI project, this class implement 2 dedicated methods:
     * include_total_tau_55_with_uncertainty
     * include_total_fine_coarse_AOT_with_uncertainty

    see function :func:`scheduler.tasks.generate_products`.

    :param file: path of the PRIOR netCDF file
    :type file: string

    .. todo::

       [solution_tile.py] If you modify the method save_prior_data() from the C++ TileNetCDFManager class, you must update this class too.
    """
    
    def __init__(self, file):
        
        # print file
        nc = netCDF4.Dataset(file)

        self.source = nc.source
        self.auxiliary_data = nc.auxiliary_data
        
        # self.b_temporal_smoothness = nc.temporal_smoothness_term == "TRUE"
        # self.b_spectral_constraint = nc.spectral_constraint_term == "TRUE"
        
        # ##################################################################################### #
        #                                       DIMENSIONS                                      # 
        # ##################################################################################### #
        
        self.nb_lines     = len(nc.dimensions['line'])
        self.nb_columns   = len(nc.dimensions['column'])
        self.nb_obs       = len(nc.dimensions['max_obs'])
        self.nb_aer_class = 0
        self.nb_bands     = 0
        
        # ##################################################################################### #
        #                             BANDS + AER. CLASSES LIST                                 # 
        # ##################################################################################### #

        self.v_aerosol_class = []
        for aer_class_as_string in nc.list_aerosol_class.split(', '):
            st_name, st_type = aer_class_as_string.split(' ')
            self.v_aerosol_class.append(Aerosol_Class(st_name, st_type[1:-1]))
            # print 'aerosol class: ', self.v_aerosol_class[-1]
            
        self.nb_aer_class = len(self.v_aerosol_class)     

        self.v_bands = nc.list_bands.split(', ')
        # print 'v_bands ', self.v_bands
        self.nb_bands = len(self.v_bands)

        # print 'bands: ', self.v_bands
        # print 'aer. classes: ', self.v_aerosol_class
        # print ''

        # ##################################################################################### #
        #                                       VARIABLES                                       # 
        # ##################################################################################### #        
        
        self.nb_obs_per_pixel = nc.variables['nb_obs'][:,:]
        self.nb_iterations    = nc.variables['nb_iterations'][:,:]
        self.cost             = nc.variables['cost'][:,:]
        self.quality_flag     = nc.variables['quality_flag'][:,:,:]
        self.error_flag       = nc.variables['error_flag'][:,:]

        groupQI = nc.groups['Quality_Information']
        self.qi_aot_entropy = groupQI.variables['aot_entropy'][:,:,:]
        self.qi_convergence = groupQI.variables['convergence'][:,:,:]
        self.qi_cost = groupQI.variables['cost'][:,:,:]
        self.qi_jacobian = groupQI.variables['jacobian'][:,:,:]
        self.qi_rpv_entropy = groupQI.variables['rpv_entropy'][:,:,:]
        
        self.acquisition_time = nc.variables['acquisition_time'][:,:,:] 
        
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
        for band_group in list(nc.groups['RPV_Parameters'].groups.values()):
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
        #                                         AOT                                           # 
        # ##################################################################################### #       

        self.tau = np.full((self.nb_bands, self.nb_aer_class, self.nb_obs, self.nb_lines, self.nb_columns), np.nan)
        self.sigma_tau = np.full((self.nb_bands, self.nb_aer_class, self.nb_obs, self.nb_lines, self.nb_columns),
                                 np.nan)
        
        i_band = 0
        for band_group in list(nc.groups['AOT'].groups.values()):
            i_aer = 0
            for aer_class_group in list(band_group.groups.values()):
                self.tau[i_band, i_aer, :, :, :]       = aer_class_group.variables['tau'][:,:,:] 
                self.sigma_tau[i_band, i_aer, :, :, :] = aer_class_group.variables['sigma_tau'][:,:,:]  
                i_aer +=1
            i_band += 1
        
        # ##################################################################################### #
        #                                      TOTAL AOT                                        # 
        # ##################################################################################### #
        
        self.total_tau       = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_total_tau = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        
        i_band = 0
        for band_group in list(nc.groups['Total_AOT'].groups.values()):
            self.total_tau[i_band, :, :, :]       = band_group.variables['total_tau'][:,:,:] 
            self.sigma_total_tau[i_band, :, :, :] = band_group.variables['sigma_total_tau'][:,:,:]  
            i_band += 1
        
        # ##################################################################################### #
        #                                       AOT 550                                         # 
        # ##################################################################################### #    
        
        self.tau_55       = np.full( (self.nb_aer_class, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_tau_55 = np.full( (self.nb_aer_class, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        
        i_aer = 0
        for aer_class_group in list(nc.groups['AOT_550'].groups.values()):
            self.tau_55[i_aer, :, :, :]       = aer_class_group.variables['tau_55'][:,:,:] 
            self.sigma_tau_55[i_aer, :, :, :] = aer_class_group.variables['sigma_tau_55'][:,:,:]  
            i_aer += 1
        
        
        # ##################################################################################### #
        #                               Aerosol Class Properties                                # 
        # ##################################################################################### #    
        
        self.single_scattering_albedo       = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_single_scattering_albedo = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        self.asymmetry_factor               = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_asymmetry_factor         = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        
        i_band = 0
        for band_group in list(nc.groups['Aerosol_Class_Properties'].groups.values()):
            self.single_scattering_albedo[i_band, :, :, :]       = band_group.variables['single_scattering_albedo'][:,:,:] 
            self.sigma_single_scattering_albedo[i_band, :, :, :] = band_group.variables['sigma_single_scattering_albedo'][:,:,:]  
            self.asymmetry_factor[i_band, :, :, :]               = band_group.variables['asymmetry_factor'][:,:,:] 
            self.sigma_asymmetry_factor[i_band, :, :, :]         = band_group.variables['sigma_asymmetry_factor'][:,:,:] 
            i_band += 1
        
        # ##################################################################################### #
        #                                       TOA BRF                                         # 
        # ##################################################################################### #    
        
        self.TOA_BRF = np.full( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.nan )
        
        i_band = 0
        for band_group in list(nc.groups['TOA_BRF'].groups.values()):
            self.TOA_BRF[i_band, :, :, :] = band_group.variables['TOA_BRF'][:,:,:] 
            i_band += 1
        
        # ##################################################################################### #
        #                                         BHR                                           # 
        # ##################################################################################### #    
          
        self.BHR       = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        self.sigma_BHR = np.full( (self.nb_bands, self.nb_lines, self.nb_columns), np.nan )
        
        i_band = 0
        for band_group in list(nc.groups['BHR'].groups.values()):
            self.BHR[i_band, :, :]       = band_group.variables['BHR'][:,:] 
            self.sigma_BHR[i_band, :, :] = band_group.variables['sigma_BHR'][:,:]  
            i_band += 1
        
        # ##################################################################################### #
        #                                         END                                           # 
        # ##################################################################################### # 
        
        # print 'finished to read Solution Tile ', file
           
        nc.close()
    
    # ======================================================================================================== # 
    def include_total_tau_55_with_uncertainty(self):
        """
        Compute and add 2 Numpy arrays (for the project AEROSOL_CCI):
         * total_tau_55
         * sigma_total_tau_55

        see function :func:`scheduler.tasks.generate_products`.
        """
        
        self.total_tau_55       = np.sum(self.tau_55, axis=0)
        self.sigma_total_tau_55 = np.sqrt(np.sum(self.sigma_tau_55**2, axis=0))
        
    # ======================================================================================================== # 
    def include_total_fine_coarse_AOT_with_uncertainty(self):
        """
        Compute and add 8 Numpy arrays (for the project AEROSOL_CCI):
         * total_fine_tau + total_fine_tau_550
         * sigma_total_fine_tau + sigma_total_fine_tau_550
         * total_coarse_tau + total_coarse_tau_550
         * sigma_total_coarse_tau + sigma_total_coarse_tau_550

        see function :func:`scheduler.tasks.generate_products`.
        """
        
        self.total_fine_tau         = np.zeros( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        self.sigma_total_fine_tau   = np.zeros( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        self.total_coarse_tau       = np.zeros( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        self.sigma_total_coarse_tau = np.zeros( (self.nb_bands, self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        
        self.total_fine_tau_550         = np.zeros( (self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        self.sigma_total_fine_tau_550   = np.zeros( (self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        self.total_coarse_tau_550       = np.zeros( (self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        self.sigma_total_coarse_tau_550 = np.zeros( (self.nb_obs, self.nb_lines, self.nb_columns), np.float)
        
        for i_aer in range(self.nb_aer_class):
            
            if self.v_aerosol_class[i_aer].is_fine_mode():
                self.total_fine_tau_550[:,:,:] = self.total_fine_tau_550[:,:,:] + self.tau_55[i_aer, :, :, :]
                self.total_fine_tau[:,:,:,:]   = self.total_fine_tau[:,:,:,:] + self.tau[:, i_aer, :, :, :]
                
                self.sigma_total_fine_tau_550[:,:,:] = self.sigma_total_fine_tau_550[:,:,:] + np.square( self.sigma_tau_55[i_aer, :, :, :] )
                self.sigma_total_fine_tau[:,:,:,:]   = self.sigma_total_fine_tau[:,:,:,:] + np.square( self.sigma_tau[:, i_aer, :, :, :] )
            
            if self.v_aerosol_class[i_aer].is_coarse_mode():
                self.total_coarse_tau_550[:,:,:] = self.total_coarse_tau_550[:,:,:] + self.tau_55[i_aer, :, :, :]
                self.total_coarse_tau[:,:,:,:]   = self.total_coarse_tau[:,:,:,:] + self.tau[:, i_aer, :, :, :]
                
                self.sigma_total_coarse_tau_550[:,:,:] = self.sigma_total_coarse_tau_550[:,:,:] + np.square( self.sigma_tau_55[i_aer, :, :, :] )
                self.sigma_total_coarse_tau[:,:,:,:]   = self.sigma_total_coarse_tau[:,:,:,:] + np.square( self.sigma_tau[:, i_aer, :, :, :] )
    
        self.sigma_total_fine_tau_550   = np.sqrt( self.sigma_total_fine_tau_550 )
        self.sigma_total_coarse_tau_550 = np.sqrt( self.sigma_total_coarse_tau_550 )
        self.sigma_total_fine_tau       = np.sqrt( self.sigma_total_fine_tau )
        self.sigma_total_coarse_tau     = np.sqrt( self.sigma_total_coarse_tau )
