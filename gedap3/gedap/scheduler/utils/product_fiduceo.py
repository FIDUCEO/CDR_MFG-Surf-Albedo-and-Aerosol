"""
  The copyrights for the GEDAP algorithm and computer codes, remain with
  Rayference SPRL as an Intellectual Property Right

  Rayference Copyright (c)

"""
"""
module product_fiduceo
-----------------------------

:author: Elisa Pinat
:date: Created on 10 August 2018

This module contains the class ProductFiduceo.
"""

import os
import uuid
import numpy as np
import netCDF4
from datetime import datetime, timedelta


class ProductFiduceo(object):

    def __init__(self):

        self.line = 5000
        self.column = 5000

        self.date_hour = None

        self.latitude = np.full((self.line, self.column), -999.0, np.float)
        self.longitude = np.full((self.line, self.column), -999.0, np.float)
        self.time = np.full((self.line, self.column), -999.0, np.float)
        self.landcover = np.full((self.line, self.column), 0, np.int)

        self.AOT = np.full((self.line, self.column), -999.0, np.float)
        self.AOT_uncertainty = np.full((self.line, self.column), -999.0, np.float)

        self.BHR = np.full((self.line, self.column), -999.0, np.float)
        self.BHR_uncertainty = np.full((self.line, self.column), -999.0, np.float)
        self.quality_flag = np.full((self.line, self.column), -999.0, np.float)

        self.nb_processed_pixels = None
        self.nb_pixels_not_enough_obs = None
        self.mean_nb_iterations = None
        self.mean_cost = None
        self.total_nb_pixels = self.line * self.column
        self.nb_pixels_outside_Earth_disk = None

        self.BHR_Q1 = None
        self.BHR_Q2 = None
        self.BHR_Q3 = None
        self.BHR_mean = None

        self.AOT_Q1 = None
        self.AOT_Q2 = None
        self.AOT_Q3 = None
        self.AOT_mean = None

        self.cost_Q1 = None
        self.cost_Q2 = None
        self.cost_Q3 = None
        self.cost_mean = None

        self.nb_iterations_Q1 = None
        self.nb_iterations_Q2 = None
        self.nb_iterations_Q3 = None
        self.nb_iterations_mean = None

    # ======================================================================================================== #

    def write_netCDF(self, file):

        nc = netCDF4.Dataset(file, "w", format="NETCDF4")
        self.date_hour = datetime.strptime(os.path.basename(file).split('-')[0], '%Y%m%d')

        # ##################################################################################### #
        #                                       GROUPS                                          #
        # ##################################################################################### #

        group_static = nc.createGroup('STATIC')
        group_aod = nc.createGroup('AOT')
        group_bhr = nc.createGroup('BHR')
        group_quality = nc.createGroup('QUALITY')
        group_statistics = nc.createGroup('STATISTICS')

        # ##################################################################################### #
        #                                       DIMENSIONS                                      #
        # ##################################################################################### #

        nc.createDimension('y_vis', self.line)
        nc.createDimension('x_vis', self.column)
        group_static.createDimension('y_vis', self.line)
        group_static.createDimension('x_vis', self.column)
        group_aod.createDimension('y_vis', self.line)
        group_aod.createDimension('x_vis', self.column)
        group_bhr.createDimension('y_vis', self.line)
        group_bhr.createDimension('x_vis', self.column)
        group_quality.createDimension('y_vis', self.line)
        group_quality.createDimension('x_vis', self.column)
        # group_statistics.createDimension('d', None)
        # ##################################################################################### #
        #                                   GLOBAL ATTRIBUTS                                    #
        # ##################################################################################### #

        lat_min = np.nanmin(np.where(self.latitude > -200.0, self.latitude, np.nan))
        lat_max = np.nanmax(np.where(self.latitude > -200.0, self.latitude, np.nan))
        lon_min = np.nanmin(np.where(self.longitude > -200.0, self.longitude, np.nan))
        lon_max = np.nanmax(np.where(self.longitude > -200.0, self.longitude, np.nan))

        time_cov_start = self.date_hour.strftime("%Y%m%dT%H0000Z")
        time_cov_end = (self.date_hour + timedelta(hours=24)).strftime("%Y%m%dT%H0000SZ") #daily product

        nc.title = "Fiduceo Beta Product Example"
        nc.id = os.path.basename(file)
        nc.version = "Beta"
        nc.tracking_id = str(uuid.uuid4())
        nc.summary = "Level 2 aerosol properties retrieved using satellite data"
        nc.source = "MVIRI"
        nc.time_coverage_start = time_cov_start
        nc.time_coverage_end = time_cov_end
        nc.sensor = "MVIRI"
        nc.platform = "MFG"
        nc.resolution = "HIGH_RES"
        nc.projection = "geostationary"
        nc.geospatial_lat_min = lat_min
        nc.geospatial_lat_max = lat_max
        nc.geospatial_lon_min = lon_min
        nc.geospatial_lon_max = lon_max

        nc.conventions = "CF-1.6"
        nc.naming_authority = "rayference.eu"
        nc.date_created = datetime.now().strftime("%Y%m%dT%H%M%SZ")
        nc.project = "The Fidelity and Uncertainty in Climate Data Records from Earth Observation (FIDUCEO)-H2020"

        nc.creator_name = "Rayference"
        nc.creator_URL = "http://rayference.eu"
        nc.creator_email = "info@rayference.eu"
        nc.keywords = "satellite,observation,atmosphere"

        # ##################################################################################### #
        #                                       STATISTICS                                      #
        # ##################################################################################### #

        pixels = group_statistics.createGroup('PIXEL_INFO')

        nb_pixels_tot = pixels.createVariable('total_nb_pixels', np.int)
        nb_pixels_processed = pixels.createVariable('nb_processed_pixels', np.int)
        nb_pixels_not_enough_obs = pixels.createVariable('nb_pixels_not_enough_obs', np.int)
        nb_pixels_outside_Earth = pixels.createVariable('nb_pixels_outside_Earth_disk', np.int)

        nb_pixels_tot[:] = self.total_nb_pixels
        nb_pixels_processed[:] = self.nb_processed_pixels
        nb_pixels_not_enough_obs[:] = self.nb_pixels_not_enough_obs
        nb_pixels_outside_Earth[:] = self.nb_pixels_outside_Earth_disk

        cost = group_statistics.createGroup('COST_INFO')
        cost_Q1 = cost.createVariable('cost_Q1', np.float)
        cost_Q2 = cost.createVariable('cost_Q2', np.float)
        cost_Q3 = cost.createVariable('cost_Q3', np.float)
        cost_mean = cost.createVariable('cost_mean', np.float)

        cost_Q1[:] = self.cost_Q1
        cost_Q2[:] = self.cost_Q2
        cost_Q3[:] = self.cost_Q3
        cost_mean[:] = self.cost_mean

        nb_iterations = group_statistics.createGroup('NB_ITERATIONS_INFO')
        nb_iterations_Q1 = nb_iterations.createVariable('nb_iterations_Q1', np.float)
        nb_iterations_Q2 = nb_iterations.createVariable('nb_iterations_Q2', np.float)
        nb_iterations_Q3 = nb_iterations.createVariable('nb_iterations_Q3', np.float)
        nb_iterations_mean = nb_iterations.createVariable('nb_iterations_mean', np.float)

        nb_iterations_Q1[:] = self.nb_iterations_Q1
        nb_iterations_Q2[:] = self.nb_iterations_Q2
        nb_iterations_Q3[:] = self.nb_iterations_Q3
        nb_iterations_mean[:] = self.nb_iterations_mean

        info_aot = group_statistics.createGroup('AOT_INFO')
        AOT_Q1 = info_aot.createVariable('AOT_Q1', np.float)
        AOT_Q2 = info_aot.createVariable('AOT_Q2', np.float)
        AOT_Q3 = info_aot.createVariable('AOT_Q3', np.float)
        AOT_mean = info_aot.createVariable('AOT_mean', np.float)

        AOT_Q1[:] = self.AOT_Q1
        AOT_Q2[:] = self.AOT_Q2
        AOT_Q3[:] = self.AOT_Q3
        AOT_mean[:] = self.AOT_mean

        info_bhr = group_statistics.createGroup('BHR_INFO')
        BHR_Q1 = info_bhr.createVariable('BHR_Q1', np.float)
        BHR_Q2 = info_bhr.createVariable('BHR_Q2', np.float)
        BHR_Q3 = info_bhr.createVariable('BHR_Q3', np.float)
        BHR_mean = info_bhr.createVariable('BHR_mean', np.float)

        BHR_Q1[:] = self.BHR_Q1
        BHR_Q2[:] = self.BHR_Q2
        BHR_Q3[:] = self.BHR_Q3
        BHR_mean[:] = self.BHR_mean

        # nc.setncatts({'total_nb_pixels': self.total_nb_pixels,
        #               'nb_processed_pixels': self.nb_processed_pixels,
        #               'nb_pixels_not_enough_obs': self.nb_pixels_not_enough_obs,
        #               'nb_pixels_outside_Earth_disk': self.nb_pixels_outside_Earth_disk,
        #               'mean_nb_iterations': self.mean_nb_iterations,
        #               'mean_cost': self.mean_cost })


        # ##################################################################################### #
        #                                       VARIABLES                                       #
        # ##################################################################################### #

        nc_latitude = group_static.createVariable('Latitude', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_latitude.long_name = "Latitude"
        nc_latitude.standard_name = "latitude"
        nc_latitude.units = "degrees_north"
        nc_latitude.valid_range = [-90.0, 90.0]

        nc_latitude[:, :] = self.latitude[:, :]

        nc_longitude = group_static.createVariable('Longitude', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_longitude.long_name = "Longitude"
        nc_longitude.standard_name = "longitude"
        nc_longitude.units = "degrees_east"
        nc_longitude.valid_range = [-180.0, 180.0]
        nc_longitude[:, :] = self.longitude[:, :]

        nc_landcover = group_static.createVariable('LandCover', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_landcover.long_name = "Land Cover Type"
        nc_landcover.standard_name = "LandCover"
        nc_landcover.units = ""
        nc_landcover.valid_range = [0, 4]
        nc_landcover[:, :] = self.landcover[:, :]

        nc_time = nc.createVariable('time', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_time.long_name = "Modified Julian Date"
        nc_time.standard_name = "time"
        nc_time.units = "MJD"
        nc_time[:, :] = self.time[:, :]

        # ================================ AOT ================================ #

        nc_AOT = group_aod.createVariable('AOT', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_AOT.long_name = "Aerosol optical thickness in the MVIRI band"
        nc_AOT.standard_name = "atmosphere_optical_thickness_due_to_ambient_aerosol"
        nc_AOT.units = "1"
        nc_AOT.valid_range = [0.0, 5.0]

        nc_AOT[:, :] = self.AOT[:, :]

        # ========================= AOD uncertainty =========================== #

        nc_AOT_uncertainty = group_aod.createVariable('AOT_uncertainty', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_AOT_uncertainty.long_name = "AOT_uncertainty"
        nc_AOT_uncertainty.units = "1"
        nc_AOT_uncertainty.valid_range = [0.0, 1000.0]

        nc_AOT_uncertainty[:, :] = self.AOT_uncertainty[:, :]

        # ========================= BHR ======================================== #

        nc_BHR = group_bhr.createVariable('BHR', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_BHR.long_name = "BHR"
        nc_BHR.units = "1"
        nc_BHR.valid_range = [0.0, 1.0]

        nc_BHR[:, :] = self.BHR[:, :]

        # ========================= BHR uncertainty ============================ #

        nc_BHR = group_bhr.createVariable('BHR_uncertainty', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_BHR.long_name = "BHR_uncertainty"
        nc_BHR.units = "1"
        nc_BHR.valid_range = [0.0, 1.0]

        nc_BHR[:, :] = self.BHR_uncertainty[:, :]

        # ========================= Quality flag ============================ #

        nc_quality = group_quality.createVariable('quality_flag', np.float, ('y_vis', 'x_vis'), fill_value=-999.0)
        nc_quality.long_name = "Quality flag"
        nc_quality.units = "1"
        nc_quality.valid_range = [0.0, 1.0]

        nc_quality[:, :] = self.quality_flag[:, :]

        # ##################################################################################### #
        #                                         END                                           #
        # ##################################################################################### #
        # print 'Finished writing Products file ', file
        nc.close()

    # ======================================================================================================== #
    def read_netCDF(self, file):

        nc = netCDF4.Dataset(file)
        self.date_hour = datetime.strptime( os.path.basename(file).split('-')[0], '%Y%m%d' )

        # ##################################################################################### #
        #                                       DIMENSIONS                                      #
        # ##################################################################################### #

        self.line = nc.dimensions['y_vis'].size
        self.column = nc.dimensions['x_vis'].size

        # ##################################################################################### #
        #                                       VARIABLES                                       #
        # ##################################################################################### #

        self.time[:, :] = nc.variables['time'][:, :]

        # ##################################################################################### #
        #                                    GROUP STATIC                                       #
        # ##################################################################################### #

        self.latitude[:, :] = nc.groups['STATIC'].variables['Latitude'][:, :]
        self.longitude[:, :] = nc.groups['STATIC'].variables['Longitude'][:, :]
        self.landcover[:, :] = nc.groups['STATIC'].variables['LandCover'][:, :]

        # ##################################################################################### #
        #                                    GROUP AOD                                          #
        # ##################################################################################### #

        self.AOT[:, :] = np.ma.masked_array(nc.groups['AOT'].variables['AOT'][:, :], nc.groups['AOT'].variables['AOT'][:,:].mask).filled(-999.0)
        self.AOT_uncertainty[:, :] = nc.groups['AOT'].variables['AOT_uncertainty'][:, :]

        # ##################################################################################### #
        #                                    GROUP BHR                                          #
        # ##################################################################################### #

        self.BHR[:, :] = np.ma.masked_array(nc.groups['BHR'].variables['BHR'][:, :], nc.groups['BHR'].variables['BHR'][:, :].mask).filled(-999.0)
        self.BHR_uncertainty[:, :] = nc.groups['BHR'].variables['BHR_uncertainty'][:, :]
        self.BHR_uncertainty[self.landcover==2] = -999.0

        # ##################################################################################### #
        #                                    GROUP QUALITY                                      #
        # ##################################################################################### #

        self.quality_flag[:, :] = np.ma.masked_array(nc.groups['QUALITY'].variables['quality_flag'][:, :],
                                                     nc.groups['QUALITY'].variables['quality_flag'][:, :].mask).filled(-999.0)
        self.quality_flag[:, :] = nc.groups['QUALITY'].variables['quality_flag'][:, :]

        # ##################################################################################### #
        #                                         END                                           #
        # ##################################################################################### #

        # print 'Finished reading Products file ', file

        nc.close()

