"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
"""
module MFG_filenames_generator
------------------------------

:author: Cedric Goossens
:date: Created on May 2017

"""

from datetime import date, datetime, timedelta
from os.path import isfile
from collections import namedtuple
from glob import glob
import os

# ==================================================================================================================== #
#                                                  NAMED TUPLES                                                        #
# ==================================================================================================================== #

# This namedtuple stores the paths of files required to build STATIC netCDF files by the C++ Tile_Maker module
# (see C++ documentation).
InputFilesStaticTile = namedtuple('InputFilesStaticTile',
                                  ['path_lat_lon_file', 'path_land_sea_mask_file', 'path_elevation_file'])

# This namedtuple stores the paths of files required to build INPUT netCDF files by the C++ Tile_Maker module
# (see C++ documentation).
InputFilesTile = namedtuple('InputFilesTile',
                            ['date_and_time', 'slot', 'path_image_file', 'path_mask_file', #'path_angles_file',
                             'path_modparams_file', 'path_latlon_file', 'path_elevation_file','path_geopotential_file',
                             'path_aerosol_elevation_file_1', 'path_aerosol_elevation_file_2'])



# ==================================================================================================================== #
#                                         MFG FILENAMES GENERATOR CLASS                                                #
# ==================================================================================================================== #
class MFGFilenamesGenerator(object):
    """
    This class is used to generate all the paths of files containing relevant input data associated with the processing
    of MVIRI image. \
    These paths are required by the classes of the C++ Tile_Maker module to build STATIC and INPUT netCDF files
    (see C++ Documentation).

    :param dict_dir_external_data: dictionary (data type, directory path) where data types are static, image, angle,
        modparams and mask.
    :type dict_dir_external_data: dictionary[string,string]
    """

    def __init__(self, dict_dir_external_data, prefix, suffix):

        #: dictionary where values are directories containing image, angle, model parameters and mask files respectively
        self.directories = dict_dir_external_data

        self.prefix = prefix
        self.suffix = suffix

        #: index of the first considered slot.
        # the third slot corresponds to the half hour 01:30 -> 02:00
        self.first_slot_day = 0
        #: index of the last considered slot.
        # the 46th slot corresponds to the half hour 23:00 -> 23:30
        self.last_slot_day = 47
        self.isFirstTileOfTheDay = True

    def generate_list_static_data_files(self, mission, ssp, projection):
        """
        Generate the list of paths of files required to build STATIC netCDF files by the C++ StaticReader class \
        (see C++ documentation/code).

        :return: a namedtuple :class:`InputFilesStaticTile` containing the list of paths required to build STATIC netCDF
            files.
        :rtype: namedtuple :class:`InputFilesStaticTile`
        """

        directory = self.directories['static']
        ssp = str(ssp).zfill(4)

        # lat_lon_file = os.path.join(dict_dir_external_data["angle"], self.platform,
        #                             'FIDUCEO_FCDR_L15_MVIRI_MET7-00.0_STATIC_v2.2_fv2.5.nc')
        # The lat_lon_file to be used in the generation of static tiles should be the one provided
        # by the user, in this case the FIDUCEO static file (the one used for the angles)

        lat_lon_file = os.path.join(
            directory,
            "Latitude_Longitude_{mission}_{ssp}_{projection}.nc".format(
                mission=mission,
                ssp=ssp,
                projection=projection))
        land_sea_mask_file = os.path.join(
            directory,
            "LandSeaMask_{mission}_{ssp}_{projection}.nc".format(
                mission=mission,
                ssp=ssp,
                projection=projection))
        elevation_file = os.path.join(
            directory,
            "Elevation_{mission}_{ssp}_{projection}.nc".format(
                mission=mission,
                ssp=ssp,
                projection=projection))

        input_files_static_tile = InputFilesStaticTile(lat_lon_file, land_sea_mask_file, elevation_file)
        return input_files_static_tile

    def build_path_ERA_INTERIM_file(self, date, logger):
        """
        Generate the path to an ERA-INTERIM file corresponding to a given month.

        :param date: date of the ERA-INTERIM file.
        :type date: date
        :return: path to the ERA-INTERIM file.
        :rtype: string
        """
        st_year_dir = date.strftime('%Y')
        st_month_dir = date.strftime('%m')

        # build path for file associated with ERA-INTERIM data
        era_interim_file = glob(os.path.join(self.directories['modparams'], st_year_dir, st_month_dir, '*.nc'))
        if len(era_interim_file)== 0:
            logger.error(
                'No ERA-INTERIM file associated to year {} month {} can be found'.format(st_year_dir,st_month_dir))
            raise RuntimeError(
                'No ERA-INTERIM file associated to year {} month {} can be found'.format(st_year_dir,st_month_dir))
        else:
            era_interim_file = era_interim_file[0]
        return era_interim_file

    def get_input_files_for_one_tile(self, date, slot, mission, ssp, projection, platform, ecmwf_constant,
                                     aerosol_height_constant, daily_cloudmask, list_tile_input_files, logger):
        """
        Generate the list of paths of files required to build INPUT netCDF files associated with a given time slot \
        (see C++ documentation/code of the Tile_Maker module). \
        The list of paths is store in a namedtuple :class:`InputFilesTile` that is appended to a list list_tile_input_files.

        :param date: date corresponding to a MVIRI image
        :type date: date
        :param slot: slot of the MVIRI image (the slot is the index number of a MVIRI observation. There are 4 slots per hour and 96 per day)
        :type slot: int
        :param list_tile_input_files: list of input files that need to be extend to include the current slot.
        :type: list[string]
        """
        slot_hour = (slot * 30) // 60
        slot_hour2 = ((slot + 1) * 30) // 60
        slot_min = (slot * 30) % 60
        slot_min2 = ((slot + 1) * 30) % 60

        input_datetime = datetime(date.year, date.month, date.day, slot_hour, slot_min)

        # strings associated with date and slot
        st_year_dir = date.strftime('%Y') + "/"
        st_month_dir = date.strftime('%m') + "/"
        st_opmission_dir = platform

        path_image_file = os.path.join(
            self.directories['image'],
            st_opmission_dir,
            st_year_dir,
            st_month_dir,
            "{prefix}{date1}{time1}_{date2}{time2}{suffix}".format(
                prefix=self.prefix['image'],
                date1=date.strftime('%Y') + date.strftime('%m') + date.strftime('%d'),
                time1='{0:02d}'.format(slot_hour) + '{0:02d}'.format(slot_min),
                date2=date.strftime('%Y') + date.strftime('%m') + date.strftime('%d'),
                time2='{0:02d}'.format(slot_hour2) + '{0:02d}'.format(slot_min2),
                suffix=self.suffix['image']))

        if daily_cloudmask:
            path_mask_file = os.path.join(
                self.directories['mask'],
                st_year_dir,
                st_month_dir,
                "{prefix}{date}{suffix}".format(
                    prefix=self.prefix['mask'],
                    date=date.strftime('%Y') + date.strftime('%m') + date.strftime('%d'),
                    suffix=self.suffix['mask']))

        else:
            path_mask_file = os.path.join(
                self.directories['mask'],
                st_year_dir,
                st_month_dir,
                "{prefix}{date}{time}{suffix}".format(
                    prefix=self.prefix['mask'],
                    date=date.strftime('%Y') + date.strftime('%m') + date.strftime('%d'),
                    time='{0:02d}'.format(slot_hour) + '{0:02d}'.format(slot_min),
                    suffix=self.suffix['mask']))




        # checks if new ERA-INTERIM data needs to be read (they are given every 6 hours each day).
        # If true, we build the path of the file associated with the ERA-INTERIM data
        path_modparams_file = ""
        if not ecmwf_constant:
            path_modparams_file = self.build_path_ERA_INTERIM_file(date, logger)


        path_geopotential_file = os.path.join(
            self.directories['static'],
            "Geopotential_{mission}_{ssp}_{projection}.nc".format(
                mission=mission,
                ssp=ssp,
                projection=projection))

        suffix_array = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

        if not aerosol_height_constant:
            month = date.month
            day = date.day
            if day >= 10:
                path_aerosol_elevation_file_1 = os.path.join(
                    self.directories['static'],
                    "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                    "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                        suffix=suffix_array[int(date.strftime('%m')) - 1],
                        mission=mission,
                        ssp=ssp,
                        projection=projection))

                if month == 12:
                    path_aerosol_elevation_file_2 = os.path.join(
                        self.directories['static'],
                        "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                        "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                            suffix=suffix_array[0],
                            mission=mission,
                            ssp=ssp,
                            projection=projection))
                else:
                    path_aerosol_elevation_file_2 = os.path.join(
                        self.directories['static'],
                        "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                        "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                            suffix=suffix_array[month],
                            mission=mission,
                            ssp=ssp,
                            projection=projection))
            else:
                if (month - 2) == -1:
                    path_aerosol_elevation_file_1 = os.path.join(
                        self.directories['static'],
                        "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                        "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                            suffix=suffix_array[11],
                            mission=mission,
                            ssp=ssp,
                            projection=projection))
                    path_aerosol_elevation_file_2 = os.path.join(
                        self.directories['static'],
                        "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                        "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                            suffix=suffix_array[month - 1],
                            mission=mission,
                            ssp=ssp,
                            projection=projection))
                else:
                    path_aerosol_elevation_file_1 = os.path.join(
                        self.directories['static'],
                        "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                        "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                            suffix=suffix_array[month - 2],
                            mission=mission,
                            ssp=ssp,
                            projection=projection))
                    path_aerosol_elevation_file_2 = os.path.join(
                        self.directories['static'],
                        "AEROSOL_HEIGHT/AEROSOL_HEIGHT_"
                        "{suffix}_{mission}_{ssp}_{projection}.nc".format(
                            suffix=suffix_array[month - 1],
                            mission=mission,
                            ssp=ssp,
                            projection=projection))
        else:
            day = 0
            suffix1 = ""
            suffix2 = ""

        path_latlon_file = os.path.join(
            self.directories['static'],
            "Latitude_Longitude_{mission}_{ssp}_{projection}.nc".format(
                mission=mission,
                ssp=str(ssp).zfill(4),
                projection=projection))

        path_elevation_file = os.path.join(
            self.directories['static'],
            "Elevation_{mission}_{ssp}_{projection}.nc".format(
                mission=mission,
                ssp=ssp,
                projection=projection))

        # if files exist, create an InputFilesTile structure and return it

        if isfile(path_image_file) and isfile(path_mask_file): #and isfile(path_angles_file)
            list_tile_input_files.append(
                InputFilesTile(input_datetime, slot, path_image_file, path_mask_file, #path_angles_file,
                               path_modparams_file, path_latlon_file, path_elevation_file, path_geopotential_file,
                               path_aerosol_elevation_file_1, path_aerosol_elevation_file_2 ))

            if self.isFirstTileOfTheDay:
                self.isFirstTileOfTheDay = False
        elif not isfile(path_image_file):
            logger.info('The {} image is missing.'.format(path_image_file))
        # elif not isfile(path_angles_file,):
        #     logger.info('The {} angle file is missing.'.format(path_angles_file,))
        elif not isfile(path_mask_file):
            logger.info('The {} mask file is missing.'.format(path_mask_file))
        elif not isfile(path_latlon_file):
            logger.info('The {} latlon file is missing.'.format(path_latlon_file))


    def generate_list_input_files(self, start_date, end_date, mission, ssp, projection, platform, ecmwf_constant,
                                  aerosol_height_constant, daily_cloudmask, logger):
        """
        Generate the list of namedtuple :class:`InputFilesTile` associated with a given accumulation period. \
        see method :func:`get_input_files_for_one_tile`.

        :param start_date: starting date of the accumulation period.
        :type start_date: datetime
        :param end_date: ending date of the accumulation period.
        :type end_date: datetime
        :return: list of namedtuple :class:`InputFilesTile`
        :rtype: list[InputFilesTile]
        """

        # list containing namedtuple InputFilesTile --> see method get_input_files_for_one_tile
        list_tile_input_files = []

        current_date = start_date

        while current_date <= end_date:
            self.isFirstTileOfTheDay = True
            for slot in range(self.first_slot_day, self.last_slot_day + 1):
                self.get_input_files_for_one_tile(current_date, slot, mission, str(ssp).zfill(4), projection, platform,
                                                  ecmwf_constant, aerosol_height_constant, daily_cloudmask,
                                                  list_tile_input_files, logger)
            current_date += timedelta(days=1)


        return list_tile_input_files
