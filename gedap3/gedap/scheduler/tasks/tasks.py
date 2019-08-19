"""
  The copyrights for the GEDAP algorithm and computer codes, remain with
  Rayference SPRL as an Intellectual Property Right

  Rayference Copyright (c)

"""

import os
import time
import subprocess
import shutil
import filecmp
import pprint

from glob import glob
from datetime import timedelta
from collections import OrderedDict

from ..setup.tile_handler import initialize_tilemaker, initialize_tileprocessor
from ..utils import tile as tile_utils
from ..utils.jdutil import jd_to_mjd, datetime_to_jd
from ..utils.static_data_tile import StaticDataTile
from ..report.logger import print_and_log


# ==================================================================================================================== #
#                                                 GEDAP FUNCTIONS                                                      #
# ==================================================================================================================== #


def update_accumulation_periods_list(ctx, period_index):
    """
    Update the list of accumulation periods and the JSON file list_accumulation_periods.json (see class :class:`scheduler.context.Context`). \
    Basically, this function must be called after finishing to generate all product files associated with a given accumulation period \
    (see function :func:`scheduler.main.run_scheduler`).

    :param period_index: index of the accumulation period for which all tiles have been processed and all product files generated.
    :type period_index: int
    """

    ctx.v_acc_periods[period_index].b_processed = True
    ctx.dump_list_accumulation_periods_in_json_file()


# ==================================================================================================================== #
def prepare_to_process_first_acc_period(ctx, logger):
    """
    First method to call when you start to run GEDAP. Uses an instance of the C++ TileMaker and TileProcessor classes.

    This function executes the following steps:
     * create STATIC netCDF files associated with all active tiles (call the C++ method make_static_data_tile of the class TileMaker)
     * create all INPUT netCDF files associated with the first accumulation to process (call the function :func:`make_input_tiles`)
     * create a PRIOR netCDF file with default values (call the C++ method create_default_prior_tile of the class TileProcessor).

    .. warning:: netCDF files associated with tiles are created only if they do not exist yet (no replacement).

    .. seealso:: :func:`scheduler.main.run_scheduler`
    """

    # ============================================= #
    #   GET FIRST UNPROCESSED ACCUMULATION PERIOD   #
    # ============================================= #

    first_acc_period = ctx.v_acc_periods[ctx.index_first_acc_period_to_process]

    logger.info("    --> create static data tiles \n" +
                "    --> create default Prior tiles \n" +
                "    --> read first ERA-INTERIM file \n" +
                "    --> create input tiles for the first accumulation period"
                       )

    tilemaker = initialize_tilemaker(ctx, logger)

    tileprocessor = initialize_tileprocessor(ctx, first_acc_period.start_date)

    # ==================================================== #
    #   CREATE DEFAULT PRIOR TILES IF DOES NOT EXIST YET   #
    # ==================================================== #

    start_date_mjd = jd_to_mjd(datetime_to_jd(first_acc_period.start_date))
    for tile in ctx.v_tiles:
        dir_prior_tiles = os.path.join(ctx.rootdir_tile_maker[tile.machine], ctx.subdir_prior_tiles)


        tile_prior_file = os.path.join(dir_prior_tiles,
                                       tile_utils.get_tile_filename("PRIOR", ctx.satellite, ctx.radiometer,
                                                                    ctx.projection,
                                                                    str(ctx.ssp).zfill(4), ctx.platform,
                                                                    first_acc_period.start_date.strftime('%Y-%m-%d'),
                                                                    first_acc_period.start_date.strftime('%H-%M'),
                                                                    tile.id))


        if not os.path.isfile(tile_prior_file) or ctx.cold_start:
            os.remove(tile_prior_file) if os.path.exists(tile_prior_file) else None
            dir_input_tiles = os.path.join(ctx.rootdir_tile_maker[tile.machine] , ctx.subdir_input_tiles)
            file_static_tile = os.path.join(dir_input_tiles,
                                            tile_utils.get_tile_filename("STATIC", ctx.satellite, ctx.radiometer,
                                                                         ctx.projection, str(ctx.ssp).zfill(4),
                                                                         ctx.platform,"","", tile.id))

            try:
                static_tile = StaticDataTile(file_static_tile)
            except IOError:
                logger.error("No static tile located at {}".format(file_static_tile))
                raise IOError("No static tile located at {}".format(file_static_tile))
            try:
                latitude = (static_tile.latitude).flatten().tolist()
                longitude = (static_tile.longitude).flatten().tolist()

                tileprocessor.create_default_prior_tile(tile_prior_file, start_date_mjd - 1, tile.get_tuple_dims(),
                                                        tile.dict_prior_rpv, tile.v_prior_aot, tile.sigma_prior_rpv,
                                                        tile.sigma_prior_aot,
                                                        tile.coeff_fine_coarse, ctx.aotFolder, latitude, longitude,
                                                        ctx.cold_start)
                default_prior_file = os.path.join(dir_prior_tiles,
                                                  "DEFAULT_" + tile_utils.get_tile_filename("PRIOR", ctx.satellite,
                                                                                            ctx.radiometer,
                                                                                            ctx.projection,
                                                                                            str(ctx.ssp).zfill(4),
                                                                                            ctx.platform, "","", tile.id))
                subprocess.check_call(['cp', tile_prior_file, default_prior_file])
            except RuntimeError:
                logger.error("Runtime error from C++ function tilemaker.create_default_prior_tile")
                raise RuntimeError('Runtime error from function tilemaker.create_default_prior_tile')


# ==================================================================================================================== #
def make_input_tiles(ctx, period_index, logger):
    """
    Create all INPUT netCDF files associated with a given accumulation period and for all active tiles (see common_init.json file). \
    Uses an instance of the C++ TileMaker class.

    .. warning:: netCDF files associated with tiles are created only if they do not exist yet (no replacement).

    .. seealso:: :func:`scheduler.main.run_scheduler`
    """
    tilemaker = initialize_tilemaker(ctx, logger)

    start_date_acc_period = ctx.v_acc_periods[period_index].start_date
    end_date_acc_period = ctx.v_acc_periods[period_index].end_date


    print(('TileMaker module starts to create input tiles for period {0} to {1}'.format(start_date_acc_period.strftime("%Y-%m-%d"),
                                                                                       end_date_acc_period.strftime("%Y-%m-%d"))))
    logger.info('TileMaker module starts to create input tiles for period {0} to {1}'.format(start_date_acc_period.strftime("%Y-%m-%d"),
                                                                                                    end_date_acc_period.strftime("%Y-%m-%d")))


    # get the list of input tiles for the given accumulation period
    list_input_tile_filenames = ctx.filenamesGenerator.generate_list_input_files(start_date_acc_period,
                                                                                 end_date_acc_period,
                                                                                 ctx.satellite, str(ctx.ssp).zfill(4),
                                                                                 ctx.projection, ctx.platform,
                                                                                 ctx.ecmwf_constant,
                                                                                 ctx.aerosol_height_constant,
                                                                                 ctx.daily_cloudmask,
                                                                                 logger)

    if len(list_input_tile_filenames) == 0:
        logger.warning('The list of input tile filenames is empty for two possible reasons.\n' 
                       '\t\t- check that all the paths to image and cloud data files are correct\n' 
                       '\t\t  (image and cloud directories should be structured as follows: path_to_folder/YYYY/MM/).\n'
                       '\t\t- no observation present for this accumulation period. Period processing disregarded.')
    # read the ERA-INTERIM file (The file contains a time dimension, that will be used for the temporal interpolation of ECMWF data)


    ## ========================= LOOP OVER SLOTS =========================
    for index, inputs in enumerate(list_input_tile_filenames):
        if index == 0 and not ctx.ecmwf_constant:
            st_date = inputs.date_and_time.strftime('%Y-%m-%d')
            st_time = inputs.date_and_time.strftime('%H-%M')
            dayofmonth = int(st_date.split('-', 2)[2])

            path_era_interim_file = ctx.filenamesGenerator.build_path_ERA_INTERIM_file(start_date_acc_period, logger)

            if inputs.slot < 12:  # assigned to ERA INTERIM slot 00:00
                ecmwf_image_slot = (dayofmonth-1) * 4
            elif inputs.slot < 24:  # assigned to ERA INTERIM slot 06:00
                ecmwf_image_slot = (dayofmonth-1) * 4 + 1
            elif inputs.slot < 36:  # assigned to ERA INTERIM slot 12:00
                ecmwf_image_slot = (dayofmonth-1) * 4 + 2
            elif inputs.slot < 48:  # assigned to ERA INTERIM slot 18:00
                ecmwf_image_slot = (dayofmonth-1) * 4 + 3

            try:
                if os.path.isfile(path_era_interim_file):
                    tilemaker.read_first_ERA_INTERIM(path_era_interim_file, ecmwf_image_slot)
                    logger.debug('end read first ERA-INTERIM file')
                else:
                    print("The path to the first ERA-INTERIM file does not exist or cannot be accessed.")
            except RuntimeError:
                logger.error(
                    "Error when trying to read the first ERA-INTERIM file {0}".format(path_era_interim_file))
                raise RuntimeError(
                    'Error when trying to read the first ERA-INTERIM file {0}'.format(path_era_interim_file))



            arr_input_tile_files, arr_input_tile_files_local, dict_input_tile_files_local = tile_utils.get_list_input_tile_files_as_dict(
                ctx, st_date, st_time)

            dict_static_tile_files = tile_utils.get_list_static_data_tile_files_as_dict(ctx)

            if len(dict_static_tile_files) == 0:
                logger.error("Empty static tile list.")
                raise ValueError("Empty static tile list.")
            try:
                if os.path.isfile(inputs.path_image_file) and os.path.isfile(
                        inputs.path_mask_file) and len(arr_input_tile_files_local) > 0:
                    if ctx.radiometer == 'MVIRI':
                        args = OrderedDict([
                            ("m_input_data_files", dict_input_tile_files_local),
                            ("date", st_date),
                            ("time", st_time),
                            ("dayOfMonth", dayofmonth),
                            ("counter", index),
                            ("slot", inputs.slot),
                            ("path_image_file", inputs.path_image_file),
                            ("path_mask_file", inputs.path_mask_file),
                            ("path_modparams_file", inputs.path_modparams_file),
                            ("path_elevation_file", inputs.path_elevation_file),
                            ("path_geopotential_file", inputs.path_geopotential_file),
                            ("path_aerosol_elevation_file_1", inputs.path_aerosol_elevation_file_1),
                            ("path_aerosol_elevation_file_2", inputs.path_aerosol_elevation_file_2),
                            ("reliability_cfc_retrieval", ctx.reliability_cfc_retrieval),
                            ("heightIsConstant", ctx.aerosol_height_constant),
                            ("useCMASKSCORE", ctx.useCMASKSCORE),
                            ("usePCFC", ctx.usePCFC),
                            ("ecmwf_constant", ctx.ecmwf_constant),
                            ("cold_start", ctx.cold_start),
                            ("daily_cloudmask", ctx.daily_cloudmask)
                        ])

                        # TODO: improve that to preserve OrderedDict ordering
                        logger.debug("tilemaker.make_input_tile arguments:\n" +
                                     pprint.pformat(dict(args)))

                        tilemaker.make_input_tile(*(item[1] for item in list(args.items())))
                    else:
                        logger.error("Radiometer not recognized")
                        raise RuntimeError("Radiometer not recognized")
                    for i in range(0, len(arr_input_tile_files)):
                        if not os.path.isfile(arr_input_tile_files_local[i]):
                            logger.info("Missing Input Tile {}, the processing of this observation is disregarded.".format(arr_input_tile_files_local[i]))
                            continue
                        if not arr_input_tile_files_local[i] == arr_input_tile_files[i]:
                            shutil.copy2(arr_input_tile_files_local[i], arr_input_tile_files[i])
                        if not filecmp.cmp(arr_input_tile_files_local[i], arr_input_tile_files[i]):
                            logger.error(
                                "Error when trying to copy input tile associated with datetime {0} "
                                "from local to destination".format(inputs.date_and_time))
                        else:
                            ##CGS 5/17: The copy was successful and the local input tiles can be removed (issue 49) this makes sure that a core dump will not occur due to insuffient memory
                            if not arr_input_tile_files_local[i] == arr_input_tile_files[i]:
                                os.remove(arr_input_tile_files_local[i])

                else:
                    if not os.path.isfile(inputs.path_image_file):
                        print_and_log("The path to the image file "
                                      "does not exist or cannot be accessed [{}]".format(inputs.path_image_file),
                                      logger.error)
                    elif not os.path.isfile(inputs.path_mask_file):
                        print_and_log("The path to the mask file "
                                      "does not exist or cannot be accessed [{}]".format(inputs.path_mask_file),
                                      logger.error)
            except RuntimeError:
                logger.error(
                    "Error when trying to create input tile associated with datetime {0}".format(inputs.date_and_time))
                raise RuntimeError(
                    'Error when trying to create input tile associated with datetime {0}'.format(inputs.date_and_time))

        else: #this else is to treat the make_input_tile in case the era interim data are considered constant
            st_date = inputs.date_and_time.strftime('%Y-%m-%d')
            st_time = inputs.date_and_time.strftime('%H-%M')
            dayofmonth = int(st_date.split('-', 2)[2])

            arr_input_tile_files, arr_input_tile_files_local, dict_input_tile_files_local = tile_utils.get_list_input_tile_files_as_dict(
                ctx, st_date, st_time)

            dict_static_tile_files = tile_utils.get_list_static_data_tile_files_as_dict(ctx)

            if len(dict_static_tile_files) == 0:
                logger.error("Empty static tile list.")
                raise ValueError("Empty static tile list.")
            try:
                if os.path.isfile(inputs.path_image_file) and os.path.isfile(
                        inputs.path_mask_file) and len(arr_input_tile_files_local) > 0:

                    if ctx.radiometer == 'MVIRI':
                        args = OrderedDict([
                            ("m_input_data_files", dict_input_tile_files_local),
                            ("date", st_date),
                            ("time", st_time),
                            ("dayOfMonth", dayofmonth),
                            ("counter", index),
                            ("slot", inputs.slot),
                            ("path_image_file", inputs.path_image_file),
                            ("path_mask_file", inputs.path_mask_file),
                            ("path_modparams_file", inputs.path_modparams_file),
                            ("path_elevation_file", inputs.path_elevation_file),
                            ("path_geopotential_file", inputs.path_geopotential_file),
                            ("path_aerosol_elevation_file_1", inputs.path_aerosol_elevation_file_1),
                            ("path_aerosol_elevation_file_2", inputs.path_aerosol_elevation_file_2),
                            ("reliability_cfc_retrieval", ctx.reliability_cfc_retrieval),
                            ("heightIsConstant", ctx.aerosol_height_constant),
                            ("useCMASKSCORE", ctx.useCMASKSCORE),
                            ("usePCFC", ctx.usePCFC),
                            ("ecmwf_constant", ctx.ecmwf_constant),
                            ("cold_start", ctx.cold_start),
                            ("daily_cloudmask", ctx.daily_cloudmask)
                        ])

                        # TODO: improve that to preserve OrderedDict ordering
                        logger.debug("tilemaker.make_input_tile arguments:\n" +
                                     pprint.pformat(dict(args)))

                        tilemaker.make_input_tile(*(item[1] for item in list(args.items())))
                    else:
                        logger.error("Radiometer not recognized")
                        raise RuntimeError("Radiometer not recognized")
                    for i in range(0, len(arr_input_tile_files)):
                        if not os.path.isfile(arr_input_tile_files_local[i]):
                            logger.info(
                                "Missing Input Tile {}, the processing of this observation is disregarded.".format(
                                    arr_input_tile_files_local[i]))
                            continue
                        if not arr_input_tile_files_local[i] == arr_input_tile_files[i]:
                            shutil.copy2(arr_input_tile_files_local[i], arr_input_tile_files[i])
                        if not filecmp.cmp(arr_input_tile_files_local[i], arr_input_tile_files[i]):
                            logger.error(
                                "Error when trying to copy input tile associated with datetime {0} "
                                "from local to destination".format(inputs.date_and_time))
                        else:
                            ##CGS 5/17: The copy was successful and the local input tiles can be removed (issue 49) this makes sure that a core dump will not occur due to insuffient memory
                            if not arr_input_tile_files_local[i] == arr_input_tile_files[i]:
                                os.remove(arr_input_tile_files_local[i])

                else:
                    if not os.path.isfile(inputs.path_image_file):
                        print_and_log("The path to the image file "
                                      "does not exist or cannot be accessed [{}]".format(inputs.path_image_file),
                                      logger.error)
                    elif not os.path.isfile(inputs.path_mask_file):
                        print_and_log("The path to the mask file "
                                      "does not exist or cannot be accessed [{}]".format(inputs.path_mask_file),
                                      logger.error)
            except RuntimeError:
                logger.error(
                    "Error when trying to create input tile associated with datetime {0}".format(inputs.date_and_time))
                raise RuntimeError(
                    'Error when trying to create input tile associated with datetime {0}'.format(inputs.date_and_time))

# ==================================================================================================================== #
def process_tiles(ctx, period_index, logger):
    """
    Process the inversion for all active tiles (see common_init.json file) for a given accumulation period. \
    Uses an instance of the C++ TileProcessor class.

    Inside the C++ method inversion of the class TileProcessor, the different steps are executed:
        #. STATIC, INPUT and PRIOR netCDF files associated with the current tile and accumulation are read.
        #. Two vectors of C structures Pixel_Input and Pixel_Output respectively are prepared \
           (Observations are filtering according to cloud mask, values of angles, ... and first guess values are generated). \
           These two structures are required to allow to pass information between the C++ and Fortran codes.
        #. The Fortran algorithm CISAR is called and perform the inversion for each pixel in parallel (OpenMP).
        #. The prior values are updated and saved in a PRIOR netCDF file
        #. The output of the inversion process is saved in a SOLUTION netCDF file.

    :param period_index: index of the accumulation period for which you want to process the inversion.
    :type period_index: int

    .. warning:: netCDF files associated with tiles are created only if they do not exist yet (no replacement).

    .. seealso:: :func:`scheduler.main.run_scheduler`
    """
    start_date_acc_period = ctx.v_acc_periods[period_index].start_date
    end_date_acc_period = ctx.v_acc_periods[period_index].end_date

    tileprocessor = initialize_tileprocessor(ctx, start_date_acc_period)
    # execute the cold_init
    # tileprocessor.cold_init()

    list_input_tile_filenames = ctx.filenamesGenerator.generate_list_input_files(start_date_acc_period,
                                                                                 end_date_acc_period,
                                                                                 ctx.satellite, str(ctx.ssp).zfill(4),
                                                                                 ctx.projection, ctx.platform,
                                                                                 ctx.ecmwf_constant,
                                                                                 ctx.aerosol_height_constant,
                                                                                 ctx.daily_cloudmask,
                                                                                 logger)

    # building metadata lists for solution tile
    list_satellite_images = []
    list_cloud_mask = []
    list_geopotential = []
    list_elevation = []
    list_aerosol_elevation=[]
    list_modparams = []
    list_auxiliary_data=[]
    for index, inputs in enumerate(list_input_tile_filenames):
        list_satellite_images.append(inputs.path_image_file)
        list_cloud_mask.append(inputs.path_mask_file)
        list_geopotential.append(inputs.path_geopotential_file)
        list_elevation.append(inputs.path_elevation_file)
        list_aerosol_elevation.append(inputs.path_aerosol_elevation_file_1)
        list_aerosol_elevation.append(inputs.path_aerosol_elevation_file_2)
        list_modparams.append(inputs.path_modparams_file)
    list_LUT = [os.path.join(ctx.dir_cisar_lut, ctx.cisar_lut_ID)]

    list_geopotential = list(set(list_geopotential))
    list_cloud_mask = list(set(list_cloud_mask))
    list_elevation = list(set(list_elevation))
    list_aerosol_elevation = list(set(list_aerosol_elevation))
    list_modparams = list(set(list_modparams))

    list_auxiliary_data = list_cloud_mask + list_geopotential +list_elevation +\
                          list_aerosol_elevation + list_modparams+ list_LUT

    ## Get the start of the next accumulation period
    ## If it's the last accumulation period, we add one day
    if period_index == len(ctx.v_acc_periods) - 1:
        start_date_next_acc_period = ctx.v_acc_periods[period_index].end_date + timedelta(days=1)
    else:
        start_date_next_acc_period = ctx.v_acc_periods[period_index + 1].start_date

        ## Check if we are inside the training period
    b_training_period = True if start_date_acc_period.year == 2007 else False

    print(('TileProcessor module starts to process tiles for period {0} to {1}'.format(start_date_acc_period.strftime("%Y-%m-%d"),
                                                                                       end_date_acc_period.strftime("%Y-%m-%d"))))
    logger.info('TileProcessor module starts to process tiles for period {0} to {1}'.format(start_date_acc_period.strftime("%Y-%m-%d"),
                                                                                       end_date_acc_period.strftime("%Y-%m-%d")))

    ## Get directories
    dir_input_tiles = os.path.join(ctx.rootdir_tile_processor, ctx.subdir_input_tiles)
    dir_prior_tiles = os.path.join(ctx.rootdir_tile_processor, ctx.subdir_prior_tiles)
    dir_solution_tiles = os.path.join(ctx.rootdir_tile_processor, ctx.subdir_solution_tiles)

    ## Get list of tiles
    v_tiles = ctx.get_list_tiles_associated_with_machine(ctx.machine_name)

    ## ========================= LOOP OVER TILES =========================
    for tile in v_tiles:

        ## Update tile dimensions + params first guess for the current Tile
        tileprocessor.prepare_to_process_new_tile(tile.get_tuple_dims(), tile.dict_prior_rpv, tile.sigma_prior_rpv,
                                                  tile.v_prior_aot, tile.sigma_prior_aot, tile.coeff_fine_coarse,
                                                  ctx.dust_prior, ctx.dust_sigma_prior, ctx.dust_first_guess,
                                                  tile.pert_first_guess, tile.min_AOT_first_guess,
                                                  tile.max_AOT_first_guess, ctx.aotFolder)

        ## filenames of prior (current + next) tiles + solution tile
        st_date = start_date_acc_period.strftime('%Y-%m-%d')
        st_time = start_date_acc_period.strftime('%H-%M')
        st_date_next = start_date_next_acc_period.strftime('%Y-%m-%d')
        st_time_next = start_date_next_acc_period.strftime('%H-%M')

        file_prior_tile = os.path.join(dir_prior_tiles,
                                       tile_utils.get_tile_filename("PRIOR", ctx.satellite, ctx.radiometer,
                                                                    ctx.projection, str(ctx.ssp).zfill(4), ctx.platform,
                                                                    st_date, st_time, tile.id))
        file_prior_tile_next = os.path.join(dir_prior_tiles,
                                            tile_utils.get_tile_filename("PRIOR", ctx.satellite, ctx.radiometer,
                                                                         ctx.projection,
                                                                         str(ctx.ssp).zfill(4), ctx.platform,
                                                                         st_date_next, st_time_next, tile.id))
        file_solution_tile = os.path.join(dir_solution_tiles,
                                          tile_utils.get_tile_filename("SOLUTION", ctx.satellite, ctx.radiometer,
                                                                       ctx.projection,
                                                                       str(ctx.ssp).zfill(4), ctx.platform,
                                                                       st_date, st_time, tile.id))



        logger.info('Prior Tile to be loaded: {}'.format(file_prior_tile))

        ## CGS 5/17: add a check to see whether the prior tile to be loaded exists, if not an error message will be logged
        if not os.path.isfile(file_prior_tile):
            logger.error("Previous prior tile does not exist. Inversion will fail.")
        logger.info('Prior Tile to be created: {}'.format(file_prior_tile_next))
        logger.info('Solution Tile to be created: {}'.format(file_solution_tile))

        ## check if the SOLUTION Tile already exists. If not, proceed to inversion.
        ## CGS 5/17: add a check to see whether the next prior tile exist, will be done at the same time as to check for the existence of the solution tile since both are created in the same c++ procedure.
        ## if one of the two does not exist, the inversion will be redone (issue 40)
        if ctx.cold_start:
            os.remove(file_solution_tile) if os.path.isfile(file_solution_tile) else None
            os.remove(file_prior_tile_next) if os.path.isfile(file_prior_tile_next) else None

        if not os.path.isfile(file_solution_tile) or not os.path.isfile(file_prior_tile_next):

            # load static tile
            static_tile_file = os.path.join(dir_input_tiles,
                                            tile_utils.get_tile_filename("STATIC", ctx.satellite, ctx.radiometer,
                                                                         ctx.projection,
                                                                         str(ctx.ssp).zfill(4), ctx.platform,
                                                                         "", "", tile.id))

            logger.info('Static Tile to be loaded: {}'.format(static_tile_file))
            ## CGS 5/17: add a check to see whether the static tile to be loaded exists, if not an error message will be logged
            if not os.path.isfile(static_tile_file):
                logger.error("Static tile " + static_tile_file + " does not exist. Inversion will fail.")

            # load input tiles
            list_input_tiles = []
            current_date = start_date_acc_period
            while current_date <= end_date_acc_period:
                st_date = current_date.strftime('%Y-%m-%d')
                filename_pattern = tile_utils.get_tile_filename("INPUT", ctx.satellite, ctx.radiometer,
                                                                ctx.projection, str(ctx.ssp).zfill(4),
                                                                ctx.platform, st_date, "*", tile.id)


                list_input_tiles.extend([file for file in glob(dir_input_tiles + filename_pattern)])
                current_date += timedelta(days=1)

            if len(list_input_tiles) == 0 :
                logger.warning('No Input Tiles associated to period from {} to {}'.format(
                    start_date_acc_period.strftime('%Y-%m-%d'), end_date_acc_period.strftime('%Y-%m-%d')))

            list_input_tiles.sort()

            list_tiles_current_acc_period = "Tiles to process :\n"
            for st_input_tile_file in list_input_tiles:
                list_tiles_current_acc_period += "\t" + st_input_tile_file + "\n"
            logger.debug(list_tiles_current_acc_period)
            # print list_tiles_current_acc_period

            start = time.clock()


            ## compute inversion
            start_acc_period_mjd = jd_to_mjd(datetime_to_jd(start_date_acc_period))
            end_acc_period_mjd = jd_to_mjd(datetime_to_jd(end_date_acc_period))
            try:
                tileprocessor.inversion(period_index, start_acc_period_mjd, end_acc_period_mjd, static_tile_file,
                                        list_input_tiles, list_satellite_images, list_auxiliary_data,
                                        file_prior_tile, file_prior_tile_next, file_solution_tile, b_training_period,
                                        ctx.line_stride, ctx.column_stride)

                ## Remove all input netCDF files that are no longer required

            except RuntimeError as e:
                logger.error(str(e))
                raise RuntimeError(str(e))

        else:
            logger.info(
                "Solution Tile " + file_solution_tile + " and Prior Tile " + file_prior_tile_next + " already exists. No inversion to perform.")


