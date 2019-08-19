"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
import os

# ==================================================================================================================== #
#                                             UTIL FUNCTIONS                                                           #
# ==================================================================================================================== #


# function to generate Tile filename
def get_tile_filename(st_type, st_satellite, st_radiometer, st_projection, st_ssp, st_platform, st_date, st_time, st_area):
    """
    Generate the filename of netCDF file used internally by GEDAP.

    .. warning:: all Python and C++ GEDAP code assumes that filenames of tile netCDF files follow a certain pattern. \
                 Be very careful if you want to modify this function.

    :param st_type: type of the netCDF file (STATIC, INPUT, PRIOR, SOLUTION)
    :type st_type: string
    :param st_satellite: name of the satellite
    :type st_satellite: string
    :param st_radiometer: name of the radiometer
    :type st_radiometer: string
    :param st_projection: name of the projection
    :type st_projection: string
    :param st_date: date as string (should be "" for STATIC netCDF files)
    :type st_date: string
    :param st_time: time as string (should be "" for STATIC netCDF files)
    :type st_time: string
    :param st_area: area associated with the tile (tile ID).
    :type st_area: string
    :return: filename of netCDF file associated with a tile.
    :rtype: string
    """

    tile_filename = st_type + "_TILE_" + st_satellite + "_" + st_radiometer + "_"  + st_projection + "_" + st_ssp

    if st_date != "" and st_time != "" and st_type != "STATIC":
        tile_filename +=  "_" + st_platform + "_DATETIME_" + st_date + "_" + st_time
    if st_date == "" and st_time == "" and st_type == "PRIOR":
        tile_filename +=  "_" + st_platform

    tile_filename += "_AREA_" + st_area + ".nc"
    return tile_filename


def get_list_static_data_tile_files_as_dict(ctx):
    """
    Generate the list of filenames of STATIC netCDF files. \
    This list is passed to the C++ function make_static_data_tile of the class TileMaker.

    :param ctx: instance of the class :class:`scheduler.context.Context`
    :return: dictionary with items (tile ID, filename of associated STATIC netCDF file).
    :rtype: dictionary[string,string]

    .. seealso:: :func:`get_tile_filename`
    """
    dict_static_tile_files = dict()
    for tile in ctx.v_tiles:
        dir_input_tiles = os.path.join(ctx.rootdir_tile_maker[tile.machine], ctx.subdir_input_tiles)
        static_tile_file = os.path.join(dir_input_tiles, get_tile_filename("STATIC", ctx.satellite, ctx.radiometer,
                                                                           ctx.projection,
                                                                           str(ctx.ssp).zfill(4), ctx.platform,
                                                                           "", "", tile.id))
        dict_static_tile_files[tile.id] = static_tile_file
    return dict_static_tile_files


def get_list_input_tile_files_as_dict(ctx, st_date, st_time):
    """
    Generate the list of filenames of INPUT netCDF files associated with a given date and time. \
    This list is passed to the C++ function make_input_tile of the class TileMaker.

    :param ctx: instance of the class :class:`scheduler.context.Context`
    :param st_date: date as string
    :param st_time: time as string
    :return: dictionary with items (tile ID, filename of associated INPUT netCDF file).
    :rtype: dictionary[string,string]

    .. seealso:: :func:`get_tile_filename`
    """

    arr_input_tile_files = []
    dict_input_tile_files_local = dict()
    arr_input_tile_files_local = []
    for tile in ctx.v_tiles:
        dir_input_tiles = os.path.join(ctx.rootdir_tile_maker[tile.machine], ctx.subdir_input_tiles)
        dir_input_tiles_local = os.path.join(ctx.rootdir_tile_maker["master"], ctx.subdir_input_tiles)
        input_tile_file = os.path.join(dir_input_tiles, get_tile_filename("INPUT", ctx.satellite,
                                                                          ctx.radiometer, ctx.projection,
                                                                          str(ctx.ssp).zfill(4), ctx.platform,
                                                                          st_date, st_time, tile.id))
        input_tile_file_local = os.path.join(dir_input_tiles_local, get_tile_filename("INPUT", ctx.satellite,
                                                                                      ctx.radiometer,
                                                                                      ctx.projection,
                                                                                      str(ctx.ssp).zfill(4),
                                                                                      ctx.platform,
                                                                                      st_date, st_time, tile.id))
        arr_input_tile_files.append(input_tile_file)
        dict_input_tile_files_local[tile.id] = input_tile_file_local
        arr_input_tile_files_local.append(input_tile_file_local)
    return arr_input_tile_files, arr_input_tile_files_local, dict_input_tile_files_local
