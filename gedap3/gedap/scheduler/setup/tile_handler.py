"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
from ...cpp import MFGTileMaker, TileProcessor
from ..report.logger import cpp_log
import os

from datetime import datetime, timedelta
from glob import glob

def initialize_tilemaker(ctx, logger):
    # Create and prepare the TileMaker
    v_tiles_as_tuples = [tile.get_tuple_dims() for tile in ctx.v_tiles]
    tilemaker = None

    if ctx.radiometer == 'MVIRI':
        tilemaker = MFGTileMaker(ctx.v_st_bands, v_tiles_as_tuples, cpp_log)
        logger.info('MFGTileMaker created')
    else:
        logger.error('Radiometer not recognized, a new TileMaker should be implemented')
    return tilemaker


def initialize_tileprocessor(ctx, start_date_acc_period):
    #finding closes LUT in time
    lut_prefix = os.path.join(ctx.lut_ID_prefix + '_' + ctx.platform + '_' + ctx.radiometer + \
                              '_' + ctx.v_st_bands[0] + '_')
    all_luts = []

    if str(start_date_acc_period.month).zfill(2) == '01':
        year_list = [start_date_acc_period.year, start_date_acc_period.year - 1]
        for year in year_list:
            lut_list = glob(os.path.join(ctx.dir_cisar_lut, lut_prefix + str(year)) + '*')
            all_luts.extend(lut_list)
    elif start_date_acc_period.month == '12':
        year_list = [start_date_acc_period.year, start_date_acc_period.year + 1]
        for year in year_list:
            lut_list = glob(os.path.join(ctx.dir_cisar_lut, lut_prefix + year + '*'))
            all_luts.extend(lut_list)
    else:
        lut_list = glob(os.path.join(ctx.dir_cisar_lut, lut_prefix + str(start_date_acc_period.year)) + '*')
        all_luts = lut_list

    available_dates = []
    indices = []

    for index, lut in enumerate(all_luts):
        date = lut.split('_')[-1].split('.')[0]
        year = date[0:4]
        starting_date = datetime(int(year), 1, 1)
        dayoftheyear = date[4:]
        actual_date = starting_date + timedelta(days=int(dayoftheyear))
        available_dates.append(actual_date)
        indices.append([actual_date, date])

    def nearest(items, pivot):
        return min(items, key=lambda x: abs(x - pivot))

    closest_date = nearest(available_dates, start_date_acc_period)

    for i in indices:
        if i[0] == closest_date:
            lut_date = i[1]

    ctx.cisar_lut_ID = os.path.join(ctx.lut_ID_prefix + '_' + ctx.platform + '_' + ctx.radiometer + \
                                     '_' + ctx.v_st_bands[0] + '_' + lut_date)

    tileprocessor = TileProcessor(ctx.dir_cisar_lut, ctx.cisar_lut_ID, ctx.cisar_aeromod_filename,
                                  ctx.max_nb_iterations,
                                  ctx.acc_period_length, ctx.v_st_bands, ctx.dict_aerosol_classes,
                                  ctx.b_temporal_smoothness, ctx.b_spectral_constraint,
                                  ctx.b_surface_spectral_constraint, cpp_log)
    return tileprocessor


    


