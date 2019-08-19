"""
  The copyrights for the GEDAP algorithm and computer codes, remain with
  Rayference SPRL as an Intellectual Property Right

  Rayference Copyright (c)

"""

import os
from gedap.scheduler.utils.solution_tile import SolutionTile
from gedap.scheduler.utils.static_data_tile import StaticDataTile
from gedap.scheduler.utils.tile import get_tile_filename
from fiduceo.cdr.writer.cdr_writer import CDRWriter
from datetime import datetime, timedelta
import numpy

def generate_products(ctx, period_index, logger):
#    numpy.seterr(divide='ignore', invalid='ignore')
    if not os.path.isdir(ctx.dir_products): os.makedirs(ctx.dir_products)

    start_date_acc_period = ctx.v_acc_periods[period_index].start_date
    end_date_acc_period = start_date_acc_period + timedelta(minutes=23*60+59)

    st_date = start_date_acc_period.strftime('%Y-%m-%d')
    st_time = start_date_acc_period.strftime('%H-%M')
        
    fiduceo_product_albedo = CDRWriter.createTemplate('ALBEDO', 5000, 5000)
    fiduceo_product_aot = CDRWriter.createTemplate('AOT', 5000, 5000)
    
    for tile in ctx.v_tiles:
        start_line = tile.first_line
        end_line = tile.last_line
        start_column = tile.first_column
        end_column = tile.last_column
        root_dir = ctx.rootdir_tile_maker[tile.machine]
        dir_input_tiles = os.path.join(root_dir, ctx.subdir_input_tiles)
        dir_solution_tiles = os.path.join(root_dir, ctx.subdir_solution_tiles)
        
        file_static_data_tile = os.path.join(dir_input_tiles,
                                             get_tile_filename("STATIC", ctx.satellite, ctx.radiometer,
                                                               ctx.projection,
                                                               str(ctx.ssp).zfill(4), ctx.platform,
                                                               "", "", tile.id))
        file_solution_tile = os.path.join(dir_solution_tiles,
                                          get_tile_filename("SOLUTION", ctx.satellite, ctx.radiometer,
                                                            ctx.projection,
                                                            str(ctx.ssp).zfill(4), ctx.platform,
                                                            st_date, st_time, tile.id))

        if os.path.isfile(file_static_data_tile) and os.path.isfile(file_solution_tile):

            # ============= Extract data from static tile =============
            static_data_tile = StaticDataTile(file_static_data_tile)
            
            fiduceo_product_aot.latitude[start_line:end_line+1, start_column:end_column+1] = \
                static_data_tile.latitude[:, :]
            fiduceo_product_aot.longitude[start_line:end_line + 1, start_column:end_column + 1] = \
                static_data_tile.longitude[:, :]

            # ============= AOT CDR  =============
            
            solution_tile = SolutionTile(file_solution_tile)

            fiduceo_product_aot.attrs['source'] = solution_tile.source
            fiduceo_product_aot.attrs['auxiliary_data'] = solution_tile.auxiliary_data

            if solution_tile.acquisition_time.shape[0] == 1:
                # time has to be expressed in seconds from 1st January 1970
                t0 = datetime(1970, 1, 1)
                st_date_noon = start_date_acc_period + timedelta(hours=12) 
                delta = st_date_noon - t0
                fiduceo_product_aot.time[start_line:end_line + 1, ] = t0 + timedelta(seconds=delta.total_seconds())
                fiduceo_product_albedo.time[start_line:end_line + 1, ] = t0 + timedelta(seconds=delta.total_seconds())
            
                fiduceo_product_aot.aot[start_line:end_line + 1, start_column:end_column + 1] = numpy.where(solution_tile.total_tau[0, 0, :, :] < 0.01, numpy.nan, solution_tile.total_tau[0, 0, :, :])
                aot = numpy.where(solution_tile.total_tau[0, 0, :, :] < 0.01, numpy.nan, solution_tile.total_tau[0, 0, :, :])
                sigma_aot = numpy.where(solution_tile.total_tau[0, 0, :, :] < 0.01, numpy.nan, solution_tile.sigma_total_tau[0, 0, :, :])
                fiduceo_product_aot.u_independent_aot[start_line:end_line + 1, start_column:end_column + 1] = sigma_aot/aot
#                fiduceo_product_aot.u_independent_aot[start_line:end_line + 1, start_column:end_column + 1] = numpy.where(solution_tile.total_tau[0, 0, :, :] < 0.01, numpy.nan, numpy.divide(solution_tile.sigma_total_tau[0, 0, :, :], solution_tile.total_tau[0, 0, :, :]))
                
                not_enough_obs = solution_tile.error_flag[:,:] == 202
                quality_cut_aot = solution_tile.total_tau[0, 0, :, :] < 0.01
                mask_aot = numpy.logical_or(not_enough_obs, quality_cut_aot)

                fiduceo_product_aot.quality_pixel_bitmask[start_line:end_line + 1, start_column:end_column + 1] = numpy.where(mask_aot, 1, 0)

                # ============= ALBEDO CDR  =============

                fiduceo_product_albedo.attrs['source'] = solution_tile.source
                fiduceo_product_albedo.attrs['auxiliary_data'] = solution_tile.auxiliary_data

                quality_cut_albedo = solution_tile.BHR[0, :, :] < 0.001
                land_cover = static_data_tile.land_cover[:,:] == 2.0
                
                mask_albedo = numpy.logical_or(quality_cut_albedo, land_cover)

                fiduceo_product_albedo.surface_albedo[start_line:end_line + 1, start_column:end_column + 1] = numpy.where(mask_albedo, numpy.nan, solution_tile.BHR[0, :, :])

                sigma_BHR = numpy.where(mask_albedo, numpy.nan, solution_tile.sigma_BHR[0, :, :])
                BHR = numpy.where(mask_albedo, numpy.nan, solution_tile.BHR[0, :, :])
                fiduceo_product_albedo.u_independent_surface_albedo[start_line:end_line + 1, start_column:end_column + 1] = sigma_BHR/BHR
                
                fiduceo_product_albedo.quality_pixel_bitmask[start_line:end_line + 1, start_column:end_column + 1] = numpy.where(mask_albedo, 1, 0)
        
        else:
            raise SystemExit('No Solution Tile and Static Tile associated to {}'.format(tile.id))

    fiduceo_product_aot.attrs['institution'] = 'Rayference SPRL'
    fiduceo_product_aot.attrs['title'] = 'MVIRI Aerosol FIDUCEO CDR'
    fiduceo_product_aot.attrs['history'] = "Created: "+ datetime.now().strftime("%a %b %H:%M:%S %Y")
    fiduceo_product_aot.attrs['references'] = ''
    fiduceo_product_aot.attrs['comment'] = ''
    fiduceo_product_aot.attrs['configuration'] = ''
    fiduceo_product_aot.attrs['time_coverage_start'] = start_date_acc_period.strftime("%Y%m%dT%H%M%SZ")
    fiduceo_product_aot.attrs['time_coverage_end'] = end_date_acc_period.strftime("%Y%m%dT%H%M%SZ")
    fiduceo_product_aot.attrs['time_coverage_duration'] = 'P1D'
    fiduceo_product_aot.attrs['time_coverage_resolution'] = 'P1D'


    platform = ctx.platform + '-' + ctx.ssp_cdr
    filename_aot_product = CDRWriter.create_file_name_CDR('AOT', ctx.radiometer, platform, start_date_acc_period, end_date_acc_period, 'L2', str(ctx.gedap_version))
    CDRWriter.write(fiduceo_product_aot, os.path.join(ctx.dir_products, filename_aot_product))                  


    fiduceo_product_albedo.attrs['institution'] = 'Rayference SPRL'
    fiduceo_product_albedo.attrs['title'] = 'MVIRI Albedo FIDUCEO CDR'
    fiduceo_product_albedo.attrs['history'] = 'Created: '+ datetime.now().strftime("%a %b %H:%M:%S %Y")
    fiduceo_product_albedo.attrs['references'] = ''
    fiduceo_product_albedo.attrs['comment'] = ''
    fiduceo_product_albedo.attrs['configuration'] = ''
    fiduceo_product_albedo.attrs['time_coverage_start'] = start_date_acc_period.strftime("%Y%m%dT%H%M%SZ")
    fiduceo_product_albedo.attrs['time_coverage_end'] = end_date_acc_period.strftime("%Y%m%dT%H%M%SZ")
    fiduceo_product_albedo.attrs['time_coverage_duration'] = 'P1D'
    fiduceo_product_albedo.attrs['time_coverage_resolution'] = 'P1D'


    filename_albedo_product = CDRWriter.create_file_name_CDR('ALBEDO', ctx.radiometer, platform, start_date_acc_period, end_date_acc_period, 'L2', str(ctx.gedap_version))
    CDRWriter.write(fiduceo_product_albedo, os.path.join(ctx.dir_products, filename_albedo_product))
