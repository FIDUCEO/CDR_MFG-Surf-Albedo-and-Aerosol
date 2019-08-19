import os
import numpy as np
from glob import glob
import netCDF4
from datetime import datetime, timedelta

import matplotlib.lines as mlines
import matplotlib.dates as md
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 2
plt.rc('lines', lw=2)


def monitoring(ctx, logger):
    platform = ctx.platform + '-' + ctx.ssp_cdr
    surface_albedo_product_file_list = glob(os.path.join(ctx.dir_products, 'FIDUCEO_CDR_ALBEDO_' + ctx.radiometer + '_'+ platform + '*'))
    aot_product_file_list = glob(os.path.join(ctx.dir_products, 'FIDUCEO_CDR_AOT_'+ ctx.radiometer + '_' + platform + '*'))

    dates_aot = []
    dates_albedo = []

    albedo_q1 = []
    albedo_q2 = []
    albedo_q3 = []
    albedo_mean = []

    aot_q1 = []
    aot_q2 = []
    aot_q3 = []
    aot_mean = []

    for file in surface_albedo_product_file_list:
        dates_albedo.append(datetime.strptime(file.split('/')[-1].split('_')[5][:8], "%Y%m%d"))

        try:
            fiduceo_product = netCDF4.Dataset(file)

            surface_albedo = fiduceo_product.variables['surface_albedo'][:]
            flatten_array_albedo = np.ndarray.flatten(surface_albedo)

            number_valid_pixels_albedo = np.count_nonzero(~np.isnan(flatten_array_albedo))
            number_not_valid_pixels_albedo = np.count_nonzero(np.isnan(flatten_array_albedo))
            
            x_albedo = flatten_array_albedo[~np.isnan(flatten_array_albedo)]
            
            albedo_q1.append(np.percentile(x_albedo, 25))
            albedo_q2.append(np.percentile(x_albedo, 50))
            albedo_q3.append(np.percentile(x_albedo, 75))
            albedo_mean.append(np.mean(x_albedo))

        except RuntimeError:
            logger.error("Error when trying to open FIDUCEO product {0}".format(file))
            raise RuntimeError('Error when trying to open FIDUCEO product {0}'.format(file))



    for file in aot_product_file_list:
        dates_aot.append(datetime.strptime(file.split('/')[-1].split('_')[5][:8], "%Y%m%d"))

        try:
            fiduceo_product = netCDF4.Dataset(file)
            
            nb_lines = fiduceo_product.dimensions['y'].size
            nb_columns = fiduceo_product.dimensions['x'].size

            aot = fiduceo_product.variables['aot'][:]
            flatten_array_aot = np.ndarray.flatten(aot)

            number_valid_pixels_aot = np.count_nonzero(~np.isnan(flatten_array_aot))
            number_not_valid_pixels_aot = np.count_nonzero(np.isnan(flatten_array_aot))

            x_aot = flatten_array_aot[~np.isnan(flatten_array_aot)]

            aot_q1.append(np.percentile(x_aot, 25))
            aot_q2.append(np.percentile(x_aot, 50))
            aot_q3.append(np.percentile(x_aot, 75))
            aot_mean.append(np.mean(x_aot))

        except RuntimeError:
            logger.error("Error when trying to open FIDUCEO product {0}".format(file))
            raise RuntimeError('Error when trying to open FIDUCEO product {0}'.format(file))

    if dates_albedo == dates_aot:
        # need one set of x only, had troubles with DateFormatter
        dates = dates_albedo
        nb_total_pixels = nb_lines * nb_columns
        # First plot displaying the ratio between nb_pixels processed and not
        # this differs for aot and surface albedo as the sea is not processed in surface albedo retrieval

        fig, ax = plt.subplots(figsize=[15,10])
        plt.subplots_adjust(bottom=0.2)
        plt.xticks( rotation=25)
        
        ax.plot(dates, (number_valid_pixels_aot/nb_total_pixels), 'go')
        ax.plot(dates, (number_valid_pixels_albedo/nb_total_pixels), 'gs')

        ax.set_ylabel('ratio', color='blue', fontsize =22)
        ax.tick_params('y', colors='blue',size = 5)
        ax.xaxis.set_ticks(dates)
        ax.set_xticklabels(dates,size=20)
        ax.xaxis.set_major_formatter(md.DateFormatter('%Y-%m-%d'))
        
        round = mlines.Line2D([], [], color='black', marker='o', linestyle='None', fillstyle='none', markersize=10, label='AOT ratio')
        square = mlines.Line2D([], [], color='black', marker='s', linestyle='None', fillstyle='none', markersize=10, label='ALBEDO ratio')

        plt.legend(handles=[round, square])
        plt.title("Timeseries of the ratio of processed pixels over total pixel nb")
        plt.grid(linestyle='dashed')
        plt.yticks(fontsize = 15)


        plot_dir = os.path.join(ctx.rootdir_tile_maker['master'], 'Plots', 'monitoring')

        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        plt.savefig(os.path.join(plot_dir,  'pixel_ratio_timeseries.png'), dpi=96, bbox_inches='tight')


        fig, ax1 = plt.subplots(figsize=[15,10])
        plt.subplots_adjust(bottom=0.2)
        plt.xticks( rotation=25)

        ax1.plot(dates, albedo_q1, 'bo', dates, albedo_q2,'bs', dates, albedo_q3, 'b^', dates, albedo_mean, 'b*', fillstyle = 'none')
        ax1.set_ylabel('albedo', color='blue', fontsize =22)    
        ax1.tick_params('y', colors='blue',size = 5)
        #ax1.xaxis.set_ticks(dates)
        #ax1.set_xticklabels(dates,size=20)
        #ax1.xaxis.set_major_formatter(md.DateFormatter('%Y-%m-%d'))

        plt.yticks(fontsize = 15) 
        
        ax2 = ax1.twinx()
        ax2.plot(dates, aot_q1, 'ro', dates_aot, aot_q2,'rs', dates, aot_q3, 'r^', dates, aot_mean, 'r*', fillstyle = 'none')
        ax2.set_ylabel('aot', color='red', fontsize =22)
        ax2.tick_params('y', colors='red', size = 5)
        
        ax2.xaxis.set_ticks(dates)
        ax2.set_xticklabels(dates,size=20)
        ax2.xaxis.set_major_formatter(md.DateFormatter('%Y-%m-%d'))

        round = mlines.Line2D([], [], color='black', marker='o', linestyle='None', fillstyle='none', markersize=10, label='25% Quantile')
        square = mlines.Line2D([], [], color='black', marker='s', linestyle='None', fillstyle='none', markersize=10, label='50% Quantile')
        triangle = mlines.Line2D([], [], color='black', marker='^', linestyle='None', fillstyle='none', markersize=10, label='75% Quantile')
        star = mlines.Line2D([], [], color='black', marker='*', linestyle='None', fillstyle='none', markersize=10, label='Mean')
        
        plt.legend(handles=[round, square, triangle, star])

        plt.title("Timeseries of different statistics for AOT and ALBEDO CDRs ")

        plt.yticks(fontsize = 15) 
        plt.grid(linestyle='dashed')
        plt.savefig(os.path.join(plot_dir,  'statistics_timeseries.png'), dpi=96, bbox_inches='tight')

        plt.show()

