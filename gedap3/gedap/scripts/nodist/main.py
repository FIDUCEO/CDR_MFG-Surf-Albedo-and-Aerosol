"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
import argparse
import os
import time
import resource

from numpy import arange

from gedap.scheduler.report.logger import logger
from gedap.scheduler.setup.context import Context
from gedap.scheduler.tasks.tasks import make_input_tiles, process_tiles
from gedap.scheduler.tasks.tasks import prepare_to_process_first_acc_period, update_accumulation_periods_list
from gedap.scheduler.tasks.monitoring import monitoring
from gedap.scheduler.tasks.generate_products_fiduceo import generate_products
from gedap import __version__ as gedap_version

def run_scheduler(ctx):
    """
    This function defines the sequence of tasks to run GEDAP
    """
    
    if ctx.task_ID == "all_processing":
        
        #prepare_to_process_first_acc_period(ctx, logger)
        nb_acc_periods = len(ctx.v_acc_periods)

        ##################################################################################### #
        #                             NON DISTRIBUTED TILE PROCESSING                           #
        # ##################################################################################### #
        for i_acc_per in arange(ctx.index_first_acc_period_to_process, nb_acc_periods):
            # print('current accumulation period: ', ctx.v_acc_periods[i_acc_per])
            # print('[make_input_tiles]     i_acc_per: ', i_acc_per)
            # make_input_tiles(ctx, i_acc_per, logger)
            # print('[process_tiles]     i_acc_per: ', i_acc_per)
            # process_tiles(ctx, i_acc_per, logger)
            print('[generate_products]     i_acc_per: ', i_acc_per)
            generate_products(ctx, i_acc_per, logger)
            # print('time elapsed', time.clock() - start)
            # update_accumulation_periods_list(ctx, i_acc_per)

        print('finished')

    if ctx.task_ID == "run_prepare_to_process_first_acc_period":
        prepare_to_process_first_acc_period(ctx, logger)
    elif ctx.task_ID ==  "run_input_tile_maker":
        i_acc_per = ctx.period_index
        make_input_tiles(ctx, i_acc_per, logger)
    elif ctx.task_ID == "run_tile_processor":
        i_acc_per = ctx.period_index
        process_tiles(ctx, i_acc_per, logger)
    elif ctx.task_ID == "run_product_generator":
        i_acc_per = ctx.period_index
        generate_products(ctx, i_acc_per, logger)
    elif ctx.task_ID == "update_acc_period_list":
        i_acc_per = ctx.period_index
        update_accumulation_periods_list(ctx, i_acc_per)
    elif ctx.task_ID == "monitoring":
        monitoring(ctx, logger)

# ================================================================================================================== #
#                                                   MAIN                                                             #
# ================================================================================================================== #
def main(args):

    # Initialize Context: read setup file (common_init) and perform sanity checks. If any error occurs in this first
    # part, GEDAP won't be able to run. Fatal error issued.
    gedap_env = os.getenv('GEDAP_DIR')
    resources_dir = None
    if args.resources_dir is None:
        resources_dir = os.path.abspath(os.path.join(gedap_env, 'resources/templates/mviri/'))
    else:
        resources_dir = os.path.abspath(args.resources_dir)
    logger.info("Using resources directory {}".format(resources_dir))
    ctx = Context(resources_dir)

    # Process remaining command line arguments and issue warning for cold start if not in batch mode
    if args.batch:  # In batch mode, -c flag has priority over common_init.json setting
        ctx.cold_start = args.cold_start
    else:
        if args.cold_start or ctx.cold_start:  # If cold start is requested in interactive mode, ask for confirmation
            cs = input('You are about to perform a COLD_START. This means that all the tiles will be overwritten, '
                'as well as the products. Is this what you really want to do? (y/N) ')
            ctx.cold_start = True if cs in ['y', 'Y', 'yes', 'Yes', 'on', 'On', 'true', 'True'] \
                else False

    ctx.v_acc_periods = []
    if not os.path.isfile(ctx.list_acc_periods_file) or ctx.cold_start:
        ctx.build_list_accumulation_periods()
    else:
        ctx.load_list_accumulation_periods_from_json_file()
        ctx.period_list_check()
    # get index of first period to process

    #: index of the first accumulation to process. See method :func:`get_index_first_acc_period_to_process`.
    ctx.index_first_acc_period_to_process = ctx.get_index_first_acc_period_to_process(logger)

    if args.all:
        ctx.task_ID = "all_processing"
    if args.process_first:
        ctx.task_ID = "run_prepare_to_process_first_acc_period"
    if args.make_input:
        ctx.task_ID = "run_input_tile_maker"
        ctx.period_index = int(args.period_index)
    if args.tile_processor:
        ctx.task_ID = "run_tile_processor"
        ctx.period_index = int(args.period_index)
    if args.product_generator:
        ctx.task_ID = "run_product_generator"
        ctx.period_index = int(args.period_index)
    if args.update_list:
        ctx.task_ID = "update_acc_period_list"
        ctx.period_index = int(args.period_index)
    if args.tile_index:
        ctx.idx_tile = args.tile_index
    if args.monitoring:
        ctx.task_ID = "monitoring"
    
    ctx.get_v_tile()
    
    run_scheduler(ctx)
    
    # prepare_to_process_first_acc_period(ctx, logger)
    # generate_products(ctx, 0, logger)


def main_cli():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="""GEDAP non-distributed application.""")
    parser.add_argument("-b", "--batch", action="store_true",
                        help="execute in batch mode (non-interactive)")
    parser.add_argument("-c", "--cold_start", action="store_true",
                        help="in batch mode, require cold start")
    parser.add_argument("-r", "--resources_dir", default=None,
                        help="path to resources directory to be used (defaults to $GEDAP_DIR/resources)")

    ########################################## parallel processing EUMETSAT #############################################
    parser.add_argument("-a", "--all", action="store_true",
                        help="pass this argument to perform the non distributed processing on one core")
    parser.add_argument("-f", "--process_first", action="store_true",
                        help="pass this argument to perform only the prepare_to_process_first_acc_period task")
    parser.add_argument("-m", "--make_input",action="store_true",
                        help="pass this argument to perform only the make_input_tiles task")
    parser.add_argument("-p", "--tile_processor", action="store_true",
                        help="pass this argument to perform only the process_tiles task")
    parser.add_argument("-g", "--product_generator", action="store_true",
                        help="Pass this argument to perform only the generate_products task")
    parser.add_argument("-u", "--update_list", action="store_true",
                        help="Pass this argument to update the list of accumulation periods")
    parser.add_argument("-i", "--period_index", default=None, 
                        help="Passing the index of the accumulation period to be processed")
    parser.add_argument("-t", "--tile_index", default=None, 
                        help="Passing the index of the tile to be processed")
    parser.add_argument("-n", "--monitoring", action='store_true',
                        help="Launch the monitoring script, to plot product statistics timeseries.")

    ####################################################################################################################


    parser.add_argument("--version", action="store_true",
                        help="display version information and exit")
    args = parser.parse_args()

    if args.version:
        print(("GEDAP {}".format(gedap_version)))
        exit(0)

    main(args)


if __name__ == '__main__':

    main_cli()
