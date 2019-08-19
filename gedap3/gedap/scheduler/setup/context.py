"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
"""
module context
--------------

:author: Alix Damman
:date: Created on 13 Apr 2016

This module contains the class Context and some other useful classes.
Reformatted by Vincent Leroy (18 Oct 2017)
"""


import os
import json
import numpy as np

from datetime import datetime, timedelta
from collections import OrderedDict

from gedap.scheduler.report.logger import logger
from gedap.scheduler.utils.read_col import read_col
from gedap.version import __version__ as gedap_version
from .MFG_filenames_generator import MFGFilenamesGenerator

from glob import glob

class AccumulationPeriod(object):
    """
    This class is used to define an accumulation period.

    :param start_date: starting date of the accumulation period
    :type start_date: datetime
    :param end_date: ending date of the accumulation period
    :type end_date: datetime
    :param b_processed: flag telling if the accumulation period had been already processed
    :type b_processed: boolean
    """

    def __init__(self, start_date, end_date, b_processed=False):
        self.start_date = start_date
        self.end_date = end_date
        self.b_processed = b_processed

    def __repr__(self):
        return 'start_date: ' + self.start_date.strftime('%Y-%m-%d') + \
               ', end_date: ' + self.end_date.strftime('%Y-%m-%d') + \
               ', already processed?: ' + str(self.b_processed)


class TileProperties(object):
    """
    This class contains all the informations setup associated with a geographical tile.

    :param id: ID of the tile as defined in JSON file common_init.json in /resources folder
    :type id: string
    :param machine: name of the machine responsible to process this machine
    :type machine: string
    :param first_line: index of the first line of the tile in the satellite image
    :type first_line: int
    :param last_line: index of the last line of the tile in the satellite image
    :type last_line: int
    :param first_column: index of the first column of the tile in the satellite image
    :type first_column: int
    :param last_column: index of the last column of the tile in the satellite image
    :type last_column: int
    :param st_region: name of the region associated with the tile (possible regions are: \
                      Euro (Europe), SAfr (South Africa), NAfr (North Africa), SAme (South America), 0 (Default))
    :type st_region: string
    :param dict_prior_rpv: dictionary of prior values associated with surface parameters \
                           (see C++ documentation/code of PriorManager class)
    :type dict_prior_rpv: dictionary[string,float]
    :param sigma_prior_rpv: sigma prior values associated with surface parameters \
                            (see C++ documentation/code of PriorManager class)
    :type sigma_prior_rpv: float
    :param v_prior_aot: list of prior values associated with AOT \
                        (see C++ documentation/code of PriorManager class)
    :type v_prior_aot: list[float]
    :param sigma_prior_aot: sigma prior values associated with AOT \
                            (see C++ documentation/code of PriorManager class)
    :type sigma_prior_aot: float
    :param coeff_fine_coarse: relative weight between fine and coarse classes \
                              (see C++ documentation/code of PriorManager class)
    :type coeff_fine_coarse: float
    :param pert_first_guess: perturbation factor used to build first guess values of surface parameters from prior \
                             (see C++ documentation/code of TileProcessor class)
    :type pert_first_guess: float
    :param min_AOT_first_guess: lower first guess value associated with AOT \
                                (see C++ documentation/code of TileProcessor class)
    :type min_AOT_first_guess: float
    :param max_AOT_first_guess: upper first guess value associated with AOT \
                                (see C++ documentation/code of TileProcessor class)
    :type max_AOT_first_guess: float
    """

    def __init__(self, id, machine, first_line, last_line, first_column, last_column,
                 dict_prior_rpv, sigma_prior_rpv, v_prior_aot, sigma_prior_aot, coeff_fine_coarse,
                 pert_first_guess, min_AOT_first_guess, max_AOT_first_guess):
        self.id = id
        self.machine = machine

        self.first_line = first_line
        self.last_line = last_line
        self.first_column = first_column
        self.last_column = last_column

        self.dict_prior_rpv = dict_prior_rpv
        self.sigma_prior_rpv = sigma_prior_rpv

        self.v_prior_aot = v_prior_aot
        self.sigma_prior_aot = sigma_prior_aot
        self.coeff_fine_coarse = coeff_fine_coarse

        self.pert_first_guess = pert_first_guess
        self.min_AOT_first_guess = min_AOT_first_guess
        self.max_AOT_first_guess = max_AOT_first_guess

    def __repr__(self):
        st_repr = 'id: ' + self.id + '\n'
        st_repr += 'machine (to process the tile): ' + self.machine + '\n'
        st_repr += 'dimensions: [' + \
                   str(self.first_line) + ':' + str(self.last_line) + ' , ' + \
                   str(self.first_column) + ':' + str(self.last_column) + ']' + '\n'
        st_repr += 'Prior RPV:       ' + str(self.dict_prior_rpv) + '\n'
        st_repr += 'Sigma Prior RPV: ' + str(self.sigma_prior_rpv) + '\n'
        st_repr += 'Prior AOT:       ' + str(self.v_prior_aot) + '\n'
        st_repr += 'Sigma Prior AOT: ' + str(self.sigma_prior_aot) + '\n'
        st_repr += 'Coefficient fine/coarse: ' + str(self.coeff_fine_coarse) + '\n'
        st_repr += 'Pert. First Guess: ' + str(self.pert_first_guess) + '\n'
        st_repr += 'Min First Guess:   ' + str(self.min_AOT_first_guess) + '\n'
        st_repr += 'Max First Guess:   ' + str(self.max_AOT_first_guess) + '\n'
        return st_repr
        
    def get_tuple_dims(self):
        """
        This method is used to communicate with C++ code.

        :return: tuple containing the ID and dimensions (first/last line/column).
        :rtype: tuple
        """
        
        return self.id, self.first_line, self.last_line, self.first_column, self.last_column
    
    def get_nb_lines_columns_and_pixels(self):
        """
        :returns: number of lines, columns and pixels of the tile.
        :rtype: int,int,int
        """
        
        nb_lines = np.abs(self.last_line - self.first_line) + 1
        nb_columns = np.abs(self.last_column - self.first_column) + 1
        nb_pixels = nb_lines * nb_columns

        return nb_lines, nb_columns, nb_pixels


class Context(object):
    """
    This class reads all setup information from the JSON file common_init.json  \
    generates and manages the list of accumulation periods and tiles, creates an instance XXX_filename_generator
    reponsible for generating all the filenames for the external data for a given satellite or radiometer XXX, and \
    select the index of the first accumulation to process (to easily restart GEDAP after a crash).

    :param dir_resources: path to the directory containing the JSON files common_init.json
                          (defaults to current working directory)
    :type dir_resources: string
    """
    
    def __init__(self, dir_resources='.'):

        self.dir_resources = dir_resources
        self.list_acc_periods_file = os.path.join(self.dir_resources, 'list_accumulation_periods.json')

        # read configuration file common_init.json
        try:
            dict_header = json.load(open(os.path.join(dir_resources, 'common_init.json'), 'r'))
        except ValueError as e:
            logger.error("Invalid json format: {}".format(e))
            raise ValueError("Invalid json format: {}".format(e))

        self.cold_start = dict_header.get('cold_start', True)

        self.task_ID = dict_header.get('task_ID', '')
        self.period_index = dict_header.get('period_index', '')

        # General info
        #: current version of GEDAP
        self.gedap_version = dict_header.get('gedap_version', gedap_version)
        #: login string to communicate with the RabbitMQ server
        #    (see code of module :mod:`scheduler.celery` and Celery documentation)
        self.rabbitmq_connect = dict_header.get('rabbitmq_connect', "amqp://gedap:gedap@gedap/rabbitmq_vhost")
        #: name of the satellite
        self.satellite = dict_header.get('mission', 'MFG')
        #: name of the radiometer
        self.radiometer = dict_header.get('radiometer', 'MVIRI')
        #: number of operational mission
        self.platform = dict_header['platform']
        #: subsatellite point: at the moment the default is the zeroth degree mission
        self.ssp = dict_header['SSP']

        if self.ssp not in [0,63]:
            logger.error("Sub-satellite point not recognized. Allowed values are 0 (0DM) and 63 (IODC).")
        
        # the notation of the subsatellite point in the cdr name generation is different than the one adopted for our tiles
        if self.ssp == 0:
            self.ssp_cdr = '00.0'
        if self.ssp == 63:
            self.ssp_cdr = '63.0'

        #: name of the projection
        self.projection = dict_header.get('projection', 'HIGH_RES')

        if self.radiometer == 'SEVIRI':
            self.v_st_bands = dict_header.get('bands', ['VIS06', 'VIS08', 'VIS16'])
        if self.radiometer == 'MVIRI':
            self.v_st_bands = dict_header.get('bands', ['VIS'])

        #else:
        #    logger.error("{} bands do not correspond to {} radiometer".format(self.v_st_bands, self.radiometer))

        # Info used to calculate the list of accumulation periods
        st_start_date = dict_header['start_date']
        st_end_date = dict_header['end_date']

        #: starting date of the whole period to process
        self.start_date = datetime.strptime(st_start_date, "%d-%m-%Y")
        #: ending date of the whole period to process
        self.end_date = datetime.strptime(st_end_date, "%d-%m-%Y")

        #: shift between two accumulation periods (in days)
        self.acc_period_shift = dict_header.get('acc_period_shift', 1)
        #: length of one accumulation periods (in days)
        self.acc_period_length = dict_header.get('acc_period_length', 1)

        #: maximum number of iterations during the descent (see CISAR documentation)
        self.max_nb_iterations = dict_header.get('max_nb_iterations', 20)

        # CISAR LUT
        # the monochromaticity correction requires the lut to vary with time.
        # therefore a check on the date has to be performed to select the closest lut available from the list
        #: Path to the directory containing the LUT (used by CISAR algorithm)
        self.dir_cisar_lut = dict_header['dir_cisar_lut']
        #: ID of the LUT to load (CISAR)
        self.lut_ID_prefix = dict_header.get('cisar_lut_ID', 'FAST_LUT')
        
        
        lut_prefix = os.path.join(self.lut_ID_prefix + '_' + self.platform + '_' + self.radiometer +\
                                         '_' + self.v_st_bands[0] + '_')
        all_luts = []

        if str(self.start_date.month).zfill(2) == '01':
            year_list = [self.start_date.year, self.start_date.year-1]
            for year in year_list:
                lut_list = glob(os.path.join(self.dir_cisar_lut, lut_prefix + str(year)) + '*' )
                all_luts.extend(lut_list)
        elif self.start_date.month == '12':
            year_list = [self.start_date.year, self.start_date.year +1]
            for year in year_list:
                lut_list = glob(os.path.join(self.dir_cisar_lut, lut_prefix + year + '*'))
                all_luts.extend(lut_list)
        else:
            lut_list = glob(os.path.join(self.dir_cisar_lut, lut_prefix + str(self.start_date.year)) + '*')
            all_luts = lut_list

        available_dates = []
        indices = []

        for index,lut in enumerate(all_luts):
            date = lut.split('_')[-1].split('.')[0]
            year = date[0:4]
            starting_date = datetime(int(year), 1, 1)
            dayoftheyear = date[4:]
            actual_date =  starting_date + timedelta(days=int(dayoftheyear))
            available_dates.append(actual_date)
            indices.append([actual_date, date])

        if len(all_luts) == 0:
            logger.error('There is no LUT available at {}. Check the path correctness.'.format(self.dir_cisar_lut))
            raise IOError('There is no LUT available at {}. Check the path correctness.'.format(self.dir_cisar_lut))

        def nearest(items, pivot):
            return min(items, key=lambda x: abs(x - pivot))

        closest_date = nearest(available_dates, self.start_date)

        for i in indices:
            if i[0] ==  closest_date:
                lut_date = i[1]

        self.cisar_lut_ID = os.path.join(self.lut_ID_prefix + '_' + self.platform + '_' + self.radiometer +\
                                         '_' + self.v_st_bands[0] + '_' + lut_date)


        for root, subdirs, files in os.walk(self.dir_cisar_lut):
            filenames = [os.path.splitext(i)[0] for i in files]
            if  self.cisar_lut_ID not in filenames:
                logger.error('The LUT {} is not located in the directory {}'.format(self.cisar_lut_ID, self.dir_cisar_lut))
                raise IOError(1110, 'The LUT {} is not located in the directory {}'.format(self.cisar_lut_ID, self.dir_cisar_lut))

        #: filename of the aeromod file (CISAR)
        self.cisar_aeromod_filename = dict_header.get('cisar_aeromod_filename', '')
        #: list of aerosol classes to include to process the inversion + their associated mode (FINE or COARSE)
        self.dict_aerosol_classes = dict_header.get('list_aerosol_classes', ["SIXSV_continental"])
        
        
        # Booleans to activate or deactivate the temporal smoothness and spectral
        # constraint terms of the cost function

        #: flag to (de)activate the temporal smoothness (crf CISAR)
        self.b_temporal_smoothness = dict_header.get('flag_temporal_smoothness', False)
        #: flag to (de)activate the spectral constraint terms (crf CISAR)
        self.b_spectral_constraint = dict_header.get('flag_spectral_constraint', False)
        #: flag to (de)activate the surface spectral constraint terms (crf CISAR)
        self.b_surface_spectral_constraint = dict_header.get('flag_surface_spectral_constraint', False)
        
        #: flag to (de)activate the aerosol height as being constant
        self.aerosol_height_constant = dict_header.get('aerosol_height_constant', False)

        #: flag to (de)activate the use of the CMASKSCORE field for the MVIRI cloud mask
        #: flag to (de)activate the use of the PCFC field for the MVIRI cloud mask

        self.useCMASKSCORE = dict_header.get('useCMASKSCORE', False)
        self.usePCFC = dict_header.get('usePCFC', False)

        #the pixel will be a cloud if the probability that not being a cloud is smaller than 20%
        #this allows a more conservative approach, increasing the number of cloudy pixels
        self.reliability_cfc_retrieval = 20
        # Stride for lines and columns (if stride_line = 2, we process only one line over two)
        #: line stride (to process the pixels located at each x lines and columns)
        self.line_stride = dict_header.get('line_stride',1)
        #: column stride (to process the pixels located at each x lines and columns)
        self.column_stride = dict_header.get('column_stride',1)
        
        # load values required to build prior and first guess

        dict_prior_fg = dict_header.get('prior_first_guess', {})

        #: list of prior values associated with surface parameters (see C++ documentation/code of PriorManager class)
        self.dict_prior_rpv = dict_prior_fg.get('prior_rpv', {'rho_0': [0.1], 'k': [0.5],
                                                              'theta': [-0.1], 'rho_c': [0.25]})
        #: sigma prior value associated with surface parameters (see C++ documentation/code of PriorManager class)
        self.sigma_prior_rpv = dict_prior_fg.get('sigma_prior_rpv', 0.5)
        #: list of prior values associated with AOT (see C++ documentation/code of PriorManager class)
        self.prior_aot = dict_prior_fg.get('prior_aot', [0.2])
        #: list of sigma prior values associated with AOT (see C++ documentation/code of PriorManager class)
        self.sigma_prior_aot = dict_prior_fg.get('sigma_prior_aot', 10.0)
        #: relative weight between fine and coarse classes (see C++ documentation/code of PriorManager class)
        self.coeff_fine_coarse = dict_prior_fg.get('coeff_fine_coarse', 0.04)

        #: perturbation factor value used to generate the first guess values from prior (see C++ documentation/code of TileProcessor class)
        self.pert_first_guess = dict_prior_fg.get('pert_first_guess', 0.4)
        #: lower first guess value associated with AOT (see C++ documentation of TileProcessor class)
        self.min_AOT_first_guess = dict_prior_fg.get('min_AOT_first_guess', 1.0e-4)
        #: upper first guess value associated with AOT (see C++ documentation of TileProcessor class)
        self.max_AOT_first_guess = dict_prior_fg.get('max_AOT_first_guess', 0.15)

        #: prior value for AOT in presence of dust (see C++ documentation/code of PriorManager class)
        self.dust_prior = dict_prior_fg.get('dust_prior', 2.5)
        #: sigma prior value for AOT in presence of dust (see C++ documentation/code of PriorManager class)
        self.dust_sigma_prior = dict_prior_fg.get('dust_sigma_prior', 1.5)
        #: first guess value for AOT in presence of dust (see C++ documentation/code of PriorManager class)
        self.dust_first_guess = dict_prior_fg.get('dust_first_guess', 3.5)


        # Load tiles and names of processing machines

        dict_tiles = dict_header['list_of_tiles']

        #: list of tiles with their properties (see :class:`TileProperties` class and common_init.json file).
        self.v_tiles = None

        #: list of processing machines associated with all active tiles (see common_init.json file).
        self.v_processing_machine_names = None

        #: load tiles dimensions + associated region and processing machine
        self.v_tiles, self.v_processing_machine_names = self.load_list_of_tiles(dict_tiles)


        #: name of the current machine
        self.machine_name = dict_header.get('machine_name', 'master')

        # Directories
        dict_dir_external_data = dict_header['dir_external_data']

        #check on existence of directories passed
        list_inputs = ["static", "image", "mask", "aot", "modparams"]
        for i in list_inputs:
            if not os.path.exists(dict_dir_external_data[i]):
                logger.error('The path to {} is not located in the directory {}'.format(i,dict_dir_external_data[i]))
                print(('The path to {} is not located in the directory {}'.format(i,dict_dir_external_data[i])))
                
        prefix = dict_header['filename_prefix']
        suffix = dict_header['filename_suffix']

        self.aotFolder = dict_dir_external_data['aot']
        if 'modparams' not in list(dict_dir_external_data.keys()):
            self.ecmwf_constant = True
        else:
            self.ecmwf_constant = False

        if self.ssp == 63:
            self.daily_cloudmask = True
        else:
            self.daily_cloudmask = False

        rootdir = dict_header.get('rootdir_tile_maker', os.getcwd())
        machine_name = self.machine_name
        self.rootdir_tile_maker = {machine_name: rootdir}

        # Check if directory exixts, otherwise create it, as well as create the subdirectories needed for GEDAP to run
        # i.e. Prior_Tiles, Input_Tiles, Solution_Tiles and also a directory to store the products.

        if not os.path.exists(self.rootdir_tile_maker[self.machine_name]):
            os.makedirs(self.rootdir_tile_maker[self.machine_name])
            os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Input_Tiles/'))
            os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Prior_Tiles/'))
            os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Solution_Tiles/'))
            os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Products/'))
        else:
            if not os.path.exists(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Input_Tiles/')) :
                os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Input_Tiles/'))
            if not os.path.exists(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Prior_Tiles/')):
                os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Prior_Tiles/'))
            if not os.path.exists(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Solution_Tiles/')):
                os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Solution_Tiles/'))
            if not os.path.exists(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Products/')):
                os.makedirs(os.path.join(self.rootdir_tile_maker[self.machine_name], 'Products/'))

        self.rootdir_tile_processor = self.rootdir_tile_maker[self.machine_name]
        self.subdir_input_tiles = 'Input_Tiles/'
        self.subdir_prior_tiles = 'Prior_Tiles/'
        self.subdir_solution_tiles = 'Solution_Tiles/'
        self.dir_products = os.path.join(self.rootdir_tile_maker[self.machine_name], 'Products/')

        # endregion


        # region Create filename generator
        # TODO: replace this with a more elegant selection mechanism
        self.filenamesGenerator = None

        if self.radiometer == 'MVIRI':
            self.filenamesGenerator = MFGFilenamesGenerator(dict_dir_external_data, prefix, suffix)
        else:
            logger.fatal("Unrecognized radiometer {}".format(self.radiometer))
            raise ValueError
        # endregion

    def get_v_tile(self):
        # Load tiles and names of processing machines                                                                                                                                                      \
                                                                                                                                                                                                            
        try:
            dict_header = json.load(open(os.path.join(self.dir_resources, 'common_init.json'), 'r'))
        except ValueError as e:
            logger.error("Invalid json format: {}".format(e))
            raise ValueError("Invalid json format: {}".format(e))
        try:
          self.idx_tile

        except AttributeError:
            dict_tiles = dict_header['list_of_tiles']
        else:
            dict_tiles={'TILE_{ti}'.format(ti=str(self.idx_tile).zfill(4)): 'A', 'ID': 'A(active)/N(not active)'}

        self.v_tiles = None
        self.v_processing_machine_names = None
        self.v_tiles, self.v_processing_machine_names = self.load_list_of_tiles(dict_tiles)


    def get_index_first_acc_period_to_process(self, logger):
        """
        Run trough the list of accumulation periods and return the index of the first "NOT PROCESSED" accumulation period.
        """
        period_index = 0
        while self.v_acc_periods[period_index].b_processed:
            period_index += 1
            if period_index==len(self.v_acc_periods):
                print('Processing completed. All periods have been processed, exiting GEDAP.')
                logger.info(logger.info('Processing completed. All periods have been processed, exiting GEDAP.'))
                exit()
        return period_index
    
    def get_list_tiles_associated_with_machine(self, machine):
        """
        Return the list of all tiles associated with a given machine (which is responsible to process)

        :param machine: name of the machine
        :type machine: string
        """
        # print "In ctx - self.v_tiles = ", self.v_tiles
        return [tile for tile in self.v_tiles if tile.machine == machine]
        
    def build_list_accumulation_periods(self):
        """
        Build the list of accumulation periods given the shift between acc. periods (in days) and their length (in days).\n
        The ith accumulation period is given by:
         * start_acc_period = start_date + i * shift
         * end_acc_period = start_acc_period + length
        """
        
        print('Building list of accumulation periods...')

        start = self.start_date
        end = self.end_date
        period_length = timedelta(days=self.acc_period_length - 1)
        period_shift = timedelta(days=self.acc_period_shift)

        if (self.end_date + timedelta(days=1) - self.start_date) < timedelta(days=self.acc_period_length):
            print(('Not enough days ({} days) to form an acquisition period ({} days)... exiting GEDAP.'.format(
                (self.end_date + timedelta(days=1) - self.start_date).days, self.acc_period_length)))
            logger.fatal('Not enough days ({} days) to form an acquisition period ({} days)... exiting GEDAP.'.format(
                (self.end_date + timedelta(days=1) - self.start_date).days, self.acc_period_length))
            exit()

        else:
            start_acc_period = start
            end_acc_period = start + period_length
            nb_acc_periods = 1
            while end_acc_period <= end:
                if end_acc_period > end:
                    end_acc_period = end
                self.v_acc_periods.append(AccumulationPeriod(start_acc_period, end_acc_period))
                nb_acc_periods += 1
                start_acc_period += period_shift
                end_acc_period += period_shift

            self.dump_list_accumulation_periods_in_json_file()

    def dump_list_accumulation_periods_in_json_file(self):
        """
        Update the list of accumulation periods in a JSON file with filename list_accumulation_periods.json
        located in resources folder.
        """
        
        index = 0
        dict_acc_periods = dict()
        for acc_period in self.v_acc_periods:
            st_status = 'PROCESSED' if acc_period.b_processed else 'NOT PROCESSED'
            dict_acc_periods[index] = acc_period.start_date.strftime('%Y-%m-%d') + '  ->  ' + acc_period.end_date.strftime('%Y-%m-%d') + ', ' + st_status
            index += 1
        
        with open(self.list_acc_periods_file, 'w') as file:
            json.dump(dict_acc_periods, file, indent=4, separators=(',', ':\t'))
        file.close()
        
    def load_list_accumulation_periods_from_json_file(self):
        """
        Load the list of accumulation periods fron the JSON file list_accumulation_periods.json located in
        resources folder.
        """

        # print 'Loading list of accumulation periods...'
        self.v_acc_periods = []
        dict_acc_periods = dict()

        with open(self.list_acc_periods_file, 'r') as file:
            try:
                dict_acc_periods = json.load(file, object_pairs_hook=OrderedDict)
            except ValueError as e:
                logger.fatal("Invalid json format in the list_accumulation_period json file, check at: {}".format(e))
                raise ValueError("Invalid json format in the list_accumulation_period json file, check at: {}".format(e))

        file.close()

        for st_acc_period in list(dict_acc_periods.values()):
            st_date = st_acc_period.split(', ')[0]
            st_start_date = st_date.split('  ->  ')[0]
            st_end_date = st_date.split('  ->  ')[1]

            st_status = st_acc_period.split(', ')[1]
            b_processed = st_status == 'PROCESSED'

            acc_period = AccumulationPeriod(datetime.strptime(st_start_date, '%Y-%m-%d'),
                                            datetime.strptime(st_end_date, '%Y-%m-%d'),
                                            b_processed)
            self.v_acc_periods.append(acc_period)


    def period_list_check(self):
        index = 0
        if self.v_acc_periods[index].start_date !=  self.start_date or ((self.v_acc_periods[-1].end_date +timedelta(days=1)) - self.end_date) > timedelta(days=self.acc_period_length):


            print('The period to be processed has changed. Cold start required or realign periods with list_accumulation_period.json.')
            logger.fatal('The period to be processed has changed. Cold start required or realign periods with list_accumulation_period.json.')
            exit()



    def load_list_of_tiles(self, dict_tiles):
        """
        Build the list of tiles.

        :param dict_tiles: dictionary (ID, tile properties as string) loaded from the file common_init.json located in resources folder.
        :type dict_tiles: dictionary
        """
        
        # print 'Loading list of tiles...'
        
        v_tiles = []
        v_processing_machine_names = []

        # Remove the first comment line in common_init
        if 'processing_machine' in dict_tiles['ID']:
            OverwriteMachineName = True
        else:
            OverwriteMachineName = False


        del dict_tiles['ID']

        for st_id, st_tile_properties in list(dict_tiles.items()):
            v_tile_properties = st_tile_properties.split(', ')

            if OverwriteMachineName:
                st_machine = v_tile_properties[0]
                is_activated = v_tile_properties[1] == 'A'
            else:
                st_machine = 'master'
                is_activated = v_tile_properties[0] == 'A'

            if is_activated:

                file_path = os.path.join(self.dir_resources, str.lower(str(self.radiometer))+'_tiles_'+ str(self.ssp).zfill(4) +'.txt')

                tile_dictionary =  read_col(file_path, "\t", i_head=0)

                if st_id not in tile_dictionary['TileID']:
                    logger.error("The TileID {} is not a valid one. Skipping this tile".format(st_id))
                    continue
                else:
                    index = tile_dictionary['TileID'].index(st_id)

                    v_prior_aot       = self.prior_aot
                    sigma_prior_aot   = self.sigma_prior_aot
                    coeff_fine_coarse = self.coeff_fine_coarse

                    tile_properties = TileProperties(st_id, st_machine, tile_dictionary['first_line'][index],
                                                     tile_dictionary['last_line'][index], tile_dictionary['first_column'][index],
                                                     tile_dictionary['last_column'][index],
                                                     self.dict_prior_rpv, self.sigma_prior_rpv, v_prior_aot, sigma_prior_aot, coeff_fine_coarse,
                                                     self.pert_first_guess, self.min_AOT_first_guess, self.max_AOT_first_guess)

                    v_tiles.append(tile_properties)
                    if st_machine not in v_processing_machine_names:
                        v_processing_machine_names.append(st_machine)

        return v_tiles, v_processing_machine_names


class List_Aerosol_Classes(object):
    """
    This class is responsible to extract aerosol class properties from aeromod files (CISAR LUT) and to store them.
    """
    
    def __init__(self):    
        self.list_aer_class_all = []   

        self.dic_aer_class = dict()   
        
        self.nb_aer_class = len(self.list_aer_class_all)
        
    def extract_list_and_properties_of_aerosol_classes(self, dir_LUT, aeromod_file):
        """
        Extract the list of aerosol classes and their properties from an aeromod file (CISAR LUT)

        :param dir_LUT: path to the directory containing the LUT and aeromod files
        :type dir_LUT: string
        :param aeromod_file: name of the aeromod file to be read
        :type aeromod_file: string
        """
        
        # Read names and properties of aerosol classes
        # Create a dictionnary
        # print "aeromod file ", aeromod_file, " exists? : ", os.path.isfile(dir_LUT + aeromod_file)
        with open(dir_LUT + aeromod_file, 'r') as f_aermod:

            # Get number of aerosol classes and bands
            self.nb_aer_class = int(f_aermod.readline().split()[0])
            nb_bands = int(f_aermod.readline().split()[0])
            f_aermod.readline()

            # Reset list and dict
            self.list_aer_class_all = []
            self.dic_aer_class = dict()

            # Extract aerosol class properties
            for i_aer in range(self.nb_aer_class):
                # Read aerosol class name
                mode, aer_class_name = f_aermod.readline().split()
                self.list_aer_class_all.append(aer_class_name)

                # Skip da_ext_coef @ 0.55
                f_aermod.readline()

                # Read w0 and g for each band
                v_w0 = []
                v_g = []
                for i_band in range(nb_bands):
                    scattering_coef = float(f_aermod.readline().split()[0])
                    extinction_coef = float(f_aermod.readline().split()[0])
                    g = float(f_aermod.readline().split()[1])
                    w0 = scattering_coef / extinction_coef
                    v_w0.append(w0)
                    v_g.append(g)

                # Create one dict for each aerosol class and add it to dict grouping all classes
                # (it's a dict of dict)
                dic_one_class = {'w0': v_w0, 'g': v_g, 'mode': mode}
                self.dic_aer_class[aer_class_name] = dic_one_class

            # print('aerosol classes:')
            # for aer_class_name, prop_aer_class in self.dic_aer_class.iteritems():
            #     print(aer_class_name, ': mode -> ', prop_aer_class['mode'],
            #           ', w0 -> ', prop_aer_class['w0'],
            #           ', g -> ', prop_aer_class['g'])
            # print ''

    def get_aer_class_index(self, aer_class_name):
        """
        Return the index of an aerosol class given its name.

        :param aer_class_name: name of the aerosol class
        :type aer_class_name: string
        """
        return self.list_aer_class_all.index(aer_class_name)


# create a default instance (~ singleton)
# TODO: check if this is a good thing to do
aerClasses = List_Aerosol_Classes()
