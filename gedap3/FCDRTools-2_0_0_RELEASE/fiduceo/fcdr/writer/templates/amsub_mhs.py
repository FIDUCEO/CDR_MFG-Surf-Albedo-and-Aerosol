import numpy as np
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu

NUM_CHANNELS = 5
BTEMPS_FILL_VALUE = -999999
SWATH_WIDTH = 90


class AMSUB_MHS:
    @staticmethod
    def add_original_variables(dataset, height, srf_size=None, corr_dx=None, corr_dy=None, lut_size=None):
        tu.add_geolocation_variables(dataset, SWATH_WIDTH, height)
        tu.add_quality_flags(dataset, SWATH_WIDTH, height)

        # btemps
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        variable.attrs["standard_name"] = "toa_brightness_temperature"
        tu.add_encoding(variable, np.int32, -999999, scale_factor=0.01)
        tu.add_units(variable, "K")
        variable.attrs["ancillary_variables"] = "chanqual qualind scanqual"
        dataset["btemps"] = variable

        # chanqual
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.int32, dims_names=["channel", "y"], fill_value=0)
        variable = Variable(["channel", "y"], default_array)
        variable.attrs["standard_name"] = "status_flag"
        variable.attrs["flag_masks"] = "1, 2, 4, 8, 16, 32"
        variable.attrs["flag_meanings"] = "some_bad_prt_temps some_bad_space_view_counts some_bad_bb_counts no_good_prt_temps no_good_space_view_counts no_good_bb_counts"
        dataset["chanqual"] = variable

        # instrtemp
        default_array = DefaultData.create_default_vector(height, np.float32, fill_value=np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_units(variable, "K")
        tu.add_encoding(variable, np.int32, DefaultData.get_default_fill_value(np.int32), scale_factor=0.01)
        variable.attrs["long_name"] = "instrument_temperature"
        dataset["instrtemp"] = variable

        # qualind
        default_array = DefaultData.create_default_vector(height, np.int32, fill_value=0)
        variable = Variable(["y"], default_array)
        variable.attrs["standard_name"] = "status_flag"
        variable.attrs["flag_masks"] = "33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648"
        variable.attrs["flag_meanings"] = "instr_status_changed first_good_clock_update no_earth_loc no_calib data_gap_precedes time_seq_error not_use_scan"
        dataset["qualind"] = variable

        # scanqual
        default_array = DefaultData.create_default_vector(height, np.int32, fill_value=0)
        variable = Variable(["y"], default_array)
        variable.attrs["standard_name"] = "status_flag"
        variable.attrs["flag_masks"] = "8, 16, 32, 64, 128, 1024, 2048, 4096, 8192, 16384, 32768, 1048576, 2097152, 4194304, 8388608"
        variable.attrs[
            "flag_meanings"] = "earth_loc_quest_ant_pos earth_loc_quest_reas earth_loc_quest_margin earth_loc_quest_time no_earth_loc_time uncalib_instr_mode uncalib_channels calib_marg_prt uncalib_bad_prt calib_few_scans uncalib_bad_time repeat_scan_times inconsistent_time time_field_bad time_field_inferred"
        dataset["scanqual"] = variable

        # scnlin
        default_array = DefaultData.create_default_vector(height, np.int32)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int32))
        variable.attrs["long_name"] = "scanline"
        dataset["scnlin"] = variable

        # scnlindy
        default_array = DefaultData.create_default_vector(height, np.int32)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int32))
        variable.attrs["long_name"] = "Acquisition day of year of scan"
        dataset["scnlindy"] = variable

        # scnlintime
        default_array = DefaultData.create_default_vector(height, np.int32)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int32))
        variable.attrs["long_name"] = "Acquisition time of scan in milliseconds since beginning of the day"
        tu.add_units(variable, "ms")
        dataset["scnlintime"] = variable

        # scnlinyr
        default_array = DefaultData.create_default_vector(height, np.int32)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int32))
        variable.attrs["long_name"] = "Acquisition year of scan"
        dataset["scnlinyr"] = variable

        # satellite_azimuth_angle
        variable = AMSUB_MHS.create_angle_variable(height, "sensor_azimuth_angle")
        dataset["satellite_azimuth_angle"] = variable

        # satellite_zenith_angle
        variable = AMSUB_MHS.create_angle_variable(height, "sensor_zenith_angle")
        dataset["satellite_zenith_angle"] = variable

        # solar_azimuth_angle
        variable = AMSUB_MHS.create_angle_variable(height, "solar_azimuth_angle")
        dataset["solar_azimuth_angle"] = variable

        # solar_zenith_angle
        variable = AMSUB_MHS.create_angle_variable(height, "solar_zenith_angle")
        dataset["solar_zenith_angle"] = variable

        # acquisition_time
        default_array = DefaultData.create_default_vector(height, np.int32)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int32))
        variable.attrs["standard_name"] = "time"
        variable.attrs["long_name"] = "Acquisition time in seconds since 1970-01-01 00:00:00"
        tu.add_units(variable, "s")
        dataset["acquisition_time"] = variable

    @staticmethod
    def get_swath_width():
        return SWATH_WIDTH

    @staticmethod
    def add_easy_fcdr_variables(dataset, height, srf_size=None, corr_dx=None, corr_dy=None, lut_size=None):
        # u_independent_btemps
        variable = AMSUB_MHS._create_3d_float_variable(height)
        tu.add_units(variable, "K")
        variable.attrs["long_name"] = "independent uncertainty per pixel"
        dataset["u_independent_btemps"] = variable

        # u_structured_btemps
        variable = AMSUB_MHS._create_3d_float_variable(height)
        tu.add_units(variable, "K")
        variable.attrs["long_name"] = "structured uncertainty per pixel"
        dataset["u_structured_btemps"] = variable

    @staticmethod
    def add_full_fcdr_variables(dataset, height):
        # u_btemps
        variable = AMSUB_MHS._create_3d_float_variable(height)
        variable.attrs["long_name"] = "total uncertainty of brightness temperature"
        tu.add_units(variable, "K")
        dataset["u_btemps"] = variable

        # u_syst_btemps
        variable = AMSUB_MHS._create_3d_float_variable(height)
        variable.attrs["long_name"] = "systematic uncertainty of brightness temperature"
        tu.add_units(variable, "K")
        dataset["u_syst_btemps"] = variable

        # u_random_btemps
        variable = AMSUB_MHS._create_3d_float_variable(height)
        variable.attrs["long_name"] = "noise on brightness temperature"
        tu.add_units(variable, "K")
        dataset["u_random_btemps"] = variable

        # u_instrtemp
        default_array = DefaultData.create_default_vector(height, np.float32, fill_value=np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "uncertainty of instrument temperature"
        tu.add_units(variable, "K")
        dataset["u_instrtemp"] = variable

        # u_latitude
        variable = AMSUB_MHS.create_angle_uncertainty_variable("latitude", height)
        dataset["u_latitude"] = variable

        # u_longitude
        variable = AMSUB_MHS.create_angle_uncertainty_variable("longitude", height)
        dataset["u_longitude"] = variable

        # u_satellite_azimuth_angle
        variable = AMSUB_MHS.create_angle_uncertainty_variable("satellite azimuth angle", height)
        dataset["u_satellite_azimuth_angle"] = variable

        # u_satellite_zenith_angle
        variable = AMSUB_MHS.create_angle_uncertainty_variable("satellite zenith angle", height)
        dataset["u_satellite_zenith_angle"] = variable

        # u_solar_azimuth_angle
        variable = AMSUB_MHS.create_angle_uncertainty_variable("solar azimuth angle", height)
        dataset["u_solar_azimuth_angle"] = variable

        # u_solar_zenith_angle
        variable = AMSUB_MHS.create_angle_uncertainty_variable("solar zenith angle", height)
        dataset["u_solar_zenith_angle"] = variable

    @staticmethod
    def add_template_key(dataset):
        dataset.attrs["template_key"] = "AMSUB_MHS"

    @staticmethod
    def create_angle_uncertainty_variable(angle_name, height):
        variable = tu.create_float_variable(SWATH_WIDTH, height, long_name="uncertainty of " + angle_name, fill_value=np.NaN)
        tu.add_units(variable, "degree")
        return variable

    @staticmethod
    def create_angle_variable(height, standard_name):
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = standard_name
        tu.add_units(variable, "degree")
        tu.add_encoding(variable, np.int32, -999999, scale_factor=0.01)
        return variable

    @staticmethod
    def _create_3d_float_variable(height):
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        return variable
