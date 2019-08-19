import numpy as np
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu

SWATH_WIDTH = 28
NUM_CHANNELS = 5
NUM_THERMISTORS = 18
ANCIL_VAL = 10
CALIB_NUMBER = 4


class SSMT2:
    @staticmethod
    def add_original_variables(dataset, height, srf_size=None):
        tu.add_geolocation_variables(dataset, SWATH_WIDTH, height)
        tu.add_quality_flags(dataset, SWATH_WIDTH, height)

        # Temperature_misc_housekeeping
        default_array = DefaultData.create_default_array(height, NUM_THERMISTORS, np.float32, dims_names=["housekeeping", "y"], fill_value=np.NaN)
        variable = Variable(["housekeeping", "y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        variable.attrs["units"] = "TODO"
        dataset["Temperature_misc_housekeeping"] = variable

        # ancil_data
        default_array = DefaultData.create_default_array(height, ANCIL_VAL, np.float64, dims_names=["ancil_val", "y"], fill_value=np.NaN)
        variable = Variable(["ancil_val", "y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Additional per scan information: year, day_of_year, secs_of_day, sat_lat, " \
                                      "sat_long, sat_alt, sat_heading, year, day_of_year, secs_of_day"
        dataset["ancil_data"] = variable

        # channel_quality_flag
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        dataset["channel_quality_flag"] = variable

        # cold_counts
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, CALIB_NUMBER, np.float32, np.NaN)
        variable = Variable(["calib_number", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["cold_counts"] = variable

        # counts_to_tb_gain
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.float32, dims_names=["channel", "y"], fill_value=np.NaN)
        variable = Variable(["channel", "y", ], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["counts_to_tb_gain"] = variable

        # counts_to_tb_offset
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.float32, dims_names=["channel", "y"], fill_value=np.NaN)
        variable = Variable(["channel", "y", ], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["counts_to_tb_offset"] = variable

        # gain_control
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.float32, dims_names=["channel", "y"], fill_value=np.NaN)
        variable = Variable(["channel", "y", ], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["gain_control"] = variable

        # tb
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        variable.attrs["standard_name"] = "toa_brightness_temperature"
        tu.add_units(variable, "K")
        dataset["tb"] = variable

        # thermal_reference
        default_array = DefaultData.create_default_vector(height, np.float32, np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        tu.add_units(variable, "TODO")
        dataset["thermal_reference"] = variable

        # warm_counts
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, CALIB_NUMBER, np.float32, np.NaN)
        variable = Variable(["calib_number", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["warm_counts"] = variable

    @staticmethod
    def get_swath_width():
        return SWATH_WIDTH

    @staticmethod
    def add_easy_fcdr_variables(dataset, height, corr_dx=None, corr_dy=None, lut_size=None):
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        tu.add_units(variable, "K")
        variable.attrs["long_name"] = "independent uncertainty per pixel"
        dataset["u_independent_tb"] = variable

        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        tu.add_units(variable, "K")
        variable.attrs["long_name"] = "structured uncertainty per pixel"
        dataset["u_structured_tb"] = variable

    @staticmethod
    def add_full_fcdr_variables(dataset, height):
        # u_Temperature_misc_housekeeping
        default_array = DefaultData.create_default_array(height, NUM_THERMISTORS, np.float32, dims_names=["housekeeping", "y"], fill_value=np.NaN)
        variable = Variable(["housekeeping", "y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        variable.attrs["units"] = "TODO"
        dataset["u_Temperature_misc_housekeeping"] = variable

        # u_cold_counts
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, CALIB_NUMBER, np.float32, np.NaN)
        variable = Variable(["calib_number", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["u_cold_counts"] = variable

        # u_counts_to_tb_gain
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.float32, dims_names=["channel", "y"], fill_value=np.NaN)
        variable = Variable(["channel", "y", ], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["u_counts_to_tb_gain"] = variable

        # u_counts_to_tb_offset
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.float32, dims_names=["channel", "y"], fill_value=np.NaN)
        variable = Variable(["channel", "y", ], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["u_counts_to_tb_offset"] = variable

        # u_gain_control
        default_array = DefaultData.create_default_array(height, NUM_CHANNELS, np.float32, dims_names=["channel", "y"], fill_value=np.NaN)
        variable = Variable(["channel", "y", ], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["u_gain_control"] = variable

        # u_tb
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, NUM_CHANNELS, np.float32, np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        tu.add_units(variable, "K")
        dataset["u_tb"] = variable

        # u_thermal_reference
        default_array = DefaultData.create_default_vector(height, np.float32, np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        tu.add_units(variable, "TODO")
        dataset["u_thermal_reference"] = variable

        # u_warm_counts
        default_array = DefaultData.create_default_array_3d(SWATH_WIDTH, height, CALIB_NUMBER, np.float32, np.NaN)
        variable = Variable(["calib_number", "y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "TODO"
        dataset["u_warm_counts"] = variable

    @staticmethod
    def add_template_key(dataset):
        dataset.attrs["template_key"] = "SSMT2"
