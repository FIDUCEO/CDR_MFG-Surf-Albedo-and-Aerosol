import numpy as np
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu
from fiduceo.fcdr.writer.templates.mviri import FULL_SIZE, IR_SIZE, IR_X_DIMENSION

SSP_X = 2
SSP_Y = 2


class MVIRI_STATIC:
    @staticmethod
    def add_original_variables(dataset, height, srf_size=None):
        # height is ignored - supplied just for interface compatibility tb 2017-07-19
        # latitude_vis
        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "latitude"
        tu.add_units(variable, "degrees_north")
        tu.add_encoding(variable, np.int16, -32768, scale_factor=0.0027466658)
        dataset["latitude_vis"] = variable

        # longitude_vis
        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "longitude"
        tu.add_units(variable, "degrees_east")
        tu.add_encoding(variable, np.int16, -32768, scale_factor=0.0054933317)
        dataset["longitude_vis"] = variable

        # latitude_ir_wv
        default_array = DefaultData.create_default_array(IR_SIZE, IR_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable([IR_X_DIMENSION, IR_X_DIMENSION], default_array)
        variable.attrs["standard_name"] = "latitude"
        tu.add_units(variable, "degrees_north")
        tu.add_encoding(variable, np.int16, -32768, scale_factor=0.0027466658)
        dataset["latitude_ir_wv"] = variable

        # longitude_ir_wv
        default_array = DefaultData.create_default_array(IR_SIZE, IR_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable([IR_X_DIMENSION, IR_X_DIMENSION], default_array)
        variable.attrs["standard_name"] = "longitude"
        tu.add_units(variable, "degrees_east")
        tu.add_encoding(variable, np.int16, -32768, scale_factor=0.0054933317)
        dataset["longitude_ir_wv"] = variable

    @staticmethod
    def add_template_key(dataset):
        dataset.attrs["template_key"] = "MVIRI_STATIC"

    @staticmethod
    def add_easy_fcdr_variables(dataset, height, corr_dx=None, corr_dy=None, lut_size=None):
        pass  # not required in this class

    @staticmethod
    def add_full_fcdr_variables(dataset, height):
        pass  # not required in this class

    @staticmethod
    def get_swath_width():
        return FULL_SIZE
