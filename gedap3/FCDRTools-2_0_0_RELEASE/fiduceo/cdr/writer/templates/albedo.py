import numpy as np
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu

CHUNKING = (500, 500)


class Albedo:

    @staticmethod
    def add_variables(dataset, width, height):
        # @todo 1 tb/tb add geolocation 2018-06-25

        tu.add_quality_flags(dataset, width, height, chunksizes=CHUNKING)

        default_array = DefaultData.create_default_vector(height, np.int32, fill_value=-1)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, -1)
        variable.attrs["standard_name"] = "time"
        variable.attrs["long_name"] = "Acquisition time in seconds since 1970-01-01 00:00:00"
        tu.add_units(variable, "s")
        dataset["time"] = variable

        default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["standard_name"] = "surface_albedo"
        variable.attrs["coordinates"] = "longitude latitude"
        tu.add_chunking(variable, CHUNKING)
        dataset["surface_albedo"] = variable

        dataset["u_independent_surface_albedo"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of surface_albedo due to independent effects")
        dataset["u_structured_surface_albedo"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of surface_albedo due to structured effects")
        dataset["u_common_surface_albedo"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of surface_albedo due to common effects")
