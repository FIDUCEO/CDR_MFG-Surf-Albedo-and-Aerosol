import numpy as np
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu

CHUNKING = (1280, 409)


class SST_ENSEMBLE:

    @staticmethod
    def add_variables(dataset, width, height, num_samples=10):
        tu.add_geolocation_variables(dataset, width, height, chunksizes=CHUNKING)
        tu.add_quality_flags(dataset, width, height, chunksizes=CHUNKING)

        default_array = DefaultData.create_default_vector(height, np.int32, fill_value=4294967295)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, 4294967295)
        variable.attrs["standard_name"] = "time"
        variable.attrs["long_name"] = "Acquisition time in seconds since 1970-01-01 00:00:00"
        tu.add_units(variable, "s")
        dataset["time"] = variable

        default_array = DefaultData.create_default_array_3d(width, height, num_samples, np.float32, fill_value=np.NaN)
        variable = Variable(["samples","y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["standard_name"] = "sea_surface_temperature"
        variable.attrs["units"] = "K"
        variable.attrs["coordinates"] = "longitude latitude"
        dataset["sst"] = variable
