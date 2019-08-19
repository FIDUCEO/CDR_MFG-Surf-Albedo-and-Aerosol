import numpy as np
from xarray import Variable, Coordinate

from fiduceo.common.writer.default_data import DefaultData

LATITUDE_UNIT = "degrees_north"
LONGITUDE_UNIT = "degrees_east"

LAT_NAME = "latitude"
LON_NAME = "longitude"


class TemplateUtil:

    @staticmethod
    def add_geolocation_variables(dataset, width, height, chunksizes=None):
        default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=np.NaN)

        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = LAT_NAME
        TemplateUtil.add_units(variable, LATITUDE_UNIT)
        TemplateUtil.add_encoding(variable, np.int16, -32768, scale_factor=0.0027466658, chunksizes=chunksizes)
        dataset[LAT_NAME] = variable

        default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = LON_NAME
        TemplateUtil.add_units(variable, LONGITUDE_UNIT)
        TemplateUtil.add_encoding(variable, np.int16, -32768, scale_factor=0.0054933317, chunksizes=chunksizes)
        dataset[LON_NAME] = variable

    @staticmethod
    def add_gridded_geolocation_variables(dataset, width, height):
        default_array = DefaultData.create_default_vector(height, np.float32, fill_value=np.NaN)
        variable = Variable(["y"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["standard_name"] = LAT_NAME
        variable.attrs["long_name"] = LAT_NAME
        variable.attrs["bounds"] = "lat_bnds"
        TemplateUtil.add_units(variable, LATITUDE_UNIT)
        dataset["lat"] = variable

        default_array = DefaultData.create_default_array(2, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "bounds"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "latitude cell boundaries"
        TemplateUtil.add_units(variable, LATITUDE_UNIT)
        dataset["lat_bnds"] = variable

        default_array = DefaultData.create_default_vector(width, np.float32, fill_value=np.NaN)
        variable = Variable(["x"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["standard_name"] = LON_NAME
        variable.attrs["long_name"] = LON_NAME
        TemplateUtil.add_units(variable, LONGITUDE_UNIT)
        variable.attrs["bounds"] = "lon_bnds"
        dataset["lon"] = variable

        default_array = DefaultData.create_default_array(2, width, np.float32, fill_value=np.NaN)
        variable = Variable(["x", "bounds"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        TemplateUtil.add_units(variable, LONGITUDE_UNIT)
        variable.attrs["long_name"] = "longitude cell boundaries"
        dataset["lon_bnds"] = variable

    @staticmethod
    def add_quality_flags(dataset, width, height, chunksizes=None, masks_append=None, meanings_append=None):
        default_array = DefaultData.create_default_array(width, height, np.uint8, fill_value=0)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "status_flag"
        TemplateUtil.add_geolocation_attribute(variable)

        masks = "1, 2, 4, 8, 16, 32, 64, 128"
        if masks_append is not None:
            masks = masks + masks_append
        variable.attrs["flag_masks"] = masks

        meanings = "invalid use_with_caution invalid_input invalid_geoloc invalid_time sensor_error padded_data incomplete_channel_data"
        if meanings_append is not None:
            meanings = meanings + meanings_append
        variable.attrs["flag_meanings"] = meanings

        if chunksizes is not None:
            TemplateUtil.add_chunking(variable, chunksizes)
        dataset["quality_pixel_bitmask"] = variable

    @staticmethod
    def create_scalar_float_variable(long_name=None, standard_name=None, units=None, fill_value=np.NaN):
        default_array = fill_value

        variable = Variable([], default_array)
        TemplateUtil.add_fill_value(variable, fill_value)

        if long_name is not None:
            variable.attrs["long_name"] = long_name

        if standard_name is not None:
            variable.attrs["standard_name"] = standard_name

        if units is not None:
            TemplateUtil.add_units(variable, units)
        return variable

    @staticmethod
    def create_float_variable(width, height, standard_name=None, long_name=None, dim_names=None, fill_value=None):
        if fill_value is None:
            default_array = DefaultData.create_default_array(width, height, np.float32)
        else:
            default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=fill_value)

        if dim_names is None:
            variable = Variable(["y", "x"], default_array)
        else:
            variable = Variable(dim_names, default_array)

        if fill_value is None:
            variable.attrs["_FillValue"] = DefaultData.get_default_fill_value(np.float32)
        else:
            variable.attrs["_FillValue"] = fill_value

        if standard_name is not None:
            variable.attrs["standard_name"] = standard_name

        if long_name is not None:
            variable.attrs["long_name"] = long_name

        return variable

    @staticmethod
    def set_unsigned(variable):
        variable.attrs["_Unsigned"] = "true"

    @staticmethod
    def add_fill_value(variable, fill_value):
        variable.attrs["_FillValue"] = fill_value

    @staticmethod
    def add_units(variable, units):
        variable.attrs["units"] = units

    @staticmethod
    def add_scale_factor(variable, scale_factor):
        variable.attrs["scale_factor"] = scale_factor

    @staticmethod
    def add_offset(variable, offset):
        variable.attrs["add_offset"] = offset

    @staticmethod
    def add_geolocation_attribute(variable):
        variable.attrs["coordinates"] = LON_NAME + " " + LAT_NAME

    @staticmethod
    def add_encoding(variable, data_type, fill_value, scale_factor=1.0, offset=0.0, chunksizes=None):
        if chunksizes is None:
            variable.encoding = dict([('dtype', data_type), ('_FillValue', fill_value), ('scale_factor', scale_factor), ('add_offset', offset)])
        else:
            variable.encoding = dict([('dtype', data_type), ('_FillValue', fill_value), ('scale_factor', scale_factor), ('add_offset', offset), ('chunksizes', chunksizes)])

    @staticmethod
    def add_chunking(variable, chunksizes):
        variable.encoding = dict([('chunksizes', chunksizes)])

    @staticmethod
    def add_coordinates(ds, channel_data=None):
        x_dim = ds.dims["x"]
        ds["x"] = Coordinate("x", np.arange(x_dim, dtype=np.uint16))

        y_dim = ds.dims["y"]
        ds["y"] = Coordinate("y", np.arange(y_dim, dtype=np.uint16))

        if "channel" in ds.dims and channel_data is None:
            channels_dim = ds.dims["channel"]
            ds["channel"] = Coordinate("channel", np.arange(channels_dim, dtype=np.uint16))

        if channel_data is not None:
            ds["channel"] = Coordinate("channel", channel_data)

    @staticmethod
    def add_correlation_matrices(dataset, num_channels):
        default_array = np.diag(np.ones(num_channels, dtype=np.float32))
        variable = Variable(["channel", "channel"], default_array)
        variable.attrs["long_name"] = "Channel_correlation_matrix_independent_effects"
        variable.attrs["units"] = "1"
        variable.attrs["valid_min"] = "-10000"
        variable.attrs["valid_max"] = "10000"
        variable.attrs["description"] = "Channel error correlation matrix for independent effects"
        TemplateUtil.add_encoding(variable, np.int16, -32768, scale_factor=0.0001)
        dataset['channel_correlation_matrix_independent'] = variable

        default_array = np.diag(np.ones(num_channels, dtype=np.float32))
        variable = Variable(["channel", "channel"], default_array)
        variable.attrs["long_name"] = "Channel_correlation_matrix_structured_effects"
        variable.attrs["units"] = "1"
        variable.attrs["valid_min"] = "-10000"
        variable.attrs["valid_max"] = "10000"
        variable.attrs["description"] = "Channel error correlation matrix for structured effects"
        TemplateUtil.add_encoding(variable, np.int16, -32768, scale_factor=0.0001)
        dataset['channel_correlation_matrix_structured'] = variable

    @staticmethod
    def add_lookup_tables(dataset, num_channels, lut_size):
        default_array = DefaultData.create_default_array(num_channels, lut_size, np.float32, fill_value=np.NaN)
        variable = Variable(["lut_size", "channel"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["description"] = "Lookup table to convert radiance to brightness temperatures"
        dataset['lookup_table_BT'] = variable

        default_array = DefaultData.create_default_array(num_channels, lut_size, np.float32, fill_value=np.NaN)
        variable = Variable(["lut_size", "channel"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["description"] = "Lookup table to convert brightness temperatures to radiance"
        dataset['lookup_table_radiance'] = variable

    @staticmethod
    def add_correlation_coefficients(dataset, num_channels, delta_x, delta_y):
        default_array = DefaultData.create_default_array(num_channels, delta_x, np.float32, fill_value=np.NaN)
        variable = Variable(["delta_x", "channel"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "cross_element_correlation_coefficients"
        variable.attrs["description"] = "Correlation coefficients per channel for scanline correlation"
        dataset['cross_element_correlation_coefficients'] = variable

        default_array = DefaultData.create_default_array(num_channels, delta_y, np.float32, fill_value=np.NaN)
        variable = Variable(["delta_y", "channel"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "cross_line_correlation_coefficients"
        variable.attrs["description"] = "Correlation coefficients per channel for inter scanline correlation"
        dataset['cross_line_correlation_coefficients'] = variable

    @staticmethod
    def create_CDR_uncertainty(width, height, description, coordinates=None, units=None):
        default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        TemplateUtil.add_fill_value(variable, np.NaN)
        variable.attrs["description"] = description

        if units is None:
            TemplateUtil.add_units(variable, "%")
        else:
            TemplateUtil.add_units(variable, units)

        if coordinates is None:
            variable.attrs["coordinates"] = "longitude latitude"
        else:
            variable.attrs["coordinates"] = coordinates

        return variable
