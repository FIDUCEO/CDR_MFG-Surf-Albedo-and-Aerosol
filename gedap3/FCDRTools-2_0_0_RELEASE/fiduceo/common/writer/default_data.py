import numpy as np
import xarray as xr


class DefaultData:
    @staticmethod
    def create_default_vector(size, dtype, fill_value=None):
        if fill_value is None:
            fill_value = DefaultData.get_default_fill_value(dtype)

        empty_array = np.full([size], fill_value, dtype)

        default_array = xr.DataArray(empty_array, dims=['y'])
        return default_array

    @staticmethod
    def create_default_array(width, height, dtype, dims_names=None, fill_value=None):
        if fill_value is None:
            fill_value = DefaultData.get_default_fill_value(dtype)

        empty_array = np.full([height, width], fill_value, dtype)

        if dims_names is not None:
            default_array = xr.DataArray(empty_array, dims=dims_names)
        else:
            default_array = xr.DataArray(empty_array, dims=['y', 'x'])

        return default_array

    @staticmethod
    def create_default_array_3d(width, height, num_channels, dtype, fill_value=None, dims_names=None):
        if fill_value is None:
            fill_value = DefaultData.get_default_fill_value(dtype)

        empty_array = np.full([num_channels, height, width], fill_value, dtype)

        if dims_names is not None:
            default_array = xr.DataArray(empty_array, dims=dims_names)
        else:
            default_array = xr.DataArray(empty_array, dims=['channel', 'y', 'x'])

        return default_array

    @staticmethod
    def create_default_array_4d(width, height, z1, z2, dtype, fill_value=None, dims_names=None):
        if fill_value is None:
            fill_value = DefaultData.get_default_fill_value(dtype)

        empty_array = np.full([z2, z1, height, width], fill_value, dtype)

        if dims_names is not None:
            default_array = xr.DataArray(empty_array, dims=dims_names)
        else:
            default_array = xr.DataArray(empty_array, dims=['ssp_y', 'ssp_x', 'y', 'x'])

        return default_array

    @staticmethod
    def get_default_fill_value(dtype):
        """
        Returns a CF conforming default fill value for the data type
        :param dtype: numpy dtype
        :return: the fill value
        :type: numpy dtype
        """
        if dtype == np.int8:
            return np.int8(-127)
        if dtype == np.uint8:
            return np.uint8(-1)
        elif dtype == np.int16:
            return np.int16(-32767)
        elif dtype == np.uint16:
            return np.uint16(-1)
        elif dtype == np.int32:
            return np.int32(-2147483647)
        elif dtype == np.uint32:
            return np.uint32(-1)
        elif dtype == np.int64:
            return np.int64(-9223372036854775806)
        elif dtype == np.float32:
            return np.float32(9.96921E36)
        elif dtype == np.float64:
            return np.float64(9.969209968386869E36)
