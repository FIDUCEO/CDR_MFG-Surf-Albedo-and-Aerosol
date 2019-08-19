import re
from math import pi

import numexpr as ne
import numpy as np
import xarray as xr
from gridtools.resampling import resample_2d


class FCDRReader:

    @classmethod
    def read(cls, file_str, drop_variables_str=None, decode_cf=True, decode_times=True, engine_str=None):
        """Read a dataset from a netCDF 3/4 or HDF file.

        Parameters
        ----------
        file_str: str
            The netCDF file path.
        drop_variables_str: str or iterable, optional
            List of variables to be dropped.
        decode_cf: bool, optional
            Whether to decode CF attributes and coordinate variables.
        decode_times: bool, optional
            Whether to decode time information (convert time coordinates to ``datetime`` objects).
        engine_str: str, optional
            Optional netCDF engine name.

        Return
        ------
        xarray.Dataset
        """
        ds = xr.open_dataset(file_str, drop_variables=drop_variables_str, decode_cf=decode_cf, decode_times=decode_times, engine=engine_str, chunks=1000000)
        return ds

    @classmethod
    def _load_virtual_variable(cls, ds, var_name):

        v_var = ds.variables[var_name]
        if v_var is not None and "virtual" in v_var.attrs:
            if not cls._is_already_loaded(v_var):
                dic = cls._create_dictionary_of_non_virtuals(ds)
                expression_ = v_var.attrs["expression"]
                biggest_variable = cls._get_biggest_variable(dic, expression_)
                dims = biggest_variable.dims

                to_extend = cls._find_used_one_dimensional_variables_to_extend(dic, dims, expression_)
                for name in to_extend:
                    dic[name] = cls._extend_1d_vertical_to_2d(dic[name], biggest_variable)

                to_interpolate = cls._find_used_tie_point_variables_to_extend(dic, expression_)
                for name in to_interpolate:
                    dic[name] = cls._interpolate_to_raster(dic[name], biggest_variable)

                expression_ = cls._replace_constants(expression_)
                values = ne.evaluate(expression_, dic)
                tmp_var = xr.Variable(dims, values)
                tmp_var.attrs = v_var.attrs
                ds._variables[var_name] = tmp_var
        else:
            raise IOError('no such virtual variable: "' + var_name + '"')

    @classmethod
    def _replace_constants(cls, expression_):
        _pattern_to_detect_pi = "\\b[Pp][Ii]\\b"
        return re.sub(_pattern_to_detect_pi, str(pi), expression_)

    @classmethod
    def _get_biggest_variable(cls, dic, expression):
        biggest_shape_product = 0
        biggest_var = None
        sorted_keys = cls._get_keys_sorted__longest_first(dic)
        for key in sorted_keys:
            if key in expression:
                expression = expression.replace(str(key), '')
                variable = dic[key]
                shape_product = cls._get_num_data_elements(variable)
                if shape_product > biggest_shape_product:
                    biggest_shape_product = shape_product
                    biggest_var = variable
        return biggest_var

    @classmethod
    def _get_num_data_elements(cls, variable):
        shape = variable.shape
        shape_len = len(shape)
        if shape_len == 0:
            return 1  # scalar variable

        num_elems = shape[0]
        for i in range(1, shape_len):
            num_elems = shape[i] * num_elems

        return num_elems

    @classmethod
    def _is_already_loaded(cls, variable):
        shape_len = len(variable.shape)
        if shape_len >= 1:
            return True

        return False

    @classmethod
    def _find_used_one_dimensional_variables_to_extend(cls, dic, biggest_dims, expression):
        one_dimensional_variables = list()
        dim_len = len(biggest_dims)
        if dim_len == 1:
            return one_dimensional_variables

        dim_of_interest = biggest_dims[dim_len - 2]
        sorted_keys = cls._get_keys_sorted__longest_first(dic)
        for key in sorted_keys:
            if key in expression:
                expression = expression.replace(str(key), '')
                variable = dic[key]
                if len(variable.shape) == 1 and variable.dims[0] == dim_of_interest:
                    one_dimensional_variables.append(key)
        return one_dimensional_variables

    @classmethod
    def _find_used_tie_point_variables_to_extend(cls, dic, expression):
        tie_point_variables = list()

        sorted_keys = cls._get_keys_sorted__longest_first(dic)
        for key in sorted_keys:
            if key in expression:
                expression = expression.replace(str(key), '')
                variable = dic[key]
                if "tie_points" in variable.attrs:
                    tie_point_variables.append(key)

        return tie_point_variables

    @classmethod
    def _get_keys_sorted__longest_first(cls, dic):
        dic_keys = dic.keys()
        return sorted(dic_keys, key=len, reverse=True)

    @classmethod
    def _extend_1d_vertical_to_2d(cls, vertical_variable, reference_var):
        shape = reference_var.shape[-2:]
        var_reshaped = np.resize(vertical_variable, shape[::-1])
        var_reshaped = np.moveaxis(var_reshaped, 0, 1)
        return xr.Variable(reference_var.dims[-2:], var_reshaped)

    @classmethod
    def _interpolate_to_raster(cls, variable, biggest_variable):
        shape = biggest_variable.shape
        full_size_array = resample_2d(variable.values, shape[1], shape[0])
        return xr.Variable(shape, full_size_array)

    @classmethod
    def _create_dictionary_of_non_virtuals(cls, ds):
        dic = {}
        for varName in ds.variables:
            if "virtual" not in ds.variables[varName].attrs:
                dic.update({varName: ds.variables[varName]})
        return dic
