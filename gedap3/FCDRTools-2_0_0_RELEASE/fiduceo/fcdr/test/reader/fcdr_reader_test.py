import unittest as ut

import numpy as np
import xarray as xr

import fiduceo.fcdr.test.test_utils as tu
from fiduceo.fcdr.reader.fcdr_reader import FCDRReader


class FCDRReaderTest(ut.TestCase):

    def setUp(self):
        self.fcdr_reader = FCDRReader()

    def test_adding_2_three_dimensional_variables(self):
        ds = xr.Dataset()
        ds['a'] = create_three_dim_variable()
        ds['b'] = create_three_dim_variable()
        v_var = create_virtual_variable("a + b")
        ds["v_var"] = v_var

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

        self.assertTrue('v_var' in ds)
        virtual_loaded = ds['v_var']
        self.assertEqual(('z', 'y', 'x'), virtual_loaded.dims)
        self.assertEqual((4, 2, 3), virtual_loaded.shape)
        self.assertEqual(v_var.attrs, virtual_loaded.attrs)

        expected = np.asarray([[[2, 4, 6], [2.2, 4.2, 6.2]], [[22, 24, 26], [22.2, 24.2, 26.2]], [[42, 44, 46], [42.2, 44.2, 46.2]], [[62, 64, 66], [62.2, 64.2, 66.2]]])
        actual = ds['v_var'].values

        self.assertEqual(type(expected), type(actual))
        tu.assert_array_equals_with_index_error_message(self, expected, actual)

    def test_adding_1_three_dimensional_variable_and_1_two_dimensional_variable(self):
        ds = xr.Dataset()
        ds['a'] = create_three_dim_variable()
        ds['b'] = create_two_dim_variable()
        v_var = create_virtual_variable("a + b")
        ds["v_var"] = v_var

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

        self.assertTrue('v_var' in ds)
        virtual_loaded = ds['v_var']
        self.assertEqual(('z', 'y', 'x'), virtual_loaded.dims)
        self.assertEqual((4, 2, 3), virtual_loaded.shape)
        self.assertEqual(v_var.attrs, virtual_loaded.attrs)

        expected = np.asarray([[[2, 4, 6], [5.1, 7.1, 9.1]], [[12, 14, 16], [15.1, 17.1, 19.1]], [[22, 24, 26], [25.1, 27.1, 29.1]], [[32, 34, 36], [35.1, 37.1, 39.1]]])
        actual = virtual_loaded.values

        self.assertEqual(type(expected), type(actual))
        tu.assert_array_equals_with_index_error_message(self, expected, actual)

    def test_adding_1_three_dimensional_variable_and_1_one_dimensional_variable(self):
        ds = xr.Dataset()
        ds['a'] = create_three_dim_variable()
        ds['b'] = create_one_dim_variable()
        v_var = create_virtual_variable("a + b")
        ds["v_var"] = v_var

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

        self.assertTrue('v_var' in ds)
        virtual_loaded = ds['v_var']
        self.assertEqual(('z', 'y', 'x'), virtual_loaded.dims)
        self.assertEqual((4, 2, 3), virtual_loaded.shape)
        self.assertEqual(v_var.attrs, virtual_loaded.attrs)

        expected = np.asarray([[[2.3, 4.3, 6.3], [2.4, 4.4, 6.4]], [[12.3, 14.3, 16.3], [12.4, 14.4, 16.4]], [[22.3, 24.3, 26.3], [22.4, 24.4, 26.4]], [[32.3, 34.3, 36.3], [32.4, 34.4, 36.4]]])
        actual = virtual_loaded.values

        self.assertEqual(type(expected), type(actual))
        tu.assert_array_equals_with_index_error_message(self, expected, actual)

    def test_adding_1_three_dimensional_variable_and_1_vertical_one_dimensional_variable(self):
        ds = xr.Dataset()
        ds['a'] = create_three_dim_variable()
        ds['b'] = create_vertical_one_dim_variable()
        v_var = create_virtual_variable("a + b")
        ds["v_var"] = v_var

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

        self.assertTrue('v_var' in ds)
        virtual_loaded = ds['v_var']
        self.assertEqual(('z', 'y', 'x'), virtual_loaded.dims)
        self.assertEqual((4, 2, 3), virtual_loaded.shape)
        self.assertEqual(v_var.attrs, virtual_loaded.attrs)

        expected = np.asarray([[[6, 7, 8], [7.1, 8.1, 9.1]], [[16, 17, 18], [17.1, 18.1, 19.1]], [[26, 27, 28], [27.1, 28.1, 29.1]], [[36, 37, 38], [37.1, 38.1, 39.1]]])
        actual = virtual_loaded.values

        self.assertEqual(type(expected), type(actual))
        tu.assert_array_equals_with_index_error_message(self, expected, actual)

    def test_multiplying_1_three_dimensional_variable_and_1_scalar_variable(self):
        ds = xr.Dataset()
        ds['a'] = create_three_dim_variable()
        ds['b'] = create_scalar_variable(2.7)
        v_var = create_virtual_variable("a * b")
        ds["v_var"] = v_var

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

        self.assertTrue('v_var' in ds)
        virtual_loaded = ds['v_var']
        self.assertEqual(('z', 'y', 'x'), virtual_loaded.dims)
        self.assertEqual((4, 2, 3), virtual_loaded.shape)
        self.assertEqual(v_var.attrs, virtual_loaded.attrs)

        expected = np.asarray(
            [[[2.7, 5.4, 8.1], [2.97, 5.67, 8.37]], [[29.7, 32.4, 35.1], [29.97, 32.67, 35.37]], [[56.7, 59.4, 62.1], [56.97, 59.67, 62.37]], [[83.7, 86.4, 89.1], [83.97, 86.67, 89.37]]])
        actual = virtual_loaded.values

        self.assertEqual(type(expected), type(actual))
        tu.assert_array_equals_with_index_error_message(self, expected, actual)

    def test_multiplying_2d_array_and_tiepoint_array(self):
        ds = xr.Dataset()
        ds['a'] = create_two_dim_ones(10)
        ds['b'] = create_two_dim_tie_points_variable()
        ds["v_var"] = create_virtual_variable("a * b")

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

    def test_prepare_virtual_variables(self):
        ds = xr.Dataset()
        ds['a'] = create_three_dim_variable()
        ds['b'] = create_two_dim_variable()
        ds["v_var"] = create_virtual_variable("a + b")

        self.fcdr_reader._load_virtual_variable(ds, 'v_var')

        self.assertEqual(3, len(ds.data_vars))
        expected = np.full([4, 2, 3], 1, np.int32)
        actual = ds['v_var'].values
        self.assertEqual(type(expected), type(actual))
        self.assertEqual(expected.shape, actual.shape)

    def test_is_already_loaded(self):
        virtual_variable = create_virtual_variable("alpha * beta")
        self.assertFalse(self.fcdr_reader._is_already_loaded(virtual_variable))

        variable = create_three_dim_variable()
        self.assertTrue(self.fcdr_reader._is_already_loaded(variable))

    def test_get_num_data_elements_virtual(self):
        virtual_variable = create_virtual_variable("alpha * beta")
        self.assertEqual(1, self.fcdr_reader._get_num_data_elements(virtual_variable))

    def test_get_num_data_elements_one_d(self):
        one_d_variable = create_one_dim_variable()
        self.assertEqual(3, self.fcdr_reader._get_num_data_elements(one_d_variable))

    def test_get_num_data_elements_two_d(self):
        two_d_variable = create_two_dim_variable()
        self.assertEqual(6, self.fcdr_reader._get_num_data_elements(two_d_variable))

    def test_get_num_data_elements_three_d(self):
        three_d_variable = create_three_dim_variable()
        self.assertEqual(24, self.fcdr_reader._get_num_data_elements(three_d_variable))


def create_three_dim_variable():
    return xr.Variable(['z', 'y', 'x'],
                       np.asarray([[[1.0, 2.0, 3.0], [1.1, 2.1, 3.1]], [[11.0, 12.0, 13.0], [11.1, 12.1, 13.1]], [[21.0, 22.0, 23.0], [21.1, 22.1, 23.1]], [[31.0, 32.0, 33.0], [31.1, 32.1, 33.1]]]))


def create_two_dim_variable():
    return xr.Variable(['y', 'x'], np.asarray([[1, 2, 3], [4, 5, 6]]))

def create_two_dim_tie_points_variable():
    variable = xr.Variable(['y_tp', 'x_tp'], np.asarray([[1, 2], [3, 4]]))
    variable.attrs["tie_points"] = "true"
    return variable

def create_two_dim_ones(size):
    return xr.Variable(['y', 'x'], np.ones((size, size)))

def create_one_dim_variable():
    return xr.Variable(['x'], np.asarray([1.3, 2.3, 3.3]))


def create_vertical_one_dim_variable():
    return xr.Variable(['y'], np.asarray([5, 6], dtype=np.int32))


def create_scalar_variable(value):
    return xr.Variable(['s'], np.asarray([value]))


def create_virtual_variable(expression):
    v_var = xr.Variable([], np.float32(1.0))
    v_var.attrs["virtual"] = "true"
    v_var.attrs["expression"] = expression
    return v_var
