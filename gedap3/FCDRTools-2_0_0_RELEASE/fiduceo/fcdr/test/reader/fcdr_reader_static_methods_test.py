import unittest as ut

import numpy as np
import xarray as xr

import fiduceo.fcdr.test.test_utils as ftu
from fiduceo.fcdr.reader.fcdr_reader import FCDRReader as R


class FCDRReaderStaticMethodsTest(ut.TestCase):

    def setUp(self):
        ds = xr.Dataset()
        ds['asd'] = xr.Variable(['z', 'y', 'x'], np.full([3, 4, 5], 13, np.int32))
        ds['lsmf'] = xr.Variable(['y', 'x'], np.full([4, 5], 14, np.int32))
        ds['bottle'] = xr.Variable(['y', 'x'], np.full([4, 5], 15, np.int32))
        ds['ka'] = xr.Variable(['x'], np.full([5], 16, np.int32))
        ds['ody'] = xr.Variable(['y'], range(17, 21))
        ds['woody'] = xr.Variable(['y'], np.full([4], 18, np.int32))
        ds['full'] = xr.Variable(['y_f', 'x_f'], np.full([10, 10], 19, np.int32))
        variable = xr.Variable(['y_tie', 'x_tie'], np.full([2, 2], 20, np.int32))
        variable.attrs["tie_points"] = "true"
        ds['tie'] = variable

        self.ds = ds
        self.dic = R._create_dictionary_of_non_virtuals(ds)

    def test_GetBiggestDimension(self):
        expression = '(asd)*[lsmf + bottle] - ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        self.assertEqual(("z", "y", "x"), biggest_variable.dims)

        expression = '[lsmf + bottle] - ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        self.assertEqual(("y", "x"), biggest_variable.dims)

        expression = 'lsmf +  ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        self.assertEqual(("y", "x"), biggest_variable.dims)

        expression = 'bottle - ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        self.assertEqual(("y", "x"), biggest_variable.dims)

        expression = 'ka * 4'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        self.assertEqual(("x",), biggest_variable.dims)

    def test_GetKeysSorted_LongestFirst(self):
        dictionary = {"a": 1, "a2": 2, "b": 3, "bb3": 4, "c2": 5}
        sorted_keys = R._get_keys_sorted__longest_first(dictionary)
        self.assertEqual(0, sorted_keys.index('bb3'))
        self.assertTrue(sorted_keys.index('c2') == 1 or sorted_keys.index('c2') == 2)
        self.assertTrue(sorted_keys.index('a2') == 1 or sorted_keys.index('a2') == 2)
        self.assertTrue(sorted_keys.index('a') == 3 or sorted_keys.index('a') == 4)
        self.assertTrue(sorted_keys.index('b') == 3 or sorted_keys.index('b') == 4)

    # @todo 3 tb/se complete that test 2018-06-28
    # def test_ExpandOneDimensionalVariables(self):
    #     pass

    def test_FindUsedOneDimensionalVerticalToExtend_emptyResult(self):
        """
        'ka' is a one dimensional variable but not y dimension.
        Only variables with 'y' dimension must be extended.
        """
        expression = 'asd * ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        dims = biggest_variable.dims
        one_d_vars = R._find_used_one_dimensional_variables_to_extend(self.dic, dims, expression)
        self.assertEqual(list(), one_d_vars)

    def test_FindUsedOneDimensionalVerticalToExtend_oneResult(self):
        """
        'ody' is a one dimensional variable with dimension 'y'
        """
        expression = 'asd * ody + ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        dims = biggest_variable.dims
        one_d_vars = R._find_used_one_dimensional_variables_to_extend(self.dic, dims, expression)
        self.assertEqual(1, len(one_d_vars))
        self.assertTrue('ody' in one_d_vars)
        variable = self.dic['ody']
        self.assertEqual(('y',), variable.dims)
        self.assertEqual(17, variable.data[0])

    def test_FindUsedOneDimensionalVerticalToExtend_twoResults(self):
        """
        'ody' and 'woody' are a one dimensional variables with dimension 'y'
        """
        expression = 'asd * (ody - woody) + ka'
        biggest_variable = R._get_biggest_variable(self.dic, expression)
        dims = biggest_variable.dims
        one_d_vars = R._find_used_one_dimensional_variables_to_extend(self.dic, dims, expression)
        self.assertEqual(2, len(one_d_vars))
        self.assertTrue('ody' in one_d_vars)
        self.assertTrue('woody' in one_d_vars)
        var_o = self.dic['ody']
        self.assertEqual(('y',), var_o.dims)
        self.assertEqual(17, var_o.data[0])
        var_w = self.dic['woody']
        self.assertEqual(('y',), var_w.dims)
        self.assertEqual(18, var_w.data[0])

    def test_find_used_tie_point_variables_none_used(self):
        expression = 'ody * woody'
        to_extend = R._find_used_tie_point_variables_to_extend(self.dic, expression)
        self.assertEqual(0, len(to_extend))

    def test_find_used_tie_point_variables_used(self):
        expression = 'full * tie'
        to_extend = R._find_used_tie_point_variables_to_extend(self.dic, expression)
        self.assertEqual(1, len(to_extend))
        self.assertEqual("tie", to_extend[0])

    def test_interpolate_to_raster(self):
        interpolated = R._interpolate_to_raster(self.ds["tie"], self.ds["full"])
        self.assertEqual((10, 10), interpolated.shape)
        self.assertAlmostEqual(20.0, interpolated.data[0, 0], 8)
        self.assertAlmostEqual(20.0, interpolated.data[4, 4], 8)
        self.assertAlmostEqual(20.0, interpolated.data[9, 9], 8)

    def test_extend_vertical_1D_variable_to_2D(self):
        vertical_variable = xr.Variable('y', [5, 6, 7])
        reference_variable = xr.Variable(('z', 'y', 'x'), [[[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6], ], [[4, 5, 6, 7], [5, 6, 7, 8], [6, 7, 8, 9], ], ])
        extended = R._extend_1d_vertical_to_2d(vertical_variable, reference_variable)
        self.assertEqual((3, 4), extended.shape)
        self.assertEqual(('y', 'x'), extended.dims)
        expected = np.asarray([[5, 5, 5, 5], [6, 6, 6, 6], [7, 7, 7, 7], ])
        ftu.assert_array_equals_with_index_error_message(self, expected, extended.data)
