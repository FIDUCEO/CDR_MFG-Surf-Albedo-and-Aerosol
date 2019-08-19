import unittest

import numpy as np
from xarray import Variable

from fiduceo.fcdr.writer.data_utility import DataUtility
from fiduceo.common.writer.default_data import DefaultData


class DataUtilityTest(unittest.TestCase):

    def test_check_scaling_ranges_int8_vector_ok(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 11.872  # -128
        default_array[1] = 12.127  # 127
        default_array[2] = np.NaN
        default_array[3] = 12.04
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('scale_factor', 0.001), ('add_offset', 12)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_int8_vector_underflow(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 11.702  # underflow
        default_array[1] = 12.127  # 127
        default_array[2] = np.NaN
        default_array[3] = 12.04
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('scale_factor', 0.001), ('add_offset', 12)])
        
        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_uint8_vector_ok(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 11.0    # 0
        default_array[1] = 16.1    # 255
        default_array[2] = np.NaN
        default_array[3] = 13.05
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.uint8), ('_FillValue', 255), ('scale_factor', 0.02), ('add_offset', 11)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_uint8_vector_overflow(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 11.0    # 0
        default_array[1] = 16.2    # overflow
        default_array[2] = np.NaN
        default_array[3] = 13.05
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.uint8), ('_FillValue', 255), ('scale_factor', 0.02), ('add_offset', 11)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_int16_array_ok(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 75.534  # 32767
        default_array[0][1] = -55.536 # -32768
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.encoding = dict([('dtype', np.int16), ('_FillValue', -32767), ('scale_factor', 0.002), ('add_offset', 10)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_int16_array_underflow(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 75.534  # 32767
        default_array[0][1] = -55.636 # -32768
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.encoding = dict([('dtype', np.int16), ('_FillValue', -32767), ('scale_factor', 0.002), ('add_offset', 10)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_uint16_array_ok(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 9       # 9
        default_array[0][1] = 205.605 # 65535
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.encoding = dict([('dtype', np.uint16), ('_FillValue', 65535), ('scale_factor', 0.003), ('add_offset', 9)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_uint16_array_overflow(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 9       # 0
        default_array[0][1] = 205.705 # overflow
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.encoding = dict([('dtype', np.uint16), ('_FillValue', 65535), ('scale_factor', 0.003), ('add_offset', 9)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_int32_vector_ok(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 21487.83647   # 2147483647
        default_array[1] = -21461.83648  # -2147483648
        default_array[2] = np.NaN
        default_array[3] = 14.04
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int32), ('_FillValue', -2147483647), ('scale_factor', 0.00001), ('add_offset', 13)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_int32_vector_underflow(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 21487.83647   # 2147483647
        default_array[1] = -21461.93648  # underflow
        default_array[2] = np.NaN
        default_array[3] = 14.04
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int32), ('_FillValue', -2147483647), ('scale_factor', 0.00001), ('add_offset', 13)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_uint32_vector_ok(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 14          # 0
        default_array[1] = 85913.3459  # 4294967295
        default_array[2] = np.NaN
        default_array[3] = 14.01
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.uint32), ('_FillValue', 4294967295), ('scale_factor', 0.00002), ('add_offset', 14)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_uint32_vector_overflow(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = 14          # 0
        default_array[1] = 85913.4459  # 4294967295
        default_array[2] = np.NaN
        default_array[3] = 14.01
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.uint32), ('_FillValue', 4294967295), ('scale_factor', 0.00002), ('add_offset', 14)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_uint32_only_NaN(self):
        default_array = DefaultData.create_default_vector(4, np.float32)
        default_array[0] = np.NaN
        default_array[1] = np.NaN
        default_array[2] = np.NaN
        default_array[3] = np.NaN
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.uint32), ('_FillValue', 4294967295), ('scale_factor', 0.00002), ('add_offset', 14)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_int16_valid_min_max_underflow(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 60  # 25000
        default_array[0][1] = 9   # underflow
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.attrs["valid_max"] = 25000
        variable.attrs["valid_min"] = 0
        variable.encoding = dict([('dtype', np.int16), ('_FillValue', -32767), ('scale_factor', 0.002), ('add_offset', 10)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_int16_valid_min_max_overflow(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 61  # overflow
        default_array[0][1] = 10  # 0
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.attrs["valid_max"] = 25000
        variable.attrs["valid_min"] = 0
        variable.encoding = dict([('dtype', np.int16), ('_FillValue', -32767), ('scale_factor', 0.002), ('add_offset', 10)])

        try:
            DataUtility.check_scaling_ranges(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass

    def test_check_scaling_ranges_int16_valid_min_max_ok(self):
        default_array = DefaultData.create_default_array(2, 2, np.float32)
        default_array[0][0] = 60  # 25000
        default_array[0][1] = 10  # 0
        default_array[1][0] = np.NaN
        default_array[1][1] = 14.06
        variable = Variable(["y", "x"], default_array)
        variable.attrs["valid_max"] = 25000
        variable.attrs["valid_min"] = 0
        variable.encoding = dict([('dtype', np.int16), ('_FillValue', -32767), ('scale_factor', 0.002), ('add_offset', 10)])

        DataUtility.check_scaling_ranges(variable)

    def test_check_scaling_ranges_int8_vector_unscaled(self):
        default_array = DefaultData.create_default_vector(4, np.int8)
        default_array[0] = 108
        default_array[1] = 109
        default_array[2] = -127
        default_array[3] = 110
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('_FillValue', -127)])

        DataUtility.check_scaling_ranges(variable)

    def test__get_scale_factor(self):
        default_array = DefaultData.create_default_vector(2, np.float32)
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('scale_factor', 0.023), ('add_offset', 12)])

        scale_factor = DataUtility._get_scale_factor(variable)
        self.assertEqual(0.023, scale_factor)

    def test__get_scale_factor_missing(self):
        default_array = DefaultData.create_default_vector(2, np.float32)
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('add_offset', 12)])

        scale_factor = DataUtility._get_scale_factor(variable)
        self.assertEqual(1.0, scale_factor)

    def test__get_add_offset(self):
        default_array = DefaultData.create_default_vector(2, np.float32)
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('scale_factor', 0.023), ('add_offset', 1.98)])

        add_offset = DataUtility._get_add_offset(variable)
        self.assertEqual(1.98, add_offset)

    def test__get_add_offset_missing(self):
        default_array = DefaultData.create_default_vector(2, np.float32)
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('scale_factor', 0.023)])

        add_offset = DataUtility._get_add_offset(variable)
        self.assertEqual(0.0, add_offset)

    def test__get_min_max(self):
        default_array = DefaultData.create_default_vector(2, np.float32)
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('dtype', np.int8), ('_FillValue', -127), ('scale_factor', 0.023)])

        min_max = DataUtility._get_min_max(variable)
        self.assertEqual(-128, min_max.min)
        self.assertEqual(127, min_max.max)

    def test__get_min_max_missing_type(self):
        default_array = DefaultData.create_default_vector(2, np.float32)
        variable = Variable(["y"], default_array)
        variable.encoding = dict([('_FillValue', -127), ('scale_factor', 0.023)])

        try:
            DataUtility._get_min_max(variable)
            self.fail("ValueError expected")
        except ValueError:
            pass