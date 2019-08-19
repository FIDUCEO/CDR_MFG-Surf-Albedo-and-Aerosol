import unittest

import numpy as np

from fiduceo.common.writer.default_data import DefaultData


class DefaultDataTest(unittest.TestCase):
    def test_create_default_vector_int(self):
        default_array = DefaultData.create_default_vector(14, np.int32)
        self.assertEqual((14,), default_array.shape)
        self.assertEqual(-2147483647, default_array.data[2])

    def test_create_default_vector_int_fill_value(self):
        default_array = DefaultData.create_default_vector(15, np.int32, fill_value=108)
        self.assertEqual((15,), default_array.shape)
        self.assertEqual(108, default_array.data[3])

    def test_get_default_fill_value(self):
        self.assertEqual(-127, DefaultData.get_default_fill_value(np.int8))
        self.assertEqual(-32767, DefaultData.get_default_fill_value(np.int16))
        self.assertEqual(np.uint16(-1), DefaultData.get_default_fill_value(np.uint16))
        self.assertEqual(-2147483647, DefaultData.get_default_fill_value(np.int32))
        self.assertEqual(-9223372036854775806, DefaultData.get_default_fill_value(np.int64))
        self.assertEqual(np.float32(9.96921E36), DefaultData.get_default_fill_value(np.float32))
        self.assertEqual(9.969209968386869E36, DefaultData.get_default_fill_value(np.float64))

