import unittest

import numpy as np
import xarray as xr
from xarray import Variable

from fiduceo.fcdr.writer.global_flags import GlobalFlags as gf
from fiduceo.fcdr.writer.templates.mviri_flag_mapper import MVIRI_FlagMapper


class MVIRI_FlagMapperTest(unittest.TestCase):

    def setUp(self):
        self.mapper = MVIRI_FlagMapper()

        self.dataset = xr.Dataset()
        global_flags_data = xr.DataArray(np.full([3, 3], 0, np.uint8), dims=['y', 'x'])
        self.dataset["quality_pixel_bitmask"] = Variable(["y", "x"], global_flags_data)

        data_quality_data = xr.DataArray(np.full([3, 3], 0, np.uint8), dims=['y', 'x'])
        self.dataset["data_quality_bitmask"] = Variable(["y", "x"], data_quality_data)

    def test_map_global_flags_noSensorFlags(self):
        self.dataset["quality_pixel_bitmask"].data[1, 0] = gf.INVALID

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[1, 0])  # invalid
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[0, 1])

    def test_map_global_flags_uncertainty_suspicious(self):
        self.dataset["quality_pixel_bitmask"].data[2, 0] = gf.USE_WITH_CAUTION

        self.dataset["data_quality_bitmask"].data[1, 0] = 1  # uncertainty_suspicious

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 0])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[2, 0])  # use_with_caution
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[0, 1])

    def test_map_global_flags_uncertainty_too_large(self):
        self.dataset["quality_pixel_bitmask"].data[0, 1] = gf.INVALID_INPUT

        self.dataset["data_quality_bitmask"].data[2, 0] = 2  # uncertainty_too_large

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[1, 0])
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[2, 0])  # use_with_caution
        self.assertEqual(4, self.dataset["quality_pixel_bitmask"].data[0, 1])  # invalid_input

    def test_map_global_flags_space_view_suspicious(self):
        self.dataset["quality_pixel_bitmask"].data[1, 1] = gf.INVALID_GEOLOC

        self.dataset["data_quality_bitmask"].data[0, 1] = 4  # space_view_suspicious

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 1])  # use_with_caution
        self.assertEqual(8, self.dataset["quality_pixel_bitmask"].data[1, 1])  # invalid_geoloc

    def test_map_global_flags_not_on_earth(self):
        self.dataset["quality_pixel_bitmask"].data[1, 1] = gf.INVALID_GEOLOC

        self.dataset["data_quality_bitmask"].data[0, 1] = 8  # not_on_earth

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])
        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[0, 1])  # invalid
        self.assertEqual(8, self.dataset["quality_pixel_bitmask"].data[1, 1])  # invalid_geoloc

    def test_map_global_flags_suspect_time(self):
        self.dataset["quality_pixel_bitmask"].data[2, 1] = gf.INVALID_TIME

        self.dataset["data_quality_bitmask"].data[1, 1] = 16  # suspect_time

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[0, 1])
        self.assertEqual(18, self.dataset["quality_pixel_bitmask"].data[1, 1])  # use_with_caution AND invalid_time
        self.assertEqual(16, self.dataset["quality_pixel_bitmask"].data[2, 1])  # invalid_time

    def test_map_global_flags_suspect_geolocation(self):
        self.dataset["quality_pixel_bitmask"].data[0, 2] = gf.SENSOR_ERROR

        self.dataset["data_quality_bitmask"].data[2, 1] = 32  # suspect_geolocation

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[1, 1])
        self.assertEqual(10, self.dataset["quality_pixel_bitmask"].data[2, 1])  # use_with_caution AND invalid_geoloc
        self.assertEqual(32, self.dataset["quality_pixel_bitmask"].data[0, 2])  # sensor_error

    def test_map_global_all_flags(self):
        self.dataset["quality_pixel_bitmask"].data[0, 0] = gf.INVALID_GEOLOC
        self.dataset["quality_pixel_bitmask"].data[0, 1] = gf.INVALID_GEOLOC
        self.dataset["quality_pixel_bitmask"].data[0, 2] = gf.INVALID_GEOLOC
        self.dataset["quality_pixel_bitmask"].data[1, 0] = gf.INVALID_GEOLOC
        self.dataset["quality_pixel_bitmask"].data[1, 1] = gf.INVALID_GEOLOC

        self.dataset["data_quality_bitmask"].data[0, 0] = 3 # uncertainty_suspicious & uncertainty_too_large
        self.dataset["data_quality_bitmask"].data[0, 1] = 20 # space_view_suspicious & suspect_time
        self.dataset["data_quality_bitmask"].data[0, 2] = 33 # suspect_geolocation & uncertainty_suspicious

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(10, self.dataset["quality_pixel_bitmask"].data[0, 0]) # use_with_caution & invalid_geoloc
        self.assertEqual(26, self.dataset["quality_pixel_bitmask"].data[0, 1])  # invalid_input & invalid_geoloc & invalid_time
        self.assertEqual(10, self.dataset["quality_pixel_bitmask"].data[0, 2])  # invalid_input & invalid_geoloc
        self.assertEqual(8, self.dataset["quality_pixel_bitmask"].data[1, 0])  # invalid_geoloc
        self.assertEqual(8, self.dataset["quality_pixel_bitmask"].data[1, 1])  # invalid_geoloc
