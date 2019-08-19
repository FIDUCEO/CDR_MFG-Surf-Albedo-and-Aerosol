import unittest

import numpy as np
import xarray as xr
from fiduceo.fcdr.writer.global_flags import GlobalFlags as gf
from xarray import Variable

from fiduceo.fcdr.writer.templates.avhrr_flag_mapper import AVHRR_FlagMapper


class AVHRR_FlagMapperTest(unittest.TestCase):

    def setUp(self):
        self.mapper = AVHRR_FlagMapper()

        self.dataset = xr.Dataset()
        global_flags_data = xr.DataArray(np.full([3, 3], 0, np.uint8), dims=['y', 'x'])
        self.dataset["quality_pixel_bitmask"] = Variable(["y", "x"], global_flags_data)

        data_quality_data = xr.DataArray(np.full([3, 3], 0, np.uint8), dims=['y', 'x'])
        self.dataset["data_quality_bitmask"] = Variable(["y", "x"], data_quality_data)

    def test_map_global_flags_noSensorFlags(self):
        self.dataset["quality_pixel_bitmask"].data[0, 0] = gf.INVALID

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[0, 0])  # invalid
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[1, 0])
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])

    def test_map_global_flags_bad_geolocation_timing_err(self):
        self.dataset["quality_pixel_bitmask"].data[1, 0] = gf.USE_WITH_CAUTION

        self.dataset["data_quality_bitmask"].data[0, 0] = 1

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 0]) # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 0]) # use_with_caution
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])

    def test_map_global_flags_bad_calibration_radiometer_err(self):
        self.dataset["quality_pixel_bitmask"].data[2, 0] = gf.INVALID_INPUT

        self.dataset["data_quality_bitmask"].data[1, 0] = 2

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[0, 0])
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 0]) # use_with_caution
        self.assertEqual(4, self.dataset["quality_pixel_bitmask"].data[2, 0]) # invalid_input

    def test_map_global_all_flags(self):
        self.dataset["quality_pixel_bitmask"].data[2, 0] = gf.INVALID_GEOLOC
        self.dataset["quality_pixel_bitmask"].data[0, 1] = gf.INVALID_GEOLOC

        self.dataset["data_quality_bitmask"].data[2, 0] = 3

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[1, 0])
        self.assertEqual(10, self.dataset["quality_pixel_bitmask"].data[2, 0]) # invalid_input & invalid_geoloc
        self.assertEqual(8, self.dataset["quality_pixel_bitmask"].data[0, 1])  # invalid_input
