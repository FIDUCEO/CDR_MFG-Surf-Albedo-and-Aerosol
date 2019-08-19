import unittest

import numpy as np
import xarray as xr
from xarray import Variable

from fiduceo.fcdr.writer.global_flags import GlobalFlags as gf
from fiduceo.fcdr.writer.templates.hirs_flag_mapper import HIRS_FlagMapper


class HIRS_FlagMapperTest(unittest.TestCase):

    def setUp(self):
        self.mapper = HIRS_FlagMapper()

        self.dataset = xr.Dataset()
        global_flags_data = xr.DataArray(np.full([3, 3], 0, np.uint8), dims=['y', 'x'])
        self.dataset["quality_pixel_bitmask"] = Variable(["y", "x"], global_flags_data)

        data_quality_data = xr.DataArray(np.full([3, 3], 0, np.uint8), dims=['y', 'x'])
        self.dataset["data_quality_bitmask"] = Variable(["y", "x"], data_quality_data)

        data_quality_data = xr.DataArray(np.full([3, ], 0, np.int32), dims=['y'])
        self.dataset["quality_scanline_bitmask"] = Variable(["y"], data_quality_data)

        data_quality_data = xr.DataArray(np.full([3, 19], 0, np.uint8), dims=['y', 'channel'])
        self.dataset["quality_channel_bitmask"] = Variable(["y", "channel"], data_quality_data)

    def test_map_global_flags_noSensorFlags(self):
        self.dataset["quality_pixel_bitmask"].data[1, 0] = gf.INVALID

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[1, 0])  # invalid
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[0, 1])

    def test_map_global_flags_suspect_mirror(self):
        self.dataset["quality_pixel_bitmask"].data[2, 1] = gf.PADDED_DATA

        self.dataset["data_quality_bitmask"].data[1, 0] = 1  # suspect_mirror

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 0])  # use_with_caution
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 0])
        self.assertEqual(64, self.dataset["quality_pixel_bitmask"].data[2, 1])

    def test_map_global_flags_suspect_time(self):
        self.dataset["quality_pixel_bitmask"].data[2, 2] = gf.INCOMPLETE_CHANNEL_DATA

        self.dataset["data_quality_bitmask"].data[1, 1] = 4  # suspect_time

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[1, 0])
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 1])  # use_with_caution
        self.assertEqual(128, self.dataset["quality_pixel_bitmask"].data[2, 2])

    def test_map_global_flags_reduced_context_scanline(self):
        self.dataset["quality_pixel_bitmask"].data[2, 0] = gf.INCOMPLETE_CHANNEL_DATA
        self.dataset["quality_pixel_bitmask"].data[0, 0] = gf.PADDED_DATA

        self.dataset["quality_scanline_bitmask"].data[1] = 536870912  # reduced_context

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(64, self.dataset["quality_pixel_bitmask"].data[0, 0])  # padded_data
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 0])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 1])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 2])  # use_with_caution
        self.assertEqual(128, self.dataset["quality_pixel_bitmask"].data[2, 0])  # incomplete_channel_data

    def test_map_global_flags_no_rself_scanline(self):
        self.dataset["quality_pixel_bitmask"].data[0, 1] = gf.PADDED_DATA
        self.dataset["quality_pixel_bitmask"].data[1, 2] = gf.SENSOR_ERROR

        self.dataset["quality_scanline_bitmask"].data[1] = 1073741824  # bad_temp_no_rself

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(64, self.dataset["quality_pixel_bitmask"].data[0, 1])  # padded_data
        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[1, 0])  # invalid
        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[1, 1])  # invalid
        self.assertEqual(33, self.dataset["quality_pixel_bitmask"].data[1, 2])  # invalid & sensor_error
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 2])

    def test_map_global_flags_do_not_use_single_channel(self):
        self.dataset["quality_pixel_bitmask"].data[1, 1] = gf.INCOMPLETE_CHANNEL_DATA
        self.dataset["quality_pixel_bitmask"].data[2, 2] = gf.INVALID

        self.dataset["quality_channel_bitmask"].data[0, 1] = 1  # do_not_use in channel 2

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 0])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 1])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 2])  # use_with_caution
        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[2, 2])  # invalid
        self.assertEqual(128, self.dataset["quality_pixel_bitmask"].data[1, 1])  # incomplete_channel_data

    def test_map_global_flags_do_not_use_all_channels(self):
        self.dataset["quality_pixel_bitmask"].data[1, 1] = gf.USE_WITH_CAUTION
        self.dataset["quality_pixel_bitmask"].data[2, 2] = gf.SENSOR_ERROR

        self.dataset["quality_channel_bitmask"].data[2, :] = 1  # do_not_use in all channels

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[2, 0])  # invalid
        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[2, 1])  # invalid
        self.assertEqual(33, self.dataset["quality_pixel_bitmask"].data[2, 2])  # invalid | sensor_error
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[0, 1])
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 1])  # use_with_caution

    def test_map_global_flags_uncertainty_suspicious_single_channel(self):
        self.dataset["quality_pixel_bitmask"].data[1, 1] = gf.INVALID_INPUT
        self.dataset["quality_pixel_bitmask"].data[2, 2] = gf.INVALID

        self.dataset["quality_channel_bitmask"].data[0, 2] = 2  # uncertainty_suspicious in channel 3

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 0])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 1])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 2])  # use_with_caution
        self.assertEqual(1, self.dataset["quality_pixel_bitmask"].data[2, 2])  # invalid
        self.assertEqual(4, self.dataset["quality_pixel_bitmask"].data[1, 1])  # invalid_input

    def test_map_global_all_flags(self):
        self.dataset["quality_pixel_bitmask"].data[2, 0] = gf.INVALID_GEOLOC
        self.dataset["quality_pixel_bitmask"].data[0, 1] = gf.INVALID_GEOLOC

        self.dataset["data_quality_bitmask"].data[2, 0] = 3  # suspect_mirror & suspect_geo

        self.dataset["quality_channel_bitmask"].data[1, :] = 2  # uncertainty_suspicious in all channels

        self.dataset["quality_scanline_bitmask"].data[0] = 536870912  # reduced_context

        self.mapper.map_global_flags(self.dataset)

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 0])  # use_with_caution
        self.assertEqual(10, self.dataset["quality_pixel_bitmask"].data[0, 1])  # use_with_caution & invalid_geoloc
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[0, 2])  # use_with_caution

        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 0])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 1])  # use_with_caution
        self.assertEqual(2, self.dataset["quality_pixel_bitmask"].data[1, 2])  # use_with_caution

        self.assertEqual(10, self.dataset["quality_pixel_bitmask"].data[2, 0])  # use_with_caution & invalid_geoloc
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 1])
        self.assertEqual(0, self.dataset["quality_pixel_bitmask"].data[2, 2])
