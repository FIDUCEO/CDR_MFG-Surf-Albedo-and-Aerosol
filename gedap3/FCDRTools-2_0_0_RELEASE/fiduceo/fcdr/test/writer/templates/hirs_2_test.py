import unittest

import numpy as np
import xarray as xr

from fiduceo.common.test.assertions import Assertions
from fiduceo.common.writer.default_data import DefaultData
from fiduceo.fcdr.test.writer.templates.hirs_assert import HIRSAssert
from fiduceo.fcdr.writer.templates.hirs_2 import HIRS2

CHUNKING_3D = (10, 512, 56)
CHUNKING_2D = (512, 56)
NUM_CHANNELS = 19


class HIRS2Test(unittest.TestCase):
    def test_add_original_variables(self):
        ha = HIRSAssert()
        ds = xr.Dataset()
        HIRS2.add_original_variables(ds, 6)

        Assertions.assert_geolocation_variables(self, ds, 56, 6, chunking=CHUNKING_2D)
        Assertions.assert_quality_flags(self, ds, 56, 6, chunking=CHUNKING_2D)

        ha.assert_bt_variable(ds, chunking=CHUNKING_3D)
        self._assert_angle_variables(ds)
        ha.assert_common_sensor_variables(ds, 102)
        ha.assert_coordinates(ds)

    def test_get_swath_width(self):
        self.assertEqual(56, HIRS2.get_swath_width())

    def test_add_easy_fcdr_variables(self):
        delta_x = 14
        delta_y = 15
        ha = HIRSAssert()
        ds = xr.Dataset()

        HIRS2.add_easy_fcdr_variables(ds, 7, lut_size=22, corr_dx=delta_x, corr_dy=delta_y)

        ha.assert_easy_fcdr_uncertainties(ds, chunking=CHUNKING_3D)

        Assertions.assert_correlation_matrices(self, ds, NUM_CHANNELS)
        Assertions.assert_lookup_tables(self, ds, NUM_CHANNELS, 22)
        Assertions.assert_correlation_coefficients(self, ds, NUM_CHANNELS, delta_x, delta_y)

    def test_add_full_fcdr_variables(self):
        # @todo 2 tb/tb add something here
        pass

    def test_add_template_key(self):
        ds = xr.Dataset()

        HIRS2.add_template_key(ds)

        self.assertEqual("HIRS2", ds.attrs["template_key"])

    def _assert_angle_variables(self, ds):
        satellite_zenith_angle = ds.variables["satellite_zenith_angle"]
        self.assertEqual((6,), satellite_zenith_angle.shape)
        self.assertTrue(np.isnan(satellite_zenith_angle.data[3]))
        self.assertEqual(np.uint16, satellite_zenith_angle.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), satellite_zenith_angle.encoding['_FillValue'])
        self.assertEqual(0.01, satellite_zenith_angle.encoding['scale_factor'])
        self.assertEqual(-180.0, satellite_zenith_angle.encoding['add_offset'])
        self.assertEqual("platform_zenith_angle", satellite_zenith_angle.attrs["standard_name"])
        self.assertEqual("degree", satellite_zenith_angle.attrs["units"])
        self.assertEqual("longitude latitude", satellite_zenith_angle.attrs["coordinates"])

        solar_azimuth_angle = ds.variables["solar_azimuth_angle"]
        self.assertEqual((6, 56), solar_azimuth_angle.shape)
        self.assertTrue(np.isnan(solar_azimuth_angle.data[4, 4]))
        self.assertEqual(np.uint16, solar_azimuth_angle.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), solar_azimuth_angle.encoding['_FillValue'])
        self.assertEqual(0.01, solar_azimuth_angle.encoding['scale_factor'])
        self.assertEqual(-180.0, solar_azimuth_angle.encoding['add_offset'])
        self.assertEqual(CHUNKING_2D, solar_azimuth_angle.encoding['chunksizes'])
        self.assertEqual("solar_azimuth_angle", solar_azimuth_angle.attrs["standard_name"])
        self.assertEqual("degree", solar_azimuth_angle.attrs["units"])
        self.assertEqual("longitude latitude", solar_azimuth_angle.attrs["coordinates"])
