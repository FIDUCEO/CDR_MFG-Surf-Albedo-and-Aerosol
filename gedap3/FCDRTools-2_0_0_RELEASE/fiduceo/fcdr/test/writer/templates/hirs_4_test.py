import unittest

import xarray as xr

from fiduceo.common.test.assertions import Assertions
from fiduceo.fcdr.test.writer.templates.hirs_assert import HIRSAssert
from fiduceo.fcdr.writer.templates.hirs_4 import HIRS4

CHUNKING_3D = (10, 512, 56)
CHUNKING_2D = (512, 56)

NUM_CHANNELS = 19


class HIRS4Test(unittest.TestCase):
    def test_add_original_variables(self):
        ha = HIRSAssert()
        ds = xr.Dataset()
        HIRS4.add_original_variables(ds, 6)

        Assertions.assert_geolocation_variables(self, ds, 56, 6, chunking=CHUNKING_2D)
        Assertions.assert_quality_flags(self, ds, 56, 6, chunking=CHUNKING_2D)

        ha.assert_bt_variable(ds, chunking=CHUNKING_3D)
        ha.assert_common_angles(ds, chunking=CHUNKING_2D)
        ha.assert_common_sensor_variables(ds, 2751)
        ha.assert_extended_quality_flags(ds)
        ha.assert_coordinates(ds)

    def test_get_swath_width(self):
        self.assertEqual(56, HIRS4.get_swath_width())

    def test_add_easy_fcdr_variables(self):
        delta_x = 16
        delta_y = 17
        ha = HIRSAssert()
        ds = xr.Dataset()
        HIRS4.add_easy_fcdr_variables(ds, 7, lut_size=24, corr_dx=delta_x, corr_dy=delta_y)

        ha.assert_easy_fcdr_uncertainties(ds, chunking=CHUNKING_3D)

        Assertions.assert_correlation_matrices(self, ds, NUM_CHANNELS)
        Assertions.assert_lookup_tables(self, ds, NUM_CHANNELS, 24)
        Assertions.assert_correlation_coefficients(self, ds, NUM_CHANNELS, delta_x, delta_y)

    def test_add_full_fcdr_variables(self):
        ha = HIRSAssert()
        ds = xr.Dataset()
        HIRS4.add_full_fcdr_variables(ds, 6)
        ha.assert_minor_frame_flags(ds)

    def test_add_template_key(self):
        ds = xr.Dataset()

        HIRS4.add_template_key(ds)

        self.assertEqual("HIRS4", ds.attrs["template_key"])
