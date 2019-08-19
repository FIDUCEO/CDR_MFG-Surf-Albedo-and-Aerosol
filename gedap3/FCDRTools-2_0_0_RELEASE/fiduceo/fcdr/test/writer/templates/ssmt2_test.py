import unittest

import numpy as np
import xarray as xr

from fiduceo.common.test.assertions import Assertions
from fiduceo.fcdr.writer.templates.ssmt2 import SSMT2


class SSMT2Test(unittest.TestCase):
    def test_add_original_variables(self):
        ds = xr.Dataset()
        SSMT2.add_original_variables(ds, 4)

        Assertions.assert_geolocation_variables(self, ds, 28, 4)
        Assertions.assert_quality_flags(self, ds, 28, 4)

        temp_misc_hk = ds.variables["Temperature_misc_housekeeping"]
        self.assertEqual((18, 4), temp_misc_hk.shape)
        self.assertTrue(np.isnan(temp_misc_hk.data[1, 3]))
        self.assertTrue(np.isnan(temp_misc_hk.attrs['_FillValue']))
        self.assertEqual("TODO", temp_misc_hk.attrs["long_name"])
        self.assertEqual("TODO", temp_misc_hk.attrs["units"])

        ancil_data = ds.variables["ancil_data"]
        self.assertEqual((10, 4), ancil_data.shape)
        self.assertTrue(np.isnan(ancil_data.data[2, 3]))
        self.assertTrue(np.isnan(ancil_data.attrs['_FillValue']))
        self.assertEqual("Additional per scan information: year, day_of_year, secs_of_day, sat_lat, sat_long, sat_alt, sat_heading, year, day_of_year, secs_of_day", ancil_data.attrs["long_name"])

        ch_qual_flag = ds.variables["channel_quality_flag"]
        self.assertEqual((5, 4, 28), ch_qual_flag.shape)
        self.assertTrue(np.isnan(ch_qual_flag.data[2, 3, 4]))
        self.assertTrue(np.isnan(ch_qual_flag.attrs['_FillValue']))
        self.assertEqual("TODO", temp_misc_hk.attrs["long_name"])

        cold_counts = ds.variables["cold_counts"]
        self.assertEqual((4, 4, 28), cold_counts.shape)
        self.assertTrue(np.isnan(cold_counts.data[3, 0, 5]))
        self.assertTrue(np.isnan(cold_counts.attrs['_FillValue']))
        self.assertEqual("TODO", cold_counts.attrs["long_name"])

        counts_to_tb_gain = ds.variables["counts_to_tb_gain"]
        self.assertEqual((5, 4), counts_to_tb_gain.shape)
        self.assertTrue(np.isnan(counts_to_tb_gain.data[2, 3]))
        self.assertTrue(np.isnan(counts_to_tb_gain.attrs['_FillValue']))
        self.assertEqual("TODO", counts_to_tb_gain.attrs["long_name"])

        counts_to_tb_offset = ds.variables["counts_to_tb_offset"]
        self.assertEqual((5, 4), counts_to_tb_offset.shape)
        self.assertTrue(np.isnan(counts_to_tb_offset.data[3, 0]))
        self.assertTrue(np.isnan(counts_to_tb_offset.attrs['_FillValue']))
        self.assertEqual("TODO", counts_to_tb_offset.attrs["long_name"])

        gain_control = ds.variables["gain_control"]
        self.assertEqual((5, 4), gain_control.shape)
        self.assertTrue(np.isnan(gain_control.data[4, 1]))
        self.assertTrue(np.isnan(gain_control.attrs['_FillValue']))
        self.assertEqual("TODO", gain_control.attrs["long_name"])

        tb = ds.variables["tb"]
        self.assertEqual((5, 4, 28), tb.shape)
        self.assertTrue(np.isnan(tb.data[2, 3, 4]))
        self.assertTrue(np.isnan(tb.attrs['_FillValue']))
        self.assertEqual("TODO", tb.attrs["long_name"])
        self.assertEqual("toa_brightness_temperature", tb.attrs["standard_name"])
        self.assertEqual("K", tb.attrs["units"])

        thermal_reference = ds.variables["thermal_reference"]
        self.assertEqual((4,), thermal_reference.shape)
        self.assertTrue(np.isnan(thermal_reference.data[3]))
        self.assertTrue(np.isnan(thermal_reference.attrs['_FillValue']))
        self.assertEqual("TODO", thermal_reference.attrs["long_name"])
        self.assertEqual("TODO", thermal_reference.attrs["units"])

        warm_counts = ds.variables["warm_counts"]
        self.assertEqual((4, 4, 28), warm_counts.shape)
        self.assertTrue(np.isnan(warm_counts.data[0, 1, 6]))
        self.assertTrue(np.isnan(warm_counts.attrs['_FillValue']))
        self.assertEqual("TODO", warm_counts.attrs["long_name"])

    def test_get_swath_width(self):
        self.assertEqual(28, SSMT2.get_swath_width())

    def test_add_easy_fcdr_variables(self):
        ds = xr.Dataset()
        SSMT2.add_easy_fcdr_variables(ds, 4)

        u_independent_tb = ds.variables["u_independent_tb"]
        self.assertEqual((5, 4, 28), u_independent_tb.shape)
        self.assertTrue(np.isnan(u_independent_tb.data[0, 3, 0]))
        self.assertTrue(np.isnan(u_independent_tb.attrs["_FillValue"]))
        self.assertEqual("independent uncertainty per pixel", u_independent_tb.attrs["long_name"])
        self.assertEqual("K", u_independent_tb.attrs["units"])

        u_structured_tb = ds.variables["u_structured_tb"]
        self.assertEqual((5, 4, 28), u_structured_tb.shape)
        self.assertTrue(np.isnan(u_structured_tb.data[1, 0, 1]))
        self.assertTrue(np.isnan(u_structured_tb.attrs["_FillValue"]))
        self.assertEqual("structured uncertainty per pixel", u_structured_tb.attrs["long_name"])
        self.assertEqual("K", u_structured_tb.attrs["units"])

    def test_add_full_fcdr_variables(self):
        ds = xr.Dataset()
        SSMT2.add_full_fcdr_variables(ds, 5)

        u_temp_misc = ds.variables["u_Temperature_misc_housekeeping"]
        self.assertEqual((18, 5), u_temp_misc.shape)
        self.assertTrue(np.isnan(u_temp_misc.data[3, 2]))
        self.assertTrue(np.isnan(u_temp_misc.attrs["_FillValue"]))
        self.assertEqual("TODO", u_temp_misc.attrs["long_name"])
        self.assertEqual("TODO", u_temp_misc.attrs["units"])

        u_cold_counts = ds.variables["u_cold_counts"]
        self.assertEqual((4, 5, 28), u_cold_counts.shape)
        self.assertTrue(np.isnan(u_cold_counts.data[1, 1, 6]))
        self.assertTrue(np.isnan(u_cold_counts.attrs['_FillValue']))
        self.assertEqual("TODO", u_cold_counts.attrs["long_name"])

        u_counts_to_tb_gain = ds.variables["u_counts_to_tb_gain"]
        self.assertEqual((5, 5), u_counts_to_tb_gain.shape)
        self.assertTrue(np.isnan(u_counts_to_tb_gain.data[2, 3]))
        self.assertTrue(np.isnan(u_counts_to_tb_gain.attrs['_FillValue']))
        self.assertEqual("TODO", u_counts_to_tb_gain.attrs["long_name"])

        u_counts_to_tb_offset = ds.variables["u_counts_to_tb_offset"]
        self.assertEqual((5, 5), u_counts_to_tb_offset.shape)
        self.assertTrue(np.isnan(u_counts_to_tb_offset.data[3, 3]))
        self.assertTrue(np.isnan(u_counts_to_tb_offset.attrs['_FillValue']))
        self.assertEqual("TODO", u_counts_to_tb_offset.attrs["long_name"])

        u_gain_control = ds.variables["u_gain_control"]
        self.assertEqual((5, 5), u_gain_control.shape)
        self.assertTrue(np.isnan(u_gain_control.data[4, 3]))
        self.assertTrue(np.isnan(u_gain_control.attrs['_FillValue']))
        self.assertEqual("TODO", u_gain_control.attrs["long_name"])

        u_tb = ds.variables["u_tb"]
        self.assertEqual((5, 5, 28), u_tb.shape)
        self.assertTrue(np.isnan(u_tb.data[3, 4, 5]))
        self.assertTrue(np.isnan(u_tb.attrs['_FillValue']))
        self.assertEqual("TODO", u_tb.attrs["long_name"])
        self.assertEqual("K", u_tb.attrs["units"])

        u_thermal_reference = ds.variables["u_thermal_reference"]
        self.assertEqual((5,), u_thermal_reference.shape)
        self.assertTrue(np.isnan(u_thermal_reference.data[3]))
        self.assertTrue(np.isnan(u_thermal_reference.attrs['_FillValue']))
        self.assertEqual("TODO", u_thermal_reference.attrs["long_name"])
        self.assertEqual("TODO", u_thermal_reference.attrs["units"])

        u_warm_counts = ds.variables["u_warm_counts"]
        self.assertEqual((4, 5, 28), u_warm_counts.shape)
        self.assertTrue(np.isnan(u_warm_counts.data[2, 3, 7]))
        self.assertTrue(np.isnan(u_warm_counts.attrs['_FillValue']))
        self.assertEqual("TODO", u_warm_counts.attrs["long_name"])

    def test_add_template_key(self):
        ds = xr.Dataset()

        SSMT2.add_template_key(ds)

        self.assertEqual("SSMT2", ds.attrs["template_key"])