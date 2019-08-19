import unittest

import numpy as np
import xarray as xr

from fiduceo.cdr.writer.templates.uth import UTH
from fiduceo.common.test.assertions import Assertions


class UTHTest(unittest.TestCase):

    def test_add_variables(self):
        ds = xr.Dataset()

        UTH.add_variables(ds, 360, 100)

        Assertions.assert_gridded_geolocation_variables(self, ds, 360, 100)
        Assertions.assert_quality_flags(self, ds, 360, 100)

        time_ranges_asc = ds.variables["time_ranges_ascend"]
        self.assertEqual((2, 100, 360), time_ranges_asc.shape)
        self.assertEqual(-1, time_ranges_asc.values[0, 66, 178])
        self.assertEqual(4294967295, time_ranges_asc.attrs["_FillValue"])
        self.assertEqual("Minimum and maximum seconds of day pixel contribution time, ascending nodes", time_ranges_asc.attrs["description"])
        self.assertEqual("s", time_ranges_asc.attrs["units"])
        self.assertEqual("lon lat", time_ranges_asc.attrs["coordinates"])

        time_ranges_desc = ds.variables["time_ranges_descend"]
        self.assertEqual((2, 100, 360), time_ranges_desc.shape)
        self.assertEqual(-1, time_ranges_desc.values[1, 67, 179])
        self.assertEqual(4294967295, time_ranges_desc.attrs["_FillValue"])
        self.assertEqual("Minimum and maximum seconds of day pixel contribution time, descending nodes", time_ranges_desc.attrs["description"])
        self.assertEqual("s", time_ranges_desc.attrs["units"])
        self.assertEqual("lon lat", time_ranges_desc.attrs["coordinates"])

        obs_count_asc = ds.variables["observation_count_ascend"]
        self.assertEqual((100, 360), obs_count_asc.shape)
        self.assertEqual(-32767, obs_count_asc.values[67, 179])
        self.assertEqual(-32767, obs_count_asc.attrs["_FillValue"])
        self.assertEqual("Number of UTH/brightness temperature observations in a grid box for ascending passes", obs_count_asc.attrs["description"])
        self.assertEqual("lon lat", obs_count_asc.attrs["coordinates"])

        obs_count_desc = ds.variables["observation_count_descend"]
        self.assertEqual((100, 360), obs_count_desc.shape)
        self.assertEqual(-32767, obs_count_desc.values[68, 180])
        self.assertEqual(-32767, obs_count_desc.attrs["_FillValue"])
        self.assertEqual("Number of UTH/brightness temperature observations in a grid box for descending passes", obs_count_desc.attrs["description"])
        self.assertEqual("lon lat", obs_count_desc.attrs["coordinates"])

        overp_count_asc = ds.variables["overpass_count_ascend"]
        self.assertEqual((100, 360), overp_count_asc.shape)
        self.assertEqual(255, overp_count_asc.values[69, 179])
        self.assertEqual(255, overp_count_asc.attrs["_FillValue"])
        self.assertEqual("Number of satellite overpasses in a grid box for ascending passes", overp_count_asc.attrs["description"])
        self.assertEqual("lon lat", overp_count_asc.attrs["coordinates"])

        overp_count_desc = ds.variables["overpass_count_descend"]
        self.assertEqual((100, 360), overp_count_desc.shape)
        self.assertEqual(255, overp_count_desc.values[70, 180])
        self.assertEqual(255, overp_count_desc.attrs["_FillValue"])
        self.assertEqual("Number of satellite overpasses in a grid box for descending passes", overp_count_desc.attrs["description"])
        self.assertEqual("lon lat", overp_count_desc.attrs["coordinates"])

        uth_asc = ds.variables["uth_ascend"]
        self.assertEqual((100, 360), uth_asc.shape)
        self.assertTrue(np.isnan(uth_asc.values[97, 198]))
        self.assertTrue(np.isnan(uth_asc.attrs["_FillValue"]))
        self.assertEqual("lon lat", uth_asc.attrs["coordinates"])
        self.assertEqual("%", uth_asc.attrs["units"])
        self.assertEqual("upper_tropospheric_humidity", uth_asc.attrs["long_name"])
        self.assertEqual("Monthly average of all UTH retrievals in a grid box for ascending passes (calculated from daily averages)", uth_asc.attrs["description"])

        uth_desc = ds.variables["uth_descend"]
        self.assertEqual((100, 360), uth_desc.shape)
        self.assertTrue(np.isnan(uth_desc.values[98, 199]))
        self.assertTrue(np.isnan(uth_desc.attrs["_FillValue"]))
        self.assertEqual("lon lat", uth_desc.attrs["coordinates"])
        self.assertEqual("%", uth_desc.attrs["units"])
        self.assertEqual("upper_tropospheric_humidity", uth_desc.attrs["long_name"])
        self.assertEqual("Monthly average of all UTH retrievals in a grid box for descending passes (calculated from daily averages)", uth_desc.attrs["description"])

        u_ind_uth_asc = ds.variables["u_independent_uth_ascend"]
        self.assertEqual((100, 360), u_ind_uth_asc.shape)
        self.assertTrue(np.isnan(u_ind_uth_asc.values[98, 199]))
        self.assertTrue(np.isnan(u_ind_uth_asc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of UTH due to independent effects for ascending passes", u_ind_uth_asc.attrs["description"])
        self.assertEqual("lon lat", u_ind_uth_asc.attrs["coordinates"])
        self.assertEqual("%", u_ind_uth_asc.attrs["units"])

        u_ind_uth_desc = ds.variables["u_independent_uth_descend"]
        self.assertEqual((100, 360), u_ind_uth_desc.shape)
        self.assertTrue(np.isnan(u_ind_uth_desc.values[99, 200]))
        self.assertTrue(np.isnan(u_ind_uth_desc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of UTH due to independent effects for descending passes", u_ind_uth_desc.attrs["description"])
        self.assertEqual("lon lat", u_ind_uth_desc.attrs["coordinates"])
        self.assertEqual("%", u_ind_uth_desc.attrs["units"])

        u_str_uth_asc = ds.variables["u_structured_uth_ascend"]
        self.assertEqual((100, 360), u_str_uth_asc.shape)
        self.assertTrue(np.isnan(u_str_uth_asc.values[0, 201]))
        self.assertTrue(np.isnan(u_str_uth_asc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of UTH due to structured effects for ascending passes", u_str_uth_asc.attrs["description"])
        self.assertEqual("lon lat", u_str_uth_asc.attrs["coordinates"])
        self.assertEqual("%", u_str_uth_asc.attrs["units"])

        u_str_uth_desc = ds.variables["u_structured_uth_descend"]
        self.assertEqual((100, 360), u_str_uth_desc.shape)
        self.assertTrue(np.isnan(u_str_uth_desc.values[1, 202]))
        self.assertTrue(np.isnan(u_str_uth_desc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of UTH due to structured effects for descending passes", u_str_uth_desc.attrs["description"])
        self.assertEqual("lon lat", u_str_uth_desc.attrs["coordinates"])
        self.assertEqual("%", u_str_uth_desc.attrs["units"])

        u_com_uth_asc = ds.variables["u_common_uth_ascend"]
        self.assertEqual((100, 360), u_com_uth_asc.shape)
        self.assertTrue(np.isnan(u_com_uth_asc.values[70, 171]))
        self.assertTrue(np.isnan(u_com_uth_asc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of UTH due to common effects for ascending passes", u_com_uth_asc.attrs["description"])
        self.assertEqual("lon lat", u_com_uth_asc.attrs["coordinates"])
        self.assertEqual("%", u_com_uth_asc.attrs["units"])

        u_com_uth_desc = ds.variables["u_common_uth_descend"]
        self.assertEqual((100, 360), u_com_uth_desc.shape)
        self.assertTrue(np.isnan(u_com_uth_desc.values[71, 172]))
        self.assertTrue(np.isnan(u_com_uth_desc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of UTH due to common effects for descending passes", u_com_uth_desc.attrs["description"])
        self.assertEqual("lon lat", u_com_uth_desc.attrs["coordinates"])
        self.assertEqual("%", u_com_uth_desc.attrs["units"])

        uth_inhom_asc = ds.variables["uth_inhomogeneity_ascend"]
        self.assertEqual((100, 360), uth_inhom_asc.shape)
        self.assertTrue(np.isnan(uth_inhom_asc.values[71, 172]))
        self.assertTrue(np.isnan(uth_inhom_asc.attrs["_FillValue"]))
        self.assertEqual("Standard deviation of all daily UTH averages which were used to calculate the monthly UTH average in a grid box for ascending passes", uth_inhom_asc.attrs["description"])
        self.assertEqual("lon lat", uth_inhom_asc.attrs["coordinates"])
        self.assertEqual("%", uth_inhom_asc.attrs["units"])

        uth_inhom_desc = ds.variables["uth_inhomogeneity_descend"]
        self.assertEqual((100, 360), uth_inhom_desc.shape)
        self.assertTrue(np.isnan(uth_inhom_desc.values[72, 173]))
        self.assertTrue(np.isnan(uth_inhom_desc.attrs["_FillValue"]))
        self.assertEqual("Standard deviation of all daily UTH averages which were used to calculate the monthly UTH average in a grid box for descending passes", uth_inhom_desc.attrs["description"])
        self.assertEqual("lon lat", uth_inhom_desc.attrs["coordinates"])
        self.assertEqual("%", uth_inhom_desc.attrs["units"])

        bt_ascend = ds.variables["BT_ascend"]
        self.assertEqual((100, 360), bt_ascend.shape)
        self.assertTrue(np.isnan(bt_ascend.values[73, 174]))
        self.assertTrue(np.isnan(bt_ascend.attrs["_FillValue"]))
        self.assertEqual("Monthly average of all brightness temperatures which were used to retrieve UTH in a grid box for ascending passes (calculated from daily averages)", bt_ascend.attrs["description"])
        self.assertEqual("lon lat", bt_ascend.attrs["coordinates"])
        self.assertEqual("K", bt_ascend.attrs["units"])
        self.assertEqual("toa_brightness_temperature", bt_ascend.attrs["standard_name"])

        bt_descend = ds.variables["BT_descend"]
        self.assertEqual((100, 360), bt_descend.shape)
        self.assertTrue(np.isnan(bt_descend.values[74, 175]))
        self.assertTrue(np.isnan(bt_descend.attrs["_FillValue"]))
        self.assertEqual("Monthly average of all brightness temperatures which were used to retrieve UTH in a grid box for descending passes (calculated from daily averages)", bt_descend.attrs["description"])
        self.assertEqual("lon lat", bt_descend.attrs["coordinates"])
        self.assertEqual("K", bt_descend.attrs["units"])
        self.assertEqual("toa_brightness_temperature", bt_descend.attrs["standard_name"])

        u_ind_bt_asc = ds.variables["u_independent_BT_ascend"]
        self.assertEqual((100, 360), u_ind_bt_asc.shape)
        self.assertTrue(np.isnan(u_ind_bt_asc.values[99, 200]))
        self.assertTrue(np.isnan(u_ind_bt_asc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of brightness temperature due to independent effects for ascending passes", u_ind_bt_asc.attrs["description"])
        self.assertEqual("lon lat", u_ind_bt_asc.attrs["coordinates"])
        self.assertEqual("K", u_ind_bt_asc.attrs["units"])

        u_ind_bt_desc = ds.variables["u_independent_BT_descend"]
        self.assertEqual((100, 360), u_ind_bt_desc.shape)
        self.assertTrue(np.isnan(u_ind_bt_desc.values[0, 201]))
        self.assertTrue(np.isnan(u_ind_bt_desc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of brightness temperature due to independent effects for descending passes", u_ind_bt_desc.attrs["description"])
        self.assertEqual("lon lat", u_ind_bt_desc.attrs["coordinates"])
        self.assertEqual("K", u_ind_bt_desc.attrs["units"])

        u_str_bt_asc = ds.variables["u_structured_BT_ascend"]
        self.assertEqual((100, 360), u_str_bt_asc.shape)
        self.assertTrue(np.isnan(u_str_bt_asc.values[1, 202]))
        self.assertTrue(np.isnan(u_str_bt_asc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of brightness temperature due to structured effects for ascending passes", u_str_bt_asc.attrs["description"])
        self.assertEqual("lon lat", u_str_bt_asc.attrs["coordinates"])
        self.assertEqual("K", u_str_bt_asc.attrs["units"])

        u_str_bt_desc = ds.variables["u_structured_BT_descend"]
        self.assertEqual((100, 360), u_str_bt_desc.shape)
        self.assertTrue(np.isnan(u_str_bt_desc.values[2, 203]))
        self.assertTrue(np.isnan(u_str_bt_desc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of brightness temperature due to structured effects for descending passes", u_str_bt_desc.attrs["description"])
        self.assertEqual("lon lat", u_str_bt_desc.attrs["coordinates"])
        self.assertEqual("K", u_str_bt_desc.attrs["units"])

        u_com_bt_asc = ds.variables["u_common_BT_ascend"]
        self.assertEqual((100, 360), u_com_bt_asc.shape)
        self.assertTrue(np.isnan(u_com_bt_asc.values[3, 204]))
        self.assertTrue(np.isnan(u_str_bt_asc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of brightness temperature due to common effects for ascending passes", u_com_bt_asc.attrs["description"])
        self.assertEqual("lon lat", u_com_bt_asc.attrs["coordinates"])
        self.assertEqual("K", u_com_bt_asc.attrs["units"])

        u_com_bt_desc = ds.variables["u_common_BT_descend"]
        self.assertEqual((100, 360), u_com_bt_desc.shape)
        self.assertTrue(np.isnan(u_com_bt_desc.values[5, 206]))
        self.assertTrue(np.isnan(u_com_bt_desc.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of brightness temperature due to common effects for descending passes", u_com_bt_desc.attrs["description"])
        self.assertEqual("lon lat", u_com_bt_desc.attrs["coordinates"])
        self.assertEqual("K", u_com_bt_desc.attrs["units"])

        bt_inhom_asc = ds.variables["BT_inhomogeneity_ascend"]
        self.assertEqual((100, 360), bt_inhom_asc.shape)
        self.assertTrue(np.isnan(bt_inhom_asc.values[72, 173]))
        self.assertTrue(np.isnan(bt_inhom_asc.attrs["_FillValue"]))
        self.assertEqual("Standard deviation of all daily brightness temperature averages which were used to calculate the monthly brightness temperature average for ascending passes", bt_inhom_asc.attrs["description"])
        self.assertEqual("lon lat", bt_inhom_asc.attrs["coordinates"])
        self.assertEqual("K", bt_inhom_asc.attrs["units"])

        bt_inhom_desc = ds.variables["BT_inhomogeneity_descend"]
        self.assertEqual((100, 360), bt_inhom_desc.shape)
        self.assertTrue(np.isnan(bt_inhom_desc.values[73, 174]))
        self.assertTrue(np.isnan(bt_inhom_desc.attrs["_FillValue"]))
        self.assertEqual("Standard deviation of all daily brightness temperature averages which were used to calculate the monthly brightness temperature average for descending passes", bt_inhom_desc.attrs["description"])
        self.assertEqual("lon lat", bt_inhom_desc.attrs["coordinates"])
        self.assertEqual("K", bt_inhom_desc.attrs["units"])

        obs_count_all_asc = ds.variables["observation_count_all_ascend"]
        self.assertEqual((100, 360), obs_count_all_asc.shape)
        self.assertEqual(-32767, obs_count_all_asc.values[68, 180])
        self.assertEqual(-32767, obs_count_all_asc.attrs["_FillValue"])
        self.assertEqual("Number of all observations in a grid box for ascending passes - no filtering done", obs_count_all_asc.attrs["description"])
        self.assertEqual("lon lat", obs_count_all_asc.attrs["coordinates"])

        obs_count_all_desc = ds.variables["observation_count_all_descend"]
        self.assertEqual((100, 360), obs_count_all_desc.shape)
        self.assertEqual(-32767, obs_count_all_desc.values[69, 181])
        self.assertEqual(-32767, obs_count_all_desc.attrs["_FillValue"])
        self.assertEqual("Number of all observations in a grid box for descending passes - no filtering done", obs_count_all_desc.attrs["description"])
        self.assertEqual("lon lat", obs_count_all_desc.attrs["coordinates"])