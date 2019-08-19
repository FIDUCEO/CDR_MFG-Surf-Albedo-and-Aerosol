import unittest

import numpy as np
import xarray as xr

from fiduceo.common.test.assertions import Assertions
from fiduceo.common.writer.default_data import DefaultData
from fiduceo.fcdr.writer.templates.amsub_mhs import AMSUB_MHS


class AMSUB_MHSTest(unittest.TestCase):
    def test_add_original_variables(self):
        ds = xr.Dataset()
        AMSUB_MHS.add_original_variables(ds, 4)

        Assertions.assert_geolocation_variables(self, ds, 90, 4)
        Assertions.assert_quality_flags(self, ds, 90, 4)

        btemps = ds.variables["btemps"]
        self.assertEqual((5, 4, 90), btemps.shape)
        self.assertTrue(np.isnan(btemps.data[0, 2, 0]))
        self.assertEqual("toa_brightness_temperature", btemps.attrs["standard_name"])
        self.assertEqual("K", btemps.attrs["units"])
        self.assertEqual("chanqual qualind scanqual", btemps.attrs["ancillary_variables"])
        self.assertEqual(np.int32, btemps.encoding['dtype'])
        self.assertEqual(-999999, btemps.encoding['_FillValue'])
        self.assertEqual(0.01, btemps.encoding['scale_factor'])
        self.assertEqual(0.0, btemps.encoding['add_offset'])

        chanqual = ds.variables["chanqual"]
        self.assertEqual((5, 4), chanqual.shape)
        self.assertEqual(0, chanqual.data[0, 3])
        self.assertEqual("status_flag", chanqual.attrs["standard_name"])
        self.assertEqual("1, 2, 4, 8, 16, 32", chanqual.attrs["flag_masks"])
        self.assertEqual("some_bad_prt_temps some_bad_space_view_counts some_bad_bb_counts no_good_prt_temps no_good_space_view_counts no_good_bb_counts", chanqual.attrs["flag_meanings"])

        instrtemp = ds.variables["instrtemp"]
        self.assertEqual((4,), instrtemp.shape)
        self.assertTrue(np.isnan(instrtemp.data[0]))
        self.assertEqual("K", instrtemp.attrs["units"])
        self.assertEqual("instrument_temperature", instrtemp.attrs["long_name"])
        self.assertEqual(np.int32, instrtemp.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), instrtemp.encoding['_FillValue'])
        self.assertEqual(0.01, instrtemp.encoding['scale_factor'])
        self.assertEqual(0.0, instrtemp.encoding['add_offset'])

        qualind = ds.variables["qualind"]
        self.assertEqual((4,), qualind.shape)
        self.assertEqual(0, qualind.data[1])
        self.assertEqual("status_flag", qualind.attrs["standard_name"])
        self.assertEqual("33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648", qualind.attrs["flag_masks"])
        self.assertEqual("instr_status_changed first_good_clock_update no_earth_loc no_calib data_gap_precedes time_seq_error not_use_scan", qualind.attrs["flag_meanings"])

        scanqual = ds.variables["scanqual"]
        self.assertEqual((4,), scanqual.shape)
        self.assertEqual(0, scanqual.data[1])
        self.assertEqual("status_flag", scanqual.attrs["standard_name"])
        self.assertEqual("8, 16, 32, 64, 128, 1024, 2048, 4096, 8192, 16384, 32768, 1048576, 2097152, 4194304, 8388608", scanqual.attrs["flag_masks"])
        self.assertEqual(
            "earth_loc_quest_ant_pos earth_loc_quest_reas earth_loc_quest_margin earth_loc_quest_time no_earth_loc_time uncalib_instr_mode uncalib_channels calib_marg_prt uncalib_bad_prt calib_few_scans uncalib_bad_time repeat_scan_times inconsistent_time time_field_bad time_field_inferred",
            scanqual.attrs["flag_meanings"])

        scnlin = ds.variables["scnlin"]
        self.assertEqual((4,), scnlin.shape)
        self.assertEqual(-2147483647, scnlin.data[2])
        self.assertEqual(-2147483647, scnlin.attrs["_FillValue"])
        self.assertEqual("scanline", scnlin.attrs["long_name"])

        scnlindy = ds.variables["scnlindy"]
        self.assertEqual((4,), scnlindy.shape)
        self.assertEqual(-2147483647, scnlindy.data[3])
        self.assertEqual(-2147483647, scnlindy.attrs["_FillValue"])
        self.assertEqual("Acquisition day of year of scan", scnlindy.attrs["long_name"])

        scnlintime = ds.variables["scnlintime"]
        self.assertEqual((4,), scnlintime.shape)
        self.assertEqual(-2147483647, scnlintime.data[0])
        self.assertEqual(-2147483647, scnlintime.attrs["_FillValue"])
        self.assertEqual("Acquisition time of scan in milliseconds since beginning of the day", scnlintime.attrs["long_name"])
        self.assertEqual("ms", scnlintime.attrs["units"])

        scnlinyr = ds.variables["scnlinyr"]
        self.assertEqual((4,), scnlinyr.shape)
        self.assertEqual(-2147483647, scnlinyr.data[1])
        self.assertEqual(-2147483647, scnlinyr.attrs["_FillValue"])
        self.assertEqual("Acquisition year of scan", scnlinyr.attrs["long_name"])

        sat_azimuth = ds.variables["satellite_azimuth_angle"]
        self.assertEqual((4, 90), sat_azimuth.shape)
        self.assertTrue(np.isnan(sat_azimuth.data[2, 1]))
        self.assertEqual("sensor_azimuth_angle", sat_azimuth.attrs["standard_name"])
        self.assertEqual("degree", sat_azimuth.attrs["units"])
        self.assertEqual(np.int32, sat_azimuth.encoding['dtype'])
        self.assertEqual(-999999, sat_azimuth.encoding['_FillValue'])
        self.assertEqual(0.01, sat_azimuth.encoding['scale_factor'])
        self.assertEqual(0.0, sat_azimuth.encoding['add_offset'])

        sat_zenith = ds.variables["satellite_zenith_angle"]
        self.assertEqual((4, 90), sat_zenith.shape)
        self.assertTrue(np.isnan(sat_zenith.data[2, 1]))
        self.assertEqual("sensor_zenith_angle", sat_zenith.attrs["standard_name"])
        self.assertEqual("degree", sat_zenith.attrs["units"])
        self.assertEqual(np.int32, sat_zenith.encoding['dtype'])
        self.assertEqual(-999999, sat_zenith.encoding['_FillValue'])
        self.assertEqual(0.01, sat_zenith.encoding['scale_factor'])
        self.assertEqual(0.0, sat_zenith.encoding['add_offset'])

        sol_azimuth = ds.variables["solar_azimuth_angle"]
        self.assertEqual((4, 90), sol_azimuth.shape)
        self.assertTrue(np.isnan(sol_azimuth.data[3, 0]))
        self.assertEqual("solar_azimuth_angle", sol_azimuth.attrs["standard_name"])
        self.assertEqual("degree", sol_azimuth.attrs["units"])
        self.assertEqual(np.int32, sol_azimuth.encoding['dtype'])
        self.assertEqual(-999999, sol_azimuth.encoding['_FillValue'])
        self.assertEqual(0.01, sol_azimuth.encoding['scale_factor'])
        self.assertEqual(0.0, sol_azimuth.encoding['add_offset'])

        sol_zenith = ds.variables["solar_zenith_angle"]
        self.assertEqual((4, 90), sol_zenith.shape)
        self.assertTrue(np.isnan(sol_zenith.data[3, 0]))
        self.assertEqual("solar_zenith_angle", sol_zenith.attrs["standard_name"])
        self.assertEqual("degree", sol_zenith.attrs["units"])
        self.assertEqual(np.int32, sol_zenith.encoding['dtype'])
        self.assertEqual(-999999, sol_zenith.encoding['_FillValue'])
        self.assertEqual(0.01, sol_zenith.encoding['scale_factor'])
        self.assertEqual(0.0, sol_zenith.encoding['add_offset'])

        acquisition_time = ds.variables["acquisition_time"]
        self.assertEqual((4,), acquisition_time.shape)
        self.assertEqual(-2147483647, acquisition_time.data[2])
        self.assertEqual(-2147483647, acquisition_time.attrs["_FillValue"])
        self.assertEqual("time", acquisition_time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", acquisition_time.attrs["long_name"])
        self.assertEqual("s", acquisition_time.attrs["units"])

    def test_get_swath_width(self):
        self.assertEqual(90, AMSUB_MHS.get_swath_width())

    def test_add_easy_fcdr_variables(self):
        ds = xr.Dataset()
        AMSUB_MHS.add_easy_fcdr_variables(ds, 4)

        u_ind_btemps = ds.variables["u_independent_btemps"]
        self.assertEqual((5, 4, 90), u_ind_btemps.shape)
        self.assertTrue(np.isnan(u_ind_btemps.data[4, 2, 35]))
        self.assertTrue(np.isnan(u_ind_btemps.attrs["_FillValue"]))
        self.assertEqual("independent uncertainty per pixel", u_ind_btemps.attrs["long_name"])
        self.assertEqual("K", u_ind_btemps.attrs["units"])

        u_struct_btemps = ds.variables["u_structured_btemps"]
        self.assertEqual((5, 4, 90), u_struct_btemps.shape)
        self.assertTrue(np.isnan(u_struct_btemps.data[0, 3, 36]))
        self.assertTrue(np.isnan(u_struct_btemps.attrs["_FillValue"]))
        self.assertEqual("structured uncertainty per pixel", u_struct_btemps.attrs["long_name"])
        self.assertEqual("K", u_struct_btemps.attrs["units"])

    def test_add_full_fcdr_variables(self):
        ds = xr.Dataset()
        AMSUB_MHS.add_full_fcdr_variables(ds, 4)

        u_btemps = ds.variables["u_btemps"]
        self.assertEqual((5, 4, 90), u_btemps.shape)
        self.assertTrue(np.isnan(u_btemps.data[3, 1, 34]))
        self.assertTrue(np.isnan(u_btemps.attrs["_FillValue"]))
        self.assertEqual("total uncertainty of brightness temperature", u_btemps.attrs["long_name"])
        self.assertEqual("K", u_btemps.attrs["units"])

        u_syst_btemps = ds.variables["u_syst_btemps"]
        self.assertEqual((5, 4, 90), u_syst_btemps.shape)
        self.assertTrue(np.isnan(u_syst_btemps.data[4, 2, 35]))
        self.assertTrue(np.isnan(u_syst_btemps.attrs["_FillValue"]))
        self.assertEqual("systematic uncertainty of brightness temperature", u_syst_btemps.attrs["long_name"])
        self.assertEqual("K", u_syst_btemps.attrs["units"])

        u_random_btemps = ds.variables["u_random_btemps"]
        self.assertEqual((5, 4, 90), u_random_btemps.shape)
        self.assertTrue(np.isnan(u_random_btemps.data[0, 3, 36]))
        self.assertTrue(np.isnan(u_random_btemps.attrs["_FillValue"]))
        self.assertEqual("noise on brightness temperature", u_random_btemps.attrs["long_name"])
        self.assertEqual("K", u_random_btemps.attrs["units"])

        u_instrtemp = ds.variables["u_instrtemp"]
        self.assertEqual((4,), u_instrtemp.shape)
        self.assertTrue(np.isnan(u_instrtemp.data[1]))
        self.assertTrue(np.isnan(u_instrtemp.attrs["_FillValue"]))
        self.assertEqual("uncertainty of instrument temperature", u_instrtemp.attrs["long_name"])
        self.assertEqual("K", u_instrtemp.attrs["units"])

        u_latitude = ds.variables["u_latitude"]
        self.assertEqual((4, 90), u_latitude.shape)
        self.assertTrue(np.isnan(u_latitude.data[2, 38]))
        self.assertTrue(np.isnan(u_latitude.attrs["_FillValue"]))
        self.assertEqual("uncertainty of latitude", u_latitude.attrs["long_name"])
        self.assertEqual("degree", u_latitude.attrs["units"])

        u_longitude = ds.variables["u_longitude"]
        self.assertEqual((4, 90), u_longitude.shape)
        self.assertTrue(np.isnan(u_longitude.data[3, 39]))
        self.assertTrue(np.isnan(u_longitude.attrs["_FillValue"]))
        self.assertEqual("uncertainty of longitude", u_longitude.attrs["long_name"])
        self.assertEqual("degree", u_longitude.attrs["units"])

        u_sat_azimuth = ds.variables["u_satellite_azimuth_angle"]
        self.assertEqual((4, 90), u_sat_azimuth.shape)
        self.assertTrue(np.isnan(u_sat_azimuth.data[0, 40]))
        self.assertTrue(np.isnan(u_sat_azimuth.attrs["_FillValue"]))
        self.assertEqual("uncertainty of satellite azimuth angle", u_sat_azimuth.attrs["long_name"])
        self.assertEqual("degree", u_sat_azimuth.attrs["units"])

        u_sat_zenith = ds.variables["u_satellite_zenith_angle"]
        self.assertEqual((4, 90), u_sat_zenith.shape)
        self.assertTrue(np.isnan(u_sat_zenith.data[1, 41]))
        self.assertTrue(np.isnan(u_sat_zenith.attrs["_FillValue"]))
        self.assertEqual("uncertainty of satellite zenith angle", u_sat_zenith.attrs["long_name"])
        self.assertEqual("degree", u_sat_zenith.attrs["units"])

        u_sol_azimuth = ds.variables["u_solar_azimuth_angle"]
        self.assertEqual((4, 90), u_sol_azimuth.shape)
        self.assertTrue(np.isnan(u_sol_azimuth.data[2, 42]))
        self.assertTrue(np.isnan(u_sol_azimuth.attrs["_FillValue"]))
        self.assertEqual("uncertainty of solar azimuth angle", u_sol_azimuth.attrs["long_name"])
        self.assertEqual("degree", u_sol_azimuth.attrs["units"])

        u_sol_zenith = ds.variables["u_solar_zenith_angle"]
        self.assertEqual((4, 90), u_sol_zenith.shape)
        self.assertTrue(np.isnan(u_sol_zenith.data[3, 43]))
        self.assertTrue(np.isnan(u_sol_zenith.attrs["_FillValue"]))
        self.assertEqual("uncertainty of solar zenith angle", u_sol_zenith.attrs["long_name"])
        self.assertEqual("degree", u_sol_zenith.attrs["units"])

    def test_add_template_key(self):
        ds = xr.Dataset()

        AMSUB_MHS.add_template_key(ds)

        self.assertEqual("AMSUB_MHS", ds.attrs["template_key"])
