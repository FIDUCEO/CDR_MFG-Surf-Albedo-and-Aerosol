import unittest

import numpy as np

from fiduceo.common.writer.default_data import DefaultData

CHUNKING_2D = (512, 56)


class HIRSAssert(unittest.TestCase):
    # this is a mean hack - there is some error in the framework when using Python 2.7 that requests this method tb 2017-05-10
    def runTest(self):
        pass

    def assert_bt_variable(self, ds, chunking=None):
        bt = ds.variables["bt"]
        self.assertEqual((19, 6, 56), bt.shape)
        self.assertTrue(np.isnan(bt.data[0, 2, 1]))
        self.assertEqual("toa_brightness_temperature", bt.attrs["standard_name"])
        self.assertEqual("Brightness temperature, NOAA/EUMETSAT calibrated", bt.attrs["long_name"])
        self.assertEqual("K", bt.attrs["units"])
        self.assertEqual("longitude latitude", bt.attrs["coordinates"])
        self.assertEqual(np.int16, bt.encoding['dtype'])
        self.assertEqual(-999, bt.encoding['_FillValue'])
        self.assertEqual(0.01, bt.encoding['scale_factor'])
        self.assertEqual(150.0, bt.encoding['add_offset'])
        if chunking is not None:
            self.assertEqual(chunking, bt.encoding['chunksizes'])

        self.assertEqual("quality_scanline_bitmask quality_channel_bitmask", bt.attrs["ancilliary_variables"])

    def assert_common_angles(self, ds, chunking=None):
        satellite_zenith_angle = ds.variables["satellite_zenith_angle"]
        self.assertEqual((6, 56), satellite_zenith_angle.shape)
        self.assertTrue(np.isnan(satellite_zenith_angle.data[2, 2]))
        self.assertEqual(np.uint16, satellite_zenith_angle.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), satellite_zenith_angle.encoding['_FillValue'])
        self.assertEqual(0.01, satellite_zenith_angle.encoding['scale_factor'])
        self.assertEqual(-180.0, satellite_zenith_angle.encoding['add_offset'])
        if chunking is not None:
            self.assertEqual(chunking, satellite_zenith_angle.encoding['chunksizes'])
        self.assertEqual("platform_zenith_angle", satellite_zenith_angle.attrs["standard_name"])
        self.assertEqual("degree", satellite_zenith_angle.attrs["units"])
        self.assertEqual("longitude latitude", satellite_zenith_angle.attrs["coordinates"])

        solar_zenith_angle = ds.variables["solar_zenith_angle"]
        self.assertEqual((6, 56), solar_zenith_angle.shape)
        self.assertTrue(np.isnan(solar_zenith_angle.data[3, 3]))
        self.assertEqual(np.uint16, solar_zenith_angle.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), solar_zenith_angle.encoding['_FillValue'])
        self.assertEqual(0.01, solar_zenith_angle.encoding['scale_factor'])
        self.assertEqual(-180.0, solar_zenith_angle.encoding['add_offset'])
        if chunking is not None:
            self.assertEqual(chunking, solar_zenith_angle.encoding['chunksizes'])
        self.assertEqual("solar_zenith_angle", solar_zenith_angle.attrs["standard_name"])
        self.assertEqual("solar_zenith_angle", solar_zenith_angle.attrs["orig_name"])
        self.assertEqual("degree", solar_zenith_angle.attrs["units"])
        self.assertEqual("longitude latitude", solar_zenith_angle.attrs["coordinates"])

        satellite_azimuth_angle = ds.variables["satellite_azimuth_angle"]
        self.assertEqual((6, 56), satellite_azimuth_angle.shape)
        self.assertTrue(np.isnan(satellite_azimuth_angle.data[5, 5]))
        self.assertEqual(np.uint16, satellite_azimuth_angle.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), satellite_azimuth_angle.encoding['_FillValue'])
        self.assertEqual(0.01, satellite_azimuth_angle.encoding['scale_factor'])
        self.assertEqual(-180.0, satellite_azimuth_angle.encoding['add_offset'])
        if chunking is not None:
            self.assertEqual(chunking, satellite_azimuth_angle.encoding['chunksizes'])
        self.assertEqual("sensor_azimuth_angle", satellite_azimuth_angle.attrs["standard_name"])
        self.assertEqual("local_azimuth_angle", satellite_azimuth_angle.attrs["long_name"])
        self.assertEqual("degree", satellite_azimuth_angle.attrs["units"])
        self.assertEqual("longitude latitude", satellite_azimuth_angle.attrs["coordinates"])

        solar_azimuth_angle = ds.variables["solar_azimuth_angle"]
        self.assertEqual((6, 56), solar_azimuth_angle.shape)
        self.assertTrue(np.isnan(solar_azimuth_angle.data[4, 4]))
        self.assertEqual(np.uint16, solar_azimuth_angle.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), solar_azimuth_angle.encoding['_FillValue'])
        self.assertEqual(0.01, solar_azimuth_angle.encoding['scale_factor'])
        self.assertEqual(-180.0, solar_azimuth_angle.encoding['add_offset'])
        if chunking is not None:
            self.assertEqual(chunking, solar_azimuth_angle.encoding['chunksizes'])
        self.assertEqual("solar_azimuth_angle", solar_azimuth_angle.attrs["standard_name"])
        self.assertEqual("degree", solar_azimuth_angle.attrs["units"])
        self.assertEqual("longitude latitude", solar_azimuth_angle.attrs["coordinates"])

    def assert_common_sensor_variables(self, ds, srf_size):
        scanline = ds.variables["scanline"]
        self.assertEqual((6,), scanline.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), scanline.data[3])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), scanline.attrs["_FillValue"])
        self.assertEqual("scanline_number", scanline.attrs["long_name"])
        self.assertEqual("count", scanline.attrs["units"])

        time = ds.variables["time"]
        self.assertEqual((6,), time.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.uint32), time.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint32), time.attrs["_FillValue"])
        self.assertEqual("time", time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", time.attrs["long_name"])
        self.assertEqual("s", time.attrs["units"])

        dq_bitmask = ds.variables["data_quality_bitmask"]
        self.assertEqual((6, 56), dq_bitmask.shape)
        self.assertEqual(0, dq_bitmask.data[0, 5])
        self.assertEqual("1, 2, 4, 8, 16", dq_bitmask.attrs["flag_masks"])
        self.assertEqual(
            "suspect_mirror suspect_geo suspect_time outlier_nos uncertainty_too_large",
            dq_bitmask.attrs["flag_meanings"])
        self.assertEqual("status_flag", dq_bitmask.attrs["standard_name"])
        self.assertEqual("longitude latitude", dq_bitmask.attrs["coordinates"])

        qual_scan_bitmask = ds.variables["quality_scanline_bitmask"]
        self.assertEqual((6,), qual_scan_bitmask.shape)
        self.assertEqual(0, qual_scan_bitmask.data[5])
        self.assertEqual("1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912 1073741824",
                         qual_scan_bitmask.attrs["flag_masks"])
        self.assertEqual(
            "do_not_use_scan time_sequence_error data_gap_preceding_scan no_calibration no_earth_location clock_update status_changed line_incomplete, time_field_bad time_field_bad_not_inf inconsistent_sequence scan_time_repeat uncalib_bad_time calib_few_scans uncalib_bad_prt calib_marginal_prt uncalib_channels uncalib_inst_mode quest_ant_black_body zero_loc bad_loc_time bad_loc_marginal bad_loc_reason bad_loc_ant reduced_context bad_temp_no_rself",
            qual_scan_bitmask.attrs["flag_meanings"])
        self.assertEqual("status_flag", qual_scan_bitmask.attrs["standard_name"])
        self.assertEqual("quality_indicator_bitfield", qual_scan_bitmask.attrs["long_name"])

        srf_weights = ds.variables["SRF_weights"]
        self.assertEqual((19, srf_size), srf_weights.shape)
        self.assertTrue(np.isnan(srf_weights.data[4, 6]))
        self.assertEqual("Spectral Response Function weights", srf_weights.attrs["long_name"])
        self.assertEqual("Per channel: weights for the relative spectral response function", srf_weights.attrs["description"])
        self.assertEqual(-32768, srf_weights.encoding['_FillValue'])
        self.assertEqual(0.000033, srf_weights.encoding['scale_factor'])

        srf_freqs = ds.variables["SRF_wavelengths"]
        self.assertEqual((19, srf_size), srf_freqs.shape)
        self.assertTrue(np.isnan(srf_freqs.data[5, 7]))
        self.assertEqual("Spectral Response Function wavelengths", srf_freqs.attrs["long_name"])
        self.assertEqual("Per channel: wavelengths for the relative spectral response function", srf_freqs.attrs["description"])
        self.assertEqual(-2147483648, srf_freqs.encoding['_FillValue'])
        self.assertEqual(0.0001, srf_freqs.encoding['scale_factor'])
        self.assertEqual("um", srf_freqs.attrs["units"])

        scnlin_map = ds.variables["scanline_map_to_origl1bfile"]
        self.assertEqual((6,), scnlin_map.shape)
        self.assertEqual(255, scnlin_map[1])
        self.assertEqual(255, scnlin_map.attrs['_FillValue'])
        self.assertEqual("Indicator of original file", scnlin_map.attrs['long_name'])
        self.assertEqual("Indicator for mapping each line to its corresponding original level 1b file. See global attribute 'source' for the filenames. 0 corresponds to 1st listed file, 1 to 2nd file.", scnlin_map.attrs['description'])

        scnlin_orig = ds.variables["scanline_origl1b"]
        self.assertEqual((6,), scnlin_orig.shape)
        self.assertEqual(-32767, scnlin_orig[2])
        self.assertEqual(-32767, scnlin_orig.attrs['_FillValue'])
        self.assertEqual("Original_Scan_line_number", scnlin_orig.attrs['long_name'])
        self.assertEqual("Original scan line numbers from corresponding l1b records", scnlin_orig.attrs['description'])

    def assert_extended_quality_flags(self, ds):
        chqualflags = ds.variables["quality_channel_bitmask"]
        self.assertEqual((6, 19), chqualflags.shape)
        self.assertEqual(0, chqualflags.data[1, 2])
        self.assertEqual("status_flag", chqualflags.attrs["standard_name"])
        self.assertEqual("channel_quality_flags_bitfield", chqualflags.attrs["long_name"])
        self.assertEqual("1, 2, 4, 8, 16", chqualflags.attrs["flag_masks"])
        self.assertEqual("do_not_use uncertainty_suspicious self_emission_fails calibration_impossible calibration_suspect", chqualflags.attrs["flag_meanings"])

    def assert_minor_frame_flags(self, ds):
        mnfrqualflags = ds.variables["mnfrqualflags"]
        self.assertEqual((6, 64), mnfrqualflags.shape)
        self.assertEqual(0, mnfrqualflags.data[2, 5])
        self.assertEqual("status_flag", mnfrqualflags.attrs["standard_name"])
        self.assertEqual("minor_frame_quality_flags_bitfield", mnfrqualflags.attrs["long_name"])

    def assert_easy_fcdr_uncertainties(self, ds, chunking=None):
        self._assert_3d_channel_variable(ds, "u_independent", "uncertainty from independent errors", chunking=chunking)
        self._assert_3d_channel_variable(ds, "u_structured", "uncertainty from structured errors", chunking=chunking)
        self._assert_3d_channel_variable(ds, "u_common", "uncertainty from common errors", chunking=chunking)

    def assert_coordinates(self, ds):
        x = ds.coords["x"]
        self.assertEqual((56,), x.shape)
        self.assertEqual(14, x[14])

        y = ds.coords["y"]
        self.assertEqual((6,), y.shape)
        self.assertEqual(5, y[5])

        channel = ds.coords["channel"]
        self.assertEqual((19,), channel.shape)
        self.assertEqual("Ch7", channel[6])

    def _assert_3d_channel_variable(self, ds, name, long_name, chunking=None):
        variable = ds.variables[name]
        self.assertEqual((19, 7, 56), variable.shape)
        self.assertTrue(np.isnan(variable.data[2, 5, 3]))
        self.assertEqual(65535, variable.encoding["_FillValue"])
        self.assertEqual(0.001, variable.encoding["scale_factor"])
        self.assertEqual(np.uint16, variable.encoding['dtype'])
        if chunking is not None:
            self.assertEqual(chunking, variable.encoding['chunksizes'])

        self.assertEqual(long_name, variable.attrs["long_name"])
        self.assertEqual("K", variable.attrs["units"])
        self.assertEqual(1, variable.attrs["valid_min"])
        self.assertEqual(65534, variable.attrs["valid_max"])
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])
