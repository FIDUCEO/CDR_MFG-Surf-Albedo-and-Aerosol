import unittest

import numpy as np
import xarray as xr

from fiduceo.common.test.assertions import Assertions
from fiduceo.common.writer.default_data import DefaultData
from fiduceo.fcdr.writer.templates.avhrr import AVHRR

CHUNKING = (1280, 409)
SRF_SIZE = 512
LUT_SIZE = 19
NUM_CHANNELS = 6


class AVHRRTest(unittest.TestCase):
    def test_add_original_variables(self):
        ds = xr.Dataset()
        AVHRR.add_original_variables(ds, 5, srf_size=SRF_SIZE)

        Assertions.assert_geolocation_variables(self, ds, 409, 5, chunking=CHUNKING)
        Assertions.assert_quality_flags(self, ds, 409, 5, chunking=CHUNKING)

        time = ds.variables["Time"]
        self.assertEqual((5,), time.shape)
        self.assertTrue(np.isnan(time.data[2]))
        self.assertTrue(np.isnan(time.attrs["_FillValue"]))
        self.assertEqual("time", time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", time.attrs["long_name"])
        self.assertEqual("s", time.attrs["units"])

        sat_zenith = ds.variables["satellite_zenith_angle"]
        self.assertEqual((5, 409), sat_zenith.shape)
        self.assertTrue(np.isnan(sat_zenith.data[0, 5]))
        self.assertEqual("sensor_zenith_angle", sat_zenith.attrs["standard_name"])
        self.assertEqual("degree", sat_zenith.attrs["units"])
        self.assertEqual(9000, sat_zenith.attrs["valid_max"])
        self.assertEqual(0, sat_zenith.attrs["valid_min"])
        self.assertEqual("longitude latitude", sat_zenith.attrs["coordinates"])
        self.assertEqual(np.int16, sat_zenith.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), sat_zenith.encoding['_FillValue'])
        self.assertEqual(0.01, sat_zenith.encoding['scale_factor'])
        self.assertEqual(0.0, sat_zenith.encoding['add_offset'])
        self.assertEqual(CHUNKING, sat_zenith.encoding["chunksizes"])

        sol_zenith = ds.variables["solar_zenith_angle"]
        self.assertEqual((5, 409), sol_zenith.shape)
        self.assertTrue(np.isnan(sol_zenith.data[0, 7]))
        self.assertEqual("solar_zenith_angle", sol_zenith.attrs["standard_name"])
        self.assertEqual("degree", sol_zenith.attrs["units"])
        self.assertEqual(18000, sol_zenith.attrs["valid_max"])
        self.assertEqual(0, sol_zenith.attrs["valid_min"])
        self.assertEqual("longitude latitude", sol_zenith.attrs["coordinates"])
        self.assertEqual(np.int16, sol_zenith.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), sol_zenith.encoding['_FillValue'])
        self.assertEqual(0.01, sol_zenith.encoding['scale_factor'])
        self.assertEqual(0.0, sol_zenith.encoding['add_offset'])
        self.assertEqual(CHUNKING, sol_zenith.encoding["chunksizes"])

        ch1 = ds.variables["Ch1"]
        self._assert_correct_refl_variable(ch1, "Channel 1 Reflectance")

        ch2 = ds.variables["Ch2"]
        self._assert_correct_refl_variable(ch2, "Channel 2 Reflectance")

        ch3a = ds.variables["Ch3a"]
        self._assert_correct_refl_variable(ch3a, "Channel 3a Reflectance")

        ch3b = ds.variables["Ch3b"]
        self._assert_correct_bt_variable(ch3b, "Channel 3b Brightness Temperature")

        ch4 = ds.variables["Ch4"]
        self._assert_correct_bt_variable(ch4, "Channel 4 Brightness Temperature")

        ch5 = ds.variables["Ch5"]
        self._assert_correct_bt_variable(ch5, "Channel 5 Brightness Temperature")

        dq_bitmask = ds.variables["data_quality_bitmask"]
        self.assertEqual((5, 409), dq_bitmask.shape)
        self.assertEqual(0, dq_bitmask.data[1, 23])
        self.assertEqual("status_flag", dq_bitmask.attrs["standard_name"])
        self.assertEqual("bitmask for quality per pixel", dq_bitmask.attrs["long_name"])
        self.assertEqual("1,2", dq_bitmask.attrs["flag_masks"])
        self.assertEqual("bad_geolocation_timing_err bad_calibration_radiometer_err", dq_bitmask.attrs["flag_meanings"])
        self.assertEqual("longitude latitude", dq_bitmask.attrs["coordinates"])

        qs_bitmask = ds.variables["quality_scanline_bitmask"]
        self.assertEqual((5,), qs_bitmask.shape)
        self.assertEqual(0, qs_bitmask.data[1])
        self.assertEqual("status_flag", qs_bitmask.attrs["standard_name"])
        self.assertEqual("bitmask for quality per scanline", qs_bitmask.attrs["long_name"])
        self.assertEqual("1,2,4,8,16,32,64", qs_bitmask.attrs["flag_masks"])
        self.assertEqual("do_not_use bad_time bad_navigation bad_calibration channel3a_present solar_contamination_failure solar_contamination", qs_bitmask.attrs["flag_meanings"])

        qc_bitmask = ds.variables["quality_channel_bitmask"]
        self.assertEqual((5, 6), qc_bitmask.shape)
        self.assertEqual(0, qc_bitmask.data[2, 1])
        self.assertEqual("status_flag", qc_bitmask.attrs["standard_name"])
        self.assertEqual("bitmask for quality per channel", qc_bitmask.attrs["long_name"])
        self.assertEqual("1,2", qc_bitmask.attrs["flag_masks"])
        self.assertEqual("bad_channel some_pixels_not_detected_2sigma", qc_bitmask.attrs["flag_meanings"])

        srf_weights = ds.variables["SRF_weights"]
        self.assertEqual((6, SRF_SIZE), srf_weights.shape)
        self.assertTrue(np.isnan(srf_weights.data[3, 5]))
        self.assertEqual("Spectral Response Function weights", srf_weights.attrs["long_name"])
        self.assertEqual("Per channel: weights for the relative spectral response function", srf_weights.attrs["description"])
        self.assertEqual(-32768, srf_weights.encoding['_FillValue'])
        self.assertEqual(0.000033, srf_weights.encoding['scale_factor'])

        srf_wls = ds.variables["SRF_wavelengths"]
        self.assertEqual((6, SRF_SIZE), srf_wls.shape)
        self.assertTrue(np.isnan(srf_wls.data[4, 6]))
        self.assertEqual("Spectral Response Function wavelengths", srf_wls.attrs["long_name"])
        self.assertEqual("Per channel: wavelengths for the relative spectral response function", srf_wls.attrs["description"])
        self.assertEqual(-2147483648, srf_wls.encoding['_FillValue'])
        self.assertEqual(0.0001, srf_wls.encoding['scale_factor'])
        self.assertEqual("um", srf_wls.attrs["units"])

        x = ds.coords["x"]
        self.assertEqual((409,), x.shape)
        self.assertEqual(13, x[13])

        y = ds.coords["y"]
        self.assertEqual((5,), y.shape)
        self.assertEqual(4, y[4])

        channel = ds.coords["channel"]
        self.assertEqual((6,), channel.shape)
        self.assertEqual("Ch5", channel[5])

        scnlin_map = ds.variables["scanline_map_to_origl1bfile"]
        self.assertEqual((5,), scnlin_map.shape)
        self.assertEqual(255, scnlin_map[0])
        self.assertEqual(255, scnlin_map.attrs['_FillValue'])
        self.assertEqual("Indicator of original file", scnlin_map.attrs['long_name'])
        self.assertEqual("Indicator for mapping each line to its corresponding original level 1b file. See global attribute 'source' for the filenames. 0 corresponds to 1st listed file, 1 to 2nd file.", scnlin_map.attrs['description'])

        scnlin_orig = ds.variables["scanline_origl1b"]
        self.assertEqual((5,), scnlin_orig.shape)
        self.assertEqual(-32767, scnlin_orig[1])
        self.assertEqual(-32767, scnlin_orig.attrs['_FillValue'])
        self.assertEqual("Original_Scan_line_number", scnlin_orig.attrs['long_name'])
        self.assertEqual("Original scan line numbers from corresponding l1b records", scnlin_orig.attrs['description'])

    def test_get_swath_width(self):
        self.assertEqual(409, AVHRR.get_swath_width())

    def test_add_easy_fcdr_variables(self):
        ds = xr.Dataset()
        AVHRR.add_easy_fcdr_variables(ds, 5, lut_size=LUT_SIZE, corr_dx=12, corr_dy=13)

        self._assert_correct_refl_uncertainty_variable(ds, "u_independent_Ch1", long_name="independent uncertainty per pixel for channel 1", units="percent", valid_min=10, valid_max=1000)
        self._assert_correct_refl_uncertainty_variable(ds, "u_structured_Ch1", long_name="structured uncertainty per pixel for channel 1", units="percent", valid_min=3, valid_max=5)
        self._assert_correct_refl_uncertainty_variable(ds, "u_common_Ch1", long_name="common uncertainty per pixel for channel 1", units="percent", valid_min=10, valid_max=1000)
        self._assert_correct_refl_uncertainty_variable(ds, "u_independent_Ch2", long_name="independent uncertainty per pixel for channel 2", units="percent", valid_min=10, valid_max=1000)
        self._assert_correct_refl_uncertainty_variable(ds, "u_structured_Ch2", long_name="structured uncertainty per pixel for channel 2", units="percent", valid_min=3, valid_max=5)
        self._assert_correct_refl_uncertainty_variable(ds, "u_common_Ch2", long_name="common uncertainty per pixel for channel 2", units="percent", valid_min=10, valid_max=1000)
        self._assert_correct_refl_uncertainty_variable(ds, "u_independent_Ch3a", long_name="independent uncertainty per pixel for channel 3a", units="percent", valid_min=10, valid_max=1000)
        self._assert_correct_refl_uncertainty_variable(ds, "u_structured_Ch3a", long_name="structured uncertainty per pixel for channel 3a", units="percent", valid_min=3, valid_max=5)
        self._assert_correct_refl_uncertainty_variable(ds, "u_common_Ch3a", long_name="common uncertainty per pixel for channel 3a", units="percent", valid_min=10, valid_max=1000)

        self._assert_correct_bt_uncertainty_variable(ds, "u_independent_Ch3b", long_name="independent uncertainty per pixel for channel 3b")
        self._assert_correct_bt_uncertainty_variable(ds, "u_structured_Ch3b", long_name="structured uncertainty per pixel for channel 3b")
        self._assert_correct_bt_uncertainty_variable(ds, "u_common_Ch3b", long_name="common uncertainty per pixel for channel 3b")
        self._assert_correct_bt_uncertainty_variable(ds, "u_independent_Ch4", long_name="independent uncertainty per pixel for channel 4")
        self._assert_correct_bt_uncertainty_variable(ds, "u_structured_Ch4", long_name="structured uncertainty per pixel for channel 4")
        self._assert_correct_bt_uncertainty_variable(ds, "u_common_Ch4", long_name="common uncertainty per pixel for channel 4")
        self._assert_correct_bt_uncertainty_variable(ds, "u_independent_Ch5", long_name="independent uncertainty per pixel for channel 5")
        self._assert_correct_bt_uncertainty_variable(ds, "u_structured_Ch5", long_name="structured uncertainty per pixel for channel 5")
        self._assert_correct_bt_uncertainty_variable(ds, "u_common_Ch5", long_name="common uncertainty per pixel for channel 5")

        Assertions.assert_correlation_matrices(self, ds, NUM_CHANNELS)
        Assertions.assert_lookup_tables(self, ds, NUM_CHANNELS, LUT_SIZE)
        Assertions.assert_correlation_coefficients(self, ds, NUM_CHANNELS, 12, 13)

    def test_add_full_fcdr_variables(self):
        ds = xr.Dataset()
        AVHRR.add_full_fcdr_variables(ds, 5)

        u_latitude = ds.variables["u_latitude"]
        self.assertEqual((5, 409), u_latitude.shape)
        self.assertTrue(np.isnan(u_latitude.data[0, 34]))
        self.assertTrue(np.isnan(u_latitude.attrs["_FillValue"]))
        self.assertEqual("uncertainty of latitude", u_latitude.attrs["long_name"])
        self.assertEqual("degree", u_latitude.attrs["units"])
        self.assertEqual("longitude latitude", u_latitude.attrs["coordinates"])
        self.assertEqual(CHUNKING, u_latitude.encoding["chunksizes"])

        u_longitude = ds.variables["u_longitude"]
        self.assertEqual((5, 409), u_longitude.shape)
        self.assertTrue(np.isnan(u_longitude.data[1, 35]))
        self.assertTrue(np.isnan(u_longitude.attrs["_FillValue"]))
        self.assertEqual("uncertainty of longitude", u_longitude.attrs["long_name"])
        self.assertEqual("degree", u_longitude.attrs["units"])
        self.assertEqual("longitude latitude", u_longitude.attrs["coordinates"])
        self.assertEqual(CHUNKING, u_longitude.encoding["chunksizes"])

        u_time = ds.variables["u_time"]
        self.assertEqual((5,), u_time.shape)
        self.assertTrue(np.isnan(u_time.data[1]))
        self.assertTrue(np.isnan(u_time.attrs["_FillValue"]))
        self.assertEqual("uncertainty of acquisition time", u_time.attrs["long_name"])
        self.assertEqual("s", u_time.attrs["units"])

        u_sat_azimuth = ds.variables["u_satellite_azimuth_angle"]
        self.assertEqual((5, 409), u_sat_azimuth.shape)
        self.assertTrue(np.isnan(u_sat_azimuth.data[2, 36]))
        self.assertTrue(np.isnan(u_sat_azimuth.attrs["_FillValue"]))
        self.assertEqual("uncertainty of satellite azimuth angle", u_sat_azimuth.attrs["long_name"])
        self.assertEqual("degree", u_sat_azimuth.attrs["units"])
        self.assertEqual("longitude latitude", u_sat_azimuth.attrs["coordinates"])
        self.assertEqual(CHUNKING, u_sat_azimuth.encoding["chunksizes"])

        u_sat_zenith = ds.variables["u_satellite_zenith_angle"]
        self.assertEqual((5, 409), u_sat_zenith.shape)
        self.assertTrue(np.isnan(u_sat_zenith.data[2, 36]))
        self.assertTrue(np.isnan(u_sat_zenith.attrs["_FillValue"]))
        self.assertEqual("uncertainty of satellite zenith angle", u_sat_zenith.attrs["long_name"])
        self.assertEqual("degree", u_sat_zenith.attrs["units"])
        self.assertEqual("longitude latitude", u_sat_zenith.attrs["coordinates"])
        self.assertEqual(CHUNKING, u_sat_zenith.encoding["chunksizes"])

        u_sol_azimuth = ds.variables["u_solar_azimuth_angle"]
        self.assertEqual((5, 409), u_sol_azimuth.shape)
        self.assertTrue(np.isnan(u_sol_azimuth.data[2, 36]))
        self.assertTrue(np.isnan(u_sol_azimuth.attrs["_FillValue"]))
        self.assertEqual("uncertainty of solar azimuth angle", u_sol_azimuth.attrs["long_name"])
        self.assertEqual("degree", u_sol_azimuth.attrs["units"])
        self.assertEqual("longitude latitude", u_sol_azimuth.attrs["coordinates"])
        self.assertEqual(CHUNKING, u_sol_azimuth.encoding["chunksizes"])

        u_sol_zenith = ds.variables["u_solar_zenith_angle"]
        self.assertEqual((5, 409), u_sol_zenith.shape)
        self.assertTrue(np.isnan(u_sol_zenith.data[2, 36]))
        self.assertTrue(np.isnan(u_sol_zenith.attrs["_FillValue"]))
        self.assertEqual("uncertainty of solar zenith angle", u_sol_zenith.attrs["long_name"])
        self.assertEqual("degree", u_sol_zenith.attrs["units"])
        self.assertEqual("longitude latitude", u_sol_zenith.attrs["coordinates"])
        self.assertEqual(CHUNKING, u_sol_zenith.encoding["chunksizes"])

        prt_c = ds.variables["PRT_C"]
        self.assertEqual((5, 3), prt_c.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), prt_c.data[3, 2])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), prt_c.attrs["_FillValue"])
        self.assertEqual("Prt counts", prt_c.attrs["long_name"])
        self.assertEqual("count", prt_c.attrs["units"])

        u_prt = ds.variables["u_prt"]
        self.assertEqual((5, 3), u_prt.shape)
        self.assertTrue(np.isnan(u_prt.data[4, 0]))
        self.assertTrue(np.isnan(u_prt.attrs["_FillValue"]))
        self.assertEqual("Uncertainty on the PRT counts", u_prt.attrs["long_name"])
        self.assertEqual("count", u_prt.attrs["units"])
        self.assertEqual("rectangle_absolute", u_prt.attrs["pixel_correlation_form"])
        self.assertEqual("pixel", u_prt.attrs["pixel_correlation_units"])
        self.assertEqual([-np.inf, np.inf], u_prt.attrs["pixel_correlation_scales"])
        self.assertEqual("rectangle_absolute", u_prt.attrs["scan_correlation_form"])
        self.assertEqual("line", u_prt.attrs["scan_correlation_units"])
        self.assertEqual([-np.inf, np.inf], u_prt.attrs["scan_correlation_scales"])
        self.assertEqual("rectangle", u_prt.attrs["pdf_shape"])
        self.assertEqual(0.1, u_prt.attrs["pdf_parameter"])

        r_ict = ds.variables["R_ICT"]
        self.assertEqual((5, 3), r_ict.shape)
        self.assertTrue(np.isnan(r_ict.data[0, 1]))
        self.assertTrue(np.isnan(r_ict.attrs["_FillValue"]))
        self.assertEqual("Radiance of the PRT", r_ict.attrs["long_name"])
        self.assertEqual("mW m^-2 sr^-1 cm", r_ict.attrs["units"])

        t_instr = ds.variables["T_instr"]
        self.assertEqual((5,), t_instr.shape)
        self.assertTrue(np.isnan(t_instr.data[4]))
        self.assertTrue(np.isnan(t_instr.attrs["_FillValue"]))
        self.assertEqual("Instrument temperature", t_instr.attrs["long_name"])
        self.assertEqual("K", t_instr.attrs["units"])

        self._assert_correct_counts_variable(ds, "Ch1_Csp", "Ch1 Space counts")
        self._assert_correct_counts_variable(ds, "Ch2_Csp", "Ch2 Space counts")
        self._assert_correct_counts_variable(ds, "Ch3a_Csp", "Ch3a Space counts")
        self._assert_correct_counts_variable(ds, "Ch3b_Csp", "Ch3b Space counts")
        self._assert_correct_counts_variable(ds, "Ch4_Csp", "Ch4 Space counts")
        self._assert_correct_counts_variable(ds, "Ch5_Csp", "Ch5 Space counts")

        self._assert_correct_counts_variable(ds, "Ch3b_Cict", "Ch3b ICT counts")
        self._assert_correct_counts_variable(ds, "Ch4_Cict", "Ch4 ICT counts")
        self._assert_correct_counts_variable(ds, "Ch5_Cict", "Ch5 ICT counts")

        self._assert_correct_counts_variable(ds, "Ch1_Ce", "Ch1 Earth counts")
        self._assert_correct_counts_variable(ds, "Ch2_Ce", "Ch2 Earth counts")
        self._assert_correct_counts_variable(ds, "Ch3a_Ce", "Ch3a Earth counts")
        self._assert_correct_counts_variable(ds, "Ch3b_Ce", "Ch3b Earth counts")
        self._assert_correct_counts_variable(ds, "Ch4_Ce", "Ch4 Earth counts")
        self._assert_correct_counts_variable(ds, "Ch5_Ce", "Ch5 Earth counts")

        self._assert_correct_counts_uncertainty_variable(ds, "Ch1_u_Csp", "Ch1 Uncertainty on space counts")
        self._assert_counts_correlation(ds, "Ch1_u_Csp")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch2_u_Csp", "Ch2 Uncertainty on space counts")
        self._assert_counts_correlation(ds, "Ch2_u_Csp")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch3a_u_Csp", "Ch3a Uncertainty on space counts")
        self._assert_counts_correlation(ds, "Ch3a_u_Csp")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch3b_u_Csp", "Ch3b Uncertainty on space counts")
        self._assert_counts_correlation(ds, "Ch3b_u_Csp")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch4_u_Csp", "Ch4 Uncertainty on space counts")
        self._assert_counts_correlation(ds, "Ch4_u_Csp")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch5_u_Csp", "Ch5 Uncertainty on space counts")
        self._assert_counts_correlation(ds, "Ch4_u_Csp")

        self._assert_correct_counts_uncertainty_variable(ds, "Ch3b_u_Cict", "Ch3b Uncertainty on ICT counts")
        self._assert_counts_correlation(ds, "Ch3b_u_Cict")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch4_u_Cict", "Ch4 Uncertainty on ICT counts")
        self._assert_counts_correlation(ds, "Ch4_u_Cict")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch5_u_Cict", "Ch5 Uncertainty on ICT counts")
        self._assert_counts_correlation(ds, "Ch5_u_Cict")

        self._assert_correct_counts_uncertainty_variable(ds, "Ch1_u_Ce", "Ch1 Uncertainty on earth counts")
        self._assert_earth_counts_pdf(ds, "Ch1_u_Ce")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch2_u_Ce", "Ch2 Uncertainty on earth counts")
        self._assert_earth_counts_pdf(ds, "Ch2_u_Ce")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch3a_u_Ce", "Ch3a Uncertainty on earth counts")
        self._assert_earth_counts_pdf(ds, "Ch3a_u_Ce")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch3b_u_Ce", "Ch3b Uncertainty on earth counts")
        self._assert_earth_counts_pdf(ds, "Ch3b_u_Ce")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch4_u_Ce", "Ch4 Uncertainty on earth counts")
        self._assert_earth_counts_pdf(ds, "Ch4_u_Ce")
        self._assert_correct_counts_uncertainty_variable(ds, "Ch5_u_Ce", "Ch5 Uncertainty on earth counts")
        self._assert_earth_counts_pdf(ds, "Ch5_u_Ce")

        self._assert_correct_refl_uncertainty_variable(ds, "Ch1_u_Refl", long_name="Ch1 Total uncertainty on toa reflectance", units="percent")
        self._assert_correct_refl_uncertainty_variable(ds, "Ch2_u_Refl", long_name="Ch2 Total uncertainty on toa reflectance", units="percent")
        self._assert_correct_refl_uncertainty_variable(ds, "Ch3a_u_Refl", long_name="Ch3a Total uncertainty on toa reflectance", units="percent")

        self._assert_correct_bt_uncertainty_variable(ds, "Ch3b_u_Bt", long_name="Ch3b Total uncertainty on brightness temperature")
        self._assert_correct_bt_uncertainty_variable(ds, "Ch4_u_Bt", long_name="Ch4 Total uncertainty on brightness temperature")
        self._assert_correct_bt_uncertainty_variable(ds, "Ch5_u_Bt", long_name="Ch5 Total uncertainty on brightness temperature")

        self._assert_correct_bt_uncertainty_variable(ds, "Ch3b_ur_Bt", long_name="Ch3b Random uncertainty on brightness temperature")
        self._assert_correct_bt_uncertainty_variable(ds, "Ch4_ur_Bt", long_name="Ch4 Random uncertainty on brightness temperature")
        self._assert_correct_bt_uncertainty_variable(ds, "Ch5_ur_Bt", long_name="Ch5 Random uncertainty on brightness temperature")

        self._assert_correct_bt_uncertainty_variable(ds, "Ch3b_us_Bt", long_name="Ch3b Systematic uncertainty on brightness temperature")
        self._assert_correct_bt_uncertainty_variable(ds, "Ch4_us_Bt", long_name="Ch4 Systematic uncertainty on brightness temperature")
        self._assert_correct_bt_uncertainty_variable(ds, "Ch5_us_Bt", long_name="Ch5 Systematic uncertainty on brightness temperature")

    def test_add_template_key(self):
        ds = xr.Dataset()

        AVHRR.add_template_key(ds)

        self.assertEqual("AVHRR", ds.attrs["template_key"])

    def _assert_earth_counts_pdf(self, ds, name):
        variable = ds.variables[name]
        self.assertEqual("digitised_gaussian", variable.attrs["pdf_shape"])

    def _assert_counts_correlation(self, ds, name):
        variable = ds.variables[name]
        self.assertEqual("rectangle_absolute", variable.attrs["pixel_correlation_form"])
        self.assertEqual("pixel", variable.attrs["pixel_correlation_units"])
        self.assertEqual([-np.inf, np.inf], variable.attrs["pixel_correlation_scales"])
        self.assertEqual("triangle_relative", variable.attrs["scan_correlation_form"])
        self.assertEqual("line", variable.attrs["scan_correlation_units"])
        self.assertEqual([-25, 25], variable.attrs["scan_correlation_scales"])
        self.assertEqual("digitised_gaussian", variable.attrs["pdf_shape"])

    def _assert_correct_counts_variable(self, ds, name, long_name):
        variable = ds.variables[name]
        self.assertEqual((5, 409), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), variable.data[3, 306])
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), variable.attrs["_FillValue"])
        self.assertEqual(long_name, variable.attrs["long_name"])
        self.assertEqual("count", variable.attrs["units"])
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])
        self.assertEqual(CHUNKING, variable.encoding["chunksizes"])

    def _assert_correct_counts_uncertainty_variable(self, ds, name, standard_name):
        variable = ds.variables[name]
        self.assertEqual((5, 409), variable.shape)
        self.assertTrue(np.isnan(variable.data[4, 307]))
        self.assertTrue(np.isnan(variable.attrs["_FillValue"]))
        self.assertEqual(standard_name, variable.attrs["long_name"])
        self.assertEqual("count", variable.attrs["units"])
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])
        self.assertEqual(CHUNKING, variable.encoding["chunksizes"])

    def _assert_correct_refl_uncertainty_variable(self, ds, name, standard_name=None, units=None, long_name=None, valid_min=None, valid_max=None):
        variable = ds.variables[name]
        self.assertEqual((5, 409), variable.shape)
        self.assertTrue(np.isnan(variable.data[4, 307]))
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])
        self.assertEqual(-32767, variable.encoding["_FillValue"])
        self.assertEqual(np.int16, variable.encoding["dtype"])
        self.assertEqual(CHUNKING, variable.encoding["chunksizes"])

        if valid_min is not None:
            self.assertEqual(valid_min, variable.attrs["valid_min"])

        if valid_max is not None:
            self.assertEqual(valid_max, variable.attrs["valid_max"])

        if standard_name is not None:
            self.assertEqual(standard_name, variable.attrs["standard_name"])

        if long_name is not None:
            self.assertEqual(long_name, variable.attrs["long_name"])

        if units is not None:
            self.assertEqual(units, variable.attrs["units"])

    def _assert_correct_bt_uncertainty_variable(self, ds, name, standard_name=None, long_name=None):
        variable = ds.variables[name]
        self.assertEqual((5, 409), variable.shape)
        self.assertTrue(np.isnan(variable.data[4, 307]))
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])
        self.assertEqual(-32767, variable.encoding["_FillValue"])
        self.assertEqual(np.int16, variable.encoding["dtype"])
        self.assertEqual(CHUNKING, variable.encoding["chunksizes"])

        if standard_name is not None:
            self.assertEqual(standard_name, variable.attrs["standard_name"])

        if long_name is not None:
            self.assertEqual(long_name, variable.attrs["long_name"])

        self.assertEqual("K", variable.attrs["units"])

    def _assert_correct_refl_variable(self, variable, long_name):
        self.assertEqual((5, 409), variable.shape)
        self.assertTrue(np.isnan(variable.data[0, 8]))
        self.assertEqual("toa_reflectance", variable.attrs["standard_name"])
        self.assertEqual(long_name, variable.attrs["long_name"])
        self.assertEqual("1", variable.attrs["units"])
        self.assertEqual(np.int16, variable.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), variable.encoding['_FillValue'])
        self.assertEqual(0.0001, variable.encoding['scale_factor'])
        self.assertEqual(0.0, variable.encoding['add_offset'])
        self.assertEqual(CHUNKING, variable.encoding["chunksizes"])
        self.assertEqual(15000, variable.attrs["valid_max"])
        self.assertEqual(0, variable.attrs["valid_min"])
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])

    def _assert_correct_bt_variable(self, variable, long_name):
        self.assertEqual((5, 409), variable.shape)
        self.assertTrue(np.isnan(variable.data[0, 8]))
        self.assertEqual("toa_brightness_temperature", variable.attrs["standard_name"])
        self.assertEqual(long_name, variable.attrs["long_name"])
        self.assertEqual(np.int16, variable.encoding['dtype'])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), variable.encoding['_FillValue'])
        self.assertEqual(0.01, variable.encoding['scale_factor'])
        self.assertEqual(273.15, variable.encoding['add_offset'])
        self.assertEqual(CHUNKING, variable.encoding["chunksizes"])
        self.assertEqual("K", variable.attrs["units"])
        self.assertEqual(10000, variable.attrs["valid_max"])
        self.assertEqual(-20000, variable.attrs["valid_min"])
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])

