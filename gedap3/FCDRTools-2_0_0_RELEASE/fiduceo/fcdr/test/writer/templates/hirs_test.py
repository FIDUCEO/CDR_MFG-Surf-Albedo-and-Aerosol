import unittest

import numpy as np
import xarray as xr

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.fcdr.writer.templates.hirs import HIRS


class HIRSTest(unittest.TestCase):
    def test_add_full_fcdr_variables(self):
        ds = xr.Dataset()
        HIRS.add_full_fcdr_variables(ds, 7)

        c_earth = ds.variables["c_earth"]
        self.assertEqual((20, 7, 56), c_earth.shape)
        self.assertEqual(65535, c_earth.data[0, 2, 3])
        self.assertEqual(65535, c_earth.attrs["_FillValue"])
        self.assertEqual("counts_earth", c_earth.attrs["long_name"])
        self.assertEqual("count", c_earth.attrs["units"])
        self.assertEqual("scnlinf quality_scanline_bitmask quality_channel_bitmask mnfrqualflags", c_earth.attrs["ancilliary_variables"])

        l_earth = ds.variables["L_earth"]
        self.assertEqual((20, 7, 56), l_earth.shape)
        self.assertTrue(np.isnan(l_earth.data[0, 2, 4]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint32), l_earth.encoding["_FillValue"])
        self.assertEqual(0.0001, l_earth.encoding["scale_factor"])
        self.assertEqual(np.uint32, l_earth.encoding["dtype"])
        self.assertEqual("toa_outgoing_inband_radiance", l_earth.attrs["standard_name"])
        self.assertEqual("W/Hz/m ** 2/sr", l_earth.attrs["units"])
        self.assertEqual("Channel radiance, NOAA/EUMETSAT calibrated", l_earth.attrs["long_name"])
        self.assertEqual("scnlinf quality_scanline_bitmask quality_channel_bitmask mnfrqualflags", l_earth.attrs["ancilliary_variables"])

        u_lat = ds.variables["u_lat"]
        self.assertEqual((7, 56), u_lat.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_lat.data[3, 3])
        self.assertEqual(65535, u_lat.encoding["_FillValue"])
        self.assertEqual(0.01, u_lat.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_lat.encoding["dtype"])
        self.assertEqual("uncertainty_latitude", u_lat.attrs["standard_name"])
        self.assertEqual("degree", u_lat.attrs["units"])

        u_lon = ds.variables["u_lon"]
        self.assertEqual((7, 56), u_lon.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_lon.data[4, 4])
        self.assertEqual(65535, u_lon.encoding["_FillValue"])
        self.assertEqual(0.01, u_lon.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_lon.encoding["dtype"])
        self.assertEqual("uncertainty_longitude", u_lon.attrs["standard_name"])
        self.assertEqual("degree", u_lon.attrs["units"])

        u_time = ds.variables["u_time"]
        self.assertEqual((7, 56), u_time.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_time.data[5, 5])
        self.assertEqual(65535, u_time.encoding["_FillValue"])
        self.assertEqual(0.01, u_time.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_time.encoding["dtype"])
        self.assertEqual("uncertainty_time", u_time.attrs["standard_name"])
        self.assertEqual("s", u_time.attrs["units"])

        u_c_earth = ds.variables["u_c_earth"]
        self.assertEqual((19, 337), u_c_earth.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_c_earth.data[6, 6])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_c_earth.attrs["_FillValue"])
        self.assertEqual("uncertainty counts for Earth views", u_c_earth.attrs["long_name"])
        self.assertEqual("count", u_c_earth.attrs["units"])
        self.assertEqual("u_c_earth_chan_corr", u_c_earth.attrs["ancilliary_variables"])
        self.assertEqual("all", u_c_earth.attrs["channels_affected"])
        self.assertEqual("C_E", u_c_earth.attrs["parameter"])
        self.assertEqual("gaussian", u_c_earth.attrs["pdf_shape"])

        u_L_earth_idependent = ds.variables["u_L_earth_independent"]
        self.assertEqual((20, 7, 56), u_L_earth_idependent.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_L_earth_idependent.data[7, 0, 7])
        self.assertEqual(4294967295, u_L_earth_idependent.encoding["_FillValue"])
        self.assertEqual(0.01, u_L_earth_idependent.encoding["scale_factor"])
        self.assertEqual(np.uint32, u_L_earth_idependent.encoding["dtype"])
        self.assertEqual("uncertainty_radiance_Earth_random", u_L_earth_idependent.attrs["standard_name"])
        self.assertEqual("mW m^-2 sr^-1 cm", u_L_earth_idependent.attrs["units"])

        u_L_earth_structured = ds.variables["u_L_earth_structured"]
        self.assertEqual((20, 7, 56), u_L_earth_structured.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_L_earth_structured.data[8, 1, 8])
        self.assertEqual(4294967295, u_L_earth_structured.encoding["_FillValue"])
        self.assertEqual(0.01, u_L_earth_structured.encoding["scale_factor"])
        self.assertEqual(np.uint32, u_L_earth_structured.encoding["dtype"])
        self.assertEqual("uncertainty_radiance_Earth_structured", u_L_earth_structured.attrs["standard_name"])
        self.assertEqual("mW m^-2 sr^-1 cm", u_L_earth_structured.attrs["units"])

        u_L_earth_sys = ds.variables["u_L_earth_systematic"]
        self.assertEqual((20, 7, 56), u_L_earth_sys.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_L_earth_sys.data[9, 2, 9])
        self.assertEqual(4294967295, u_L_earth_sys.encoding["_FillValue"])
        self.assertEqual(0.01, u_L_earth_sys.encoding["scale_factor"])
        self.assertEqual(np.uint32, u_L_earth_sys.encoding["dtype"])
        self.assertEqual("uncertainty_radiance_Earth_systematic", u_L_earth_sys.attrs["standard_name"])
        self.assertEqual("mW m^-2 sr^-1 cm", u_L_earth_sys.attrs["units"])

        u_L_earth_total = ds.variables["u_L_earth_total"]
        self.assertEqual((20, 7, 56), u_L_earth_total.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_L_earth_total.data[10, 3, 10])
        self.assertEqual(4294967295, u_L_earth_total.encoding["_FillValue"])
        self.assertEqual(0.01, u_L_earth_total.encoding["scale_factor"])
        self.assertEqual(np.uint32, u_L_earth_total.encoding["dtype"])
        self.assertEqual("uncertainty_radiance_Earth_total", u_L_earth_total.attrs["standard_name"])
        self.assertEqual("mW m^-2 sr^-1 cm", u_L_earth_total.attrs["units"])

        S_u_L_earth = ds.variables["S_u_L_earth"]
        self.assertEqual((20, 20), S_u_L_earth.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), S_u_L_earth.data[11, 4])
        self.assertEqual(4294967295, S_u_L_earth.encoding["_FillValue"])
        self.assertEqual(0.01, S_u_L_earth.encoding["scale_factor"])
        self.assertEqual(np.uint32, S_u_L_earth.encoding["dtype"])
        self.assertEqual("covariance_radiance_Earth", S_u_L_earth.attrs["standard_name"])

        u_bt_random = ds.variables["u_bt_random"]
        self.assertEqual((19, 7, 56), u_bt_random.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_bt_random.data[13, 6, 13])
        self.assertEqual(65535, u_bt_random.encoding["_FillValue"])
        self.assertEqual(0.01, u_bt_random.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_bt_random.encoding["dtype"])
        self.assertEqual("uncertainty_bt_random", u_bt_random.attrs["standard_name"])
        self.assertEqual("K", u_bt_random.attrs["units"])

        u_bt_structured = ds.variables["u_bt_structured"]
        self.assertEqual((19, 7, 56), u_bt_structured.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_bt_structured.data[14, 0, 14])
        self.assertEqual(65535, u_bt_structured.encoding["_FillValue"])
        self.assertEqual(0.01, u_bt_structured.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_bt_structured.encoding["dtype"])
        self.assertEqual("uncertainty_bt_structured", u_bt_structured.attrs["standard_name"])
        self.assertEqual("K", u_bt_structured.attrs["units"])

        u_bt_sys = ds.variables["u_bt_systematic"]
        self.assertEqual((19, 7, 56), u_bt_sys.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_bt_sys.data[15, 1, 15])
        self.assertEqual(65535, u_bt_sys.encoding["_FillValue"])
        self.assertEqual(0.01, u_bt_sys.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_bt_sys.encoding["dtype"])
        self.assertEqual("uncertainty_bt_systematic", u_bt_sys.attrs["standard_name"])
        self.assertEqual("K", u_bt_sys.attrs["units"])

        u_bt_total = ds.variables["u_bt_total"]
        self.assertEqual((19, 7, 56), u_bt_total.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_bt_total.data[15, 1, 15])
        self.assertEqual(65535, u_bt_total.encoding["_FillValue"])
        self.assertEqual(0.01, u_bt_total.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_bt_total.encoding["dtype"])
        self.assertEqual("uncertainty_bt_total", u_bt_total.attrs["standard_name"])
        self.assertEqual("K", u_bt_total.attrs["units"])

        S_bt = ds.variables["S_bt"]
        self.assertEqual((20, 20), S_bt.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), S_bt.data[12, 3])
        self.assertEqual(65535, S_bt.encoding["_FillValue"])
        self.assertEqual(0.01, S_bt.encoding["scale_factor"])
        self.assertEqual(np.uint16, S_bt.encoding["dtype"])
        self.assertEqual("covariance_brightness_temperature", S_bt.attrs["standard_name"])

        l1b_calcof = ds.variables["l1b_calcof"]
        self.assertEqual((3, 7), l1b_calcof.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), l1b_calcof.data[0, 2])
        self.assertEqual(-2147483647, l1b_calcof.encoding["_FillValue"])
        self.assertEqual(0.01, l1b_calcof.encoding["scale_factor"])
        self.assertEqual(np.int32, l1b_calcof.encoding["dtype"])
        self.assertEqual("calibration_coefficients", l1b_calcof.attrs["standard_name"])

        self._assert_line_int32_variable(ds, "navigation_status", standard_name="status_flag", long_name="Navigation status bit field", orig_name="hrs_navstat")

        variable = self._assert_line_uint16_variable(ds, "platform_altitude", long_name="Platform altitude", orig_name="hrs_scalti")
        self.assertEqual("km", variable.attrs["units"])

        variable = self._assert_line_scaled_int16_variable(ds, "platform_pitch_angle", long_name="Platform pitch angle", orig_name="hrs_pitchang")
        self.assertEqual("degree", variable.attrs["units"])

        variable = self._assert_line_scaled_int16_variable(ds, "platform_roll_angle", long_name="Platform roll angle", orig_name="hrs_rollang")
        self.assertEqual("degree", variable.attrs["units"])

        variable = self._assert_line_scaled_int16_variable(ds, "platform_yaw_angle", long_name="Platform yaw angle", orig_name="hrs_yawang")
        self.assertEqual("degree", variable.attrs["units"])

        self._assert_line_int32_variable(ds, "quality_flags", standard_name="status_flag", long_name="Quality indicator bit field", orig_name="hrs_qualind")

        scan_angles = ds.variables["scan_angles"]
        self.assertEqual((7, 168), scan_angles.shape)
        self.assertTrue(np.isnan(scan_angles.data[4, 18]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), scan_angles.encoding["_FillValue"])
        self.assertEqual(0.01, scan_angles.encoding["scale_factor"])
        self.assertEqual(np.uint16, scan_angles.encoding["dtype"])
        self.assertEqual("Scan angles", scan_angles.attrs["long_name"])
        self.assertEqual("hrs_ang", scan_angles.attrs["orig_name"])
        self.assertEqual("degree", scan_angles.attrs["units"])

        self._assert_line_int16_variable(ds, "l1b_scanline_number", long_name="scanline number", orig_name="hrs_scnlin")
        self._assert_line_int8_variable(ds, "scanline_position", long_name="Scanline position number in 32 second cycle", orig_name="hrs_scnpos")

        sec_o_cal_coeff = ds.variables["l1b_second_original_calibration_coefficients"]
        self.assertEqual((7, 60), sec_o_cal_coeff.shape)
        self.assertTrue(np.isnan(sec_o_cal_coeff.data[4, 18]))
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), sec_o_cal_coeff.encoding["_FillValue"])
        self.assertEqual(0.01, sec_o_cal_coeff.encoding["scale_factor"])
        self.assertEqual(np.int32, sec_o_cal_coeff.encoding["dtype"])
        self.assertEqual("Second original calibration coefficients (unsorted)", sec_o_cal_coeff.attrs["long_name"])
        self.assertEqual("hrs_scalcof", sec_o_cal_coeff.attrs["orig_name"])

        self._assert_line_counts_variable(ds, "Tc_baseplate", "temperature_baseplate_counts")
        self._assert_line_counts_variable(ds, "Tc_ch", "temperature_coolerhousing_counts")
        self._assert_line_counts_variable(ds, "Tc_elec", "temperature_electronics_counts")
        self._assert_line_counts_variable(ds, "Tc_fsr", "temperature_first_stage_radiator_counts")
        self._assert_line_counts_variable(ds, "Tc_fwh", "temperature_filter_wheel_housing_counts")
        self._assert_line_counts_variable(ds, "Tc_fwm", "temperature_filter_wheel_monitor_counts")
        self._assert_line_counts_variable(ds, "Tc_icct", "temperature_internal_cold_calibration_target_counts")
        self._assert_line_counts_variable(ds, "Tc_iwct", "temperature_internal_warm_calibration_target_counts")
        self._assert_line_counts_variable(ds, "Tc_patch_exp", "temperature_patch_expanded_scale_counts")
        self._assert_line_counts_variable(ds, "Tc_patch_full", "temperature_patch_full_range_counts")
        self._assert_line_counts_variable(ds, "Tc_tlscp_prim", "temperature_telescope_primary_counts")
        self._assert_line_counts_variable(ds, "Tc_tlscp_sec", "temperature_telescope_secondary_counts")
        self._assert_line_counts_variable(ds, "Tc_tlscp_tert", "temperature_telescope_tertiary_counts")
        self._assert_line_counts_variable(ds, "Tc_scanmirror", "temperature_scanmirror_counts")
        self._assert_line_counts_variable(ds, "Tc_scanmotor", "temperature_scanmotor_counts")

        self._assert_line_counts_uncertainty_variable_uint16(ds, "u_Tc_baseplate", "uncertainty_temperature_baseplate_counts")
        self._assert_line_counts_uncertainty_variable_uint16(ds, "u_Tc_ch", "uncertainty_temperature_coolerhousing_counts")
        self._assert_line_counts_uncertainty_variable_uint16(ds, "u_Tc_elec", "uncertainty_temperature_electronics_counts")
        self._assert_line_counts_uncertainty_variable_uint16(ds, "u_Tc_fsr", "uncertainty_temperature_first_stage_radiator_counts")
        self._assert_line_counts_uncertainty_variable_uint16(ds, "u_Tc_fwh", "uncertainty_temperature_filter_wheel_housing_counts")
        self._assert_line_counts_uncertainty_variable_uint16(ds, "u_Tc_fwm", "uncertainty_temperature_filter_wheel_monitor_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_icct", "uncertainty_temperature_internal_cold_calibration_target_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_iwct", "uncertainty_temperature_internal_warm_calibration_target_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_patch_exp", "uncertainty_temperature_patch_expanded_scale_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_patch_full", "uncertainty_temperature_patch_full_range_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_tlscp_prim", "uncertainty_temperature_telescope_primary_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_tlscp_sec", "uncertainty_temperature_telescope_secondary_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_tlscp_tert", "uncertainty_temperature_telescope_tertiary_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_scanmirror", "uncertainty_temperature_scanmirror_counts")
        self._assert_line_counts_uncertainty_variable_uint32(ds, "u_Tc_scanmotor", "uncertainty_temperature_scanmotor_counts")

        u_solar_za = ds.variables["u_sol_za"]
        self.assertEqual((7, 56), u_solar_za.shape)
        self.assertEqual(-999.0, u_solar_za.data[4, 4])
        self.assertEqual(-999.0, u_solar_za.encoding["_FillValue"])
        self.assertEqual(0.01, u_solar_za.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_solar_za.encoding["dtype"])
        self.assertEqual("uncertainty_solar_zenith_angle", u_solar_za.attrs["standard_name"])
        self.assertEqual("degree", u_solar_za.attrs["units"])

        u_sol_aa = ds.variables["u_sol_aa"]
        self.assertEqual((7, 56), u_sol_aa.shape)
        self.assertEqual(-999.0, u_sol_aa.data[5, 5])
        self.assertEqual(-999.0, u_sol_aa.encoding["_FillValue"])
        self.assertEqual(0.01, u_sol_aa.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_sol_aa.encoding["dtype"])
        self.assertEqual("uncertainty_solar_azimuth_angle", u_sol_aa.attrs["standard_name"])
        self.assertEqual("degree", u_sol_aa.attrs["units"])

        u_sat_za = ds.variables["u_sat_za"]
        self.assertEqual((7, 56), u_sat_za.shape)
        self.assertEqual(-999.0, u_sat_za.data[5, 5])
        self.assertEqual(-999.0, u_sat_za.encoding["_FillValue"])
        self.assertEqual(0.01, u_sat_za.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_sat_za.encoding["dtype"])
        self.assertEqual("uncertainty_satellite_zenith_angle", u_sat_za.attrs["standard_name"])
        self.assertEqual("degree", u_sat_za.attrs["units"])

        u_sat_aa = ds.variables["u_sat_aa"]
        self.assertEqual((7, 56), u_sat_aa.shape)
        self.assertEqual(-999.0, u_sat_aa.data[6, 6])
        self.assertEqual(-999.0, u_sat_aa.encoding["_FillValue"])
        self.assertEqual(0.01, u_sat_aa.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_sat_aa.encoding["dtype"])
        self.assertEqual("uncertainty_local_azimuth_angle", u_sat_aa.attrs["standard_name"])
        self.assertEqual("degree", u_sat_aa.attrs["units"])

        u_c_earth_chan_corr = ds.variables["u_c_earth_chan_corr"]
        self.assertEqual((19, 19), u_c_earth_chan_corr.shape)
        self.assertTrue(np.isnan(u_c_earth_chan_corr.data[11, 14]))
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), u_c_earth_chan_corr.encoding["_FillValue"])
        self.assertEqual(0.01, u_c_earth_chan_corr.encoding["scale_factor"])
        self.assertEqual(np.int16, u_c_earth_chan_corr.encoding["dtype"])
        self.assertEqual("u_c_earth channel correlations", u_c_earth_chan_corr.attrs["long_name"])

        u_c_space = ds.variables["u_c_space"]
        self.assertEqual((19, 337), u_c_space.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_c_space.data[(12, 15)])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_c_space.attrs["_FillValue"])
        self.assertEqual("C_s", u_c_space.attrs["parameter"])
        self.assertEqual("gaussian", u_c_space.attrs["pdf_shape"])
        self.assertEqual(0.005, u_c_space.attrs["scale_factor"])
        self.assertEqual("count", u_c_space.attrs["units"])
        self.assertEqual("u_c_space_chan_corr", u_c_space.attrs["ancilliary_variables"])

        u_c_space_chan_corr = ds.variables["u_c_space_chan_corr"]
        self.assertEqual((19, 19), u_c_space_chan_corr.shape)
        self.assertTrue(np.isnan(u_c_space_chan_corr.data[11, 14]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_c_space_chan_corr.encoding["_FillValue"])
        self.assertEqual(0.01, u_c_space_chan_corr.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_c_space_chan_corr.encoding["dtype"])
        self.assertEqual("u_c_space channel correlations", u_c_space_chan_corr.attrs["long_name"])

        u_earthshine = ds.variables["u_Earthshine"]
        self.assertEqual((7, 19), u_earthshine.shape)
        self.assertTrue(np.isnan(u_earthshine.data[3, 5]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_earthshine.encoding["_FillValue"])
        self.assertEqual(0.01, u_earthshine.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_earthshine.encoding["dtype"])

        u_o_Re = ds.variables["u_O_Re"]
        self.assertEqual((7, 19), u_o_Re.shape)
        self.assertTrue(np.isnan(u_o_Re.data[4, 6]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_o_Re.encoding["_FillValue"])
        self.assertEqual(0.01, u_o_Re.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_o_Re.encoding["dtype"])

        u_o_TIWCT = ds.variables["u_O_TIWCT"]
        self.assertEqual((7,), u_o_TIWCT.shape)
        self.assertTrue(np.isnan(u_o_TIWCT.data[5]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_o_TIWCT.encoding["_FillValue"])
        self.assertEqual(0.01, u_o_TIWCT.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_o_TIWCT.encoding["dtype"])

        u_o_TPRT = ds.variables["u_O_TPRT"]
        self.assertEqual((7,), u_o_TPRT.shape)
        self.assertEqual(65535, u_o_TPRT.data[6])
        self.assertEqual(65535, u_o_TPRT.attrs["_FillValue"])
        self.assertEqual("all", u_o_TPRT.attrs["channels_affected"])
        self.assertEqual("rectangle_absolute", u_o_TPRT.attrs["pixel_correlation_form"])
        self.assertEqual("pixel", u_o_TPRT.attrs["pixel_correlation_units"])
        self.assertEqual([-np.inf, np.inf], u_o_TPRT.attrs["pixel_correlation_scales"])
        self.assertEqual("rectangle_absolute", u_o_TPRT.attrs["scan_correlation_form"])
        self.assertEqual("line", u_o_TPRT.attrs["scan_correlation_units"])
        self.assertEqual([-np.inf, np.inf], u_o_TPRT.attrs["scan_correlation_scales"])
        self.assertEqual("rectangle_absolute", u_o_TPRT.attrs["image_correlation_form"])
        self.assertEqual("images", u_o_TPRT.attrs["image_correlation_units"])
        self.assertEqual([-np.inf, np.inf], u_o_TPRT.attrs["image_correlation_scales"])
        self.assertEqual("O_TPRT", u_o_TPRT.attrs["parameter"])
        self.assertEqual("gaussian", u_o_TPRT.attrs["pdf_shape"])
        self.assertEqual(0.01, u_o_TPRT.attrs["scale_factor"])
        self.assertEqual("O_TPRT", u_o_TPRT.attrs["short_name"])
        self.assertEqual("K", u_o_TPRT.attrs["units"])
        self.assertEqual("u_O_TPRT_chan_corr", u_o_TPRT.attrs["ancilliary_variables"])

        u_Rself = ds.variables["u_Rself"]
        self.assertEqual((7, 19), u_Rself.shape)
        self.assertTrue(np.isnan(u_Rself.data[0, 8]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_Rself.encoding["_FillValue"])
        self.assertEqual(0.01, u_Rself.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_Rself.encoding["dtype"])

        u_srf_calib = ds.variables["u_SRF_calib"]
        self.assertEqual((7, 19), u_srf_calib.shape)
        self.assertTrue(np.isnan(u_srf_calib.data[1, 9]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_srf_calib.encoding["_FillValue"])
        self.assertEqual(0.01, u_srf_calib.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_srf_calib.encoding["dtype"])

        u_d_prt = ds.variables["u_d_PRT"]
        self.assertEqual((5, 4), u_d_prt.shape)
        self.assertTrue(np.isnan(u_d_prt.data[2,0]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_d_prt.encoding["_FillValue"])
        self.assertEqual(0.01, u_d_prt.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_d_prt.encoding["dtype"])

        u_electronics = ds.variables["u_electronics"]
        self.assertEqual((7, 19), u_electronics.shape)
        self.assertTrue(np.isnan(u_electronics.data[3,10]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_electronics.encoding["_FillValue"])
        self.assertEqual(0.01, u_electronics.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_electronics.encoding["dtype"])

        u_nonlinearity = ds.variables["u_nonlinearity"]
        self.assertEqual((19, ), u_nonlinearity.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), u_nonlinearity.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_nonlinearity.encoding["_FillValue"])
        self.assertEqual(0.01, u_nonlinearity.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_nonlinearity.encoding["dtype"])

        u_periodic_noise = ds.variables["u_periodic_noise"]
        self.assertEqual((7, 19), u_periodic_noise.shape)
        self.assertTrue(np.isnan(u_periodic_noise.data[4,11]))
        self.assertEqual(DefaultData.get_default_fill_value(np.uint16), u_periodic_noise.encoding["_FillValue"])
        self.assertEqual(0.01, u_periodic_noise.encoding["scale_factor"])
        self.assertEqual(np.uint16, u_periodic_noise.encoding["dtype"])

        emissivity = ds.variables["emissivity"]
        self.assertEqual((), emissivity.shape)
        self.assertTrue(np.isnan(emissivity.data))
        self.assertTrue(np.isnan(emissivity.attrs["_FillValue"]))
        self.assertEqual("emissivity", emissivity.attrs["long_name"])
        self.assertEqual("1", emissivity.attrs["units"])

        temp_corr_slope = ds.variables["temp_corr_slope"]
        self.assertEqual((), temp_corr_slope.shape)
        self.assertTrue(np.isnan(temp_corr_slope.data))
        self.assertTrue(np.isnan(temp_corr_slope.attrs["_FillValue"]))
        self.assertEqual("Slope for effective temperature correction", temp_corr_slope.attrs["long_name"])
        self.assertEqual("1", temp_corr_slope.attrs["units"])

        temp_corr_offset = ds.variables["temp_corr_offset"]
        self.assertEqual((), temp_corr_offset.shape)
        self.assertTrue(np.isnan(temp_corr_offset.data))
        self.assertTrue(np.isnan(temp_corr_offset.attrs["_FillValue"]))
        self.assertEqual("Offset for effective temperature correction", temp_corr_offset.attrs["long_name"])
        self.assertEqual("1", temp_corr_offset.attrs["units"])

        scnlintime = ds.variables["scnlintime"]
        self.assertEqual((7,), scnlintime.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), scnlintime.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), scnlintime.attrs["_FillValue"])
        self.assertEqual("time", scnlintime.attrs["standard_name"])
        self.assertEqual("Scan line time of day", scnlintime.attrs["long_name"])
        self.assertEqual("hrs_scnlintime", scnlintime.attrs["orig_name"])
        self.assertEqual("ms", scnlintime.attrs["units"])

        scnlinf = ds.variables["scnlinf"]
        self.assertEqual((7,), scnlinf.shape)
        self.assertEqual(0, scnlinf.data[4])
        self.assertEqual("16384, 32768", scnlinf.attrs["flag_masks"])
        self.assertEqual("clock_drift_correction southbound_data", scnlinf.attrs["flag_meanings"])
        self.assertEqual("status_flag", scnlinf.attrs["standard_name"])
        self.assertEqual("scanline_bitfield", scnlinf.attrs["long_name"])

        scantype = ds.variables["scantype"]
        self.assertEqual((7,), scantype.shape)
        self.assertEqual(0, scantype.data[5])
        self.assertEqual("0, 1, 2, 3", scantype.attrs["flag_values"])
        self.assertEqual("earth_view space_view cold_bb_view main_bb_view", scantype.attrs["flag_meanings"])
        self.assertEqual("status_flag", scantype.attrs["standard_name"])
        self.assertEqual("scantype_bitfield", scantype.attrs["long_name"])

    def _assert_2d_temperature_variable(self, ds, name, long_name, orig_name):
        variable = ds.variables[name]
        self.assertEqual((7, 5), variable.shape)
        self.assertTrue(np.isnan(variable.data[6, 4]))
        self.assertTrue(np.isnan(variable.attrs["_FillValue"]))
        self.assertEqual(long_name, variable.attrs["long_name"])
        self.assertEqual(orig_name, variable.attrs["orig_name"])
        self.assertEqual("K", variable.attrs["units"])

    def _assert_3d_temperature_variable(self, ds, name, long_name, orig_name):
        variable = ds.variables[name]
        self.assertEqual((4, 7, 5), variable.shape)
        self.assertTrue(np.isnan(variable.data[2, 5, 3]))
        self.assertTrue(np.isnan(variable.attrs["_FillValue"]))
        self.assertEqual(long_name, variable.attrs["long_name"])
        self.assertEqual(orig_name, variable.attrs["orig_name"])
        self.assertEqual("K", variable.attrs["units"])

    def _assert_line_counts_variable(self, ds, name, standard_name):
        variable = self._assert_line_int32_variable(ds, name, standard_name)
        self.assertEqual("count", variable.attrs["units"])

    def _assert_line_int8_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int8), variable.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.int8), variable.attrs["_FillValue"])
        self._assert_name_attributes(variable, standard_name, long_name, orig_name)

        return variable

    def _assert_line_int16_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), variable.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.int16), variable.attrs["_FillValue"])
        self._assert_name_attributes(variable, standard_name, long_name, orig_name)

        return variable

    def _assert_line_int32_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), variable.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.int32), variable.attrs["_FillValue"])
        self._assert_name_attributes(variable, standard_name, long_name, orig_name)

        return variable

    def _assert_line_uint16_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), variable.data[4])
        self.assertEqual(65535, variable.encoding["_FillValue"])
        self.assertEqual(0.01, variable.encoding["scale_factor"])
        self.assertEqual(np.uint16, variable.encoding["dtype"])
        self._assert_name_attributes(variable, standard_name, long_name, orig_name)

        return variable

    def _assert_line_scaled_int16_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), variable.data[4])
        self.assertEqual(-32767, variable.encoding["_FillValue"])
        self.assertEqual(0.01, variable.encoding["scale_factor"])
        self.assertEqual(np.int16, variable.encoding["dtype"])
        self._assert_name_attributes(variable, standard_name, long_name, orig_name)

        return variable

    def _assert_name_attributes(self, variable, standard_name, long_name, orig_name):
        if standard_name is not None:
            self.assertEqual(standard_name, variable.attrs["standard_name"])
        if long_name is not None:
            self.assertEqual(long_name, variable.attrs["long_name"])
        if orig_name is not None:
            self.assertEqual(orig_name, variable.attrs["orig_name"])

    def _assert_line_uint32_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        self.assertEqual(DefaultData.get_default_fill_value(np.float32), variable.data[4])
        self.assertEqual(DefaultData.get_default_fill_value(np.uint32), variable.encoding["_FillValue"])
        self.assertEqual(0.01, variable.encoding["scale_factor"])
        self.assertEqual(np.uint32, variable.encoding["dtype"])
        self._assert_name_attributes(variable, standard_name, long_name, orig_name)

        return variable

    def _assert_line_temperature_variable(self, ds, name, long_name, orig_name=None, fill_value=None):
        variable = self._assert_line_float_variable(ds, name, long_name=long_name, orig_name=orig_name, fill_value=fill_value)
        self.assertEqual("K", variable.attrs["units"])

    def _assert_line_counts_uncertainty_variable_uint16(self, ds, name, standard_name):
        variable = self._assert_line_uint16_variable(ds, name, standard_name=standard_name)
        self.assertEqual("count", variable.attrs["units"])

    def _assert_line_counts_uncertainty_variable_uint32(self, ds, name, standard_name):
        variable = self._assert_line_uint32_variable(ds, name, standard_name=standard_name)
        self.assertEqual("count", variable.attrs["units"])

    def _assert_line_angle_variable(self, ds, name, long_name=None, orig_name=None, fill_value=None):
        variable = self._assert_line_float_variable(ds, name, long_name=long_name, orig_name=orig_name, fill_value=fill_value)
        self.assertEqual("degree", variable.attrs["units"])

    def _assert_line_float_variable(self, ds, name, standard_name=None, long_name=None, orig_name=None, fill_value=None):
        variable = ds.variables[name]
        self.assertEqual((7,), variable.shape)
        if fill_value is None:
            self.assertEqual(DefaultData.get_default_fill_value(np.float32), variable.data[4])
            self.assertEqual(DefaultData.get_default_fill_value(np.float32), variable.attrs["_FillValue"])
        elif np.isnan(fill_value):
            self.assertTrue(np.isnan(variable.data[4]))
            self.assertTrue(np.isnan(variable.attrs["_FillValue"]))
        else:
            self.assertEqual(fill_value, variable.data[4])
            self.assertEqual(fill_value, variable.attrs["_FillValue"])

        self._assert_name_attributes(variable, standard_name, long_name, orig_name)
        return variable

