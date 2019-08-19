import datetime
import unittest

from fiduceo.common.test.assertions import Assertions
from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter


class FCDRWriterTest(unittest.TestCase):
    def testCreateTemplateEasy_AMSUB(self):
        ds = FCDRWriter.createTemplateEasy('AMSUB', 2561)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(19, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self._verify_amsub_specific_variables(ds)

        # easy FCDR variables
        self.assertIsNotNone(ds.variables["u_independent_btemps"])
        self.assertIsNotNone(ds.variables["u_structured_btemps"])

    def testCreateTemplateFull_AMSUB(self):
        ds = FCDRWriter.createTemplateFull('AMSUB', 2562)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(27, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self._verify_amsub_specific_variables(ds)

        # full FCDR variables
        self.assertIsNotNone(ds.variables["u_btemps"])
        self.assertIsNotNone(ds.variables["u_syst_btemps"])
        self.assertIsNotNone(ds.variables["u_random_btemps"])
        self.assertIsNotNone(ds.variables["u_instrtemp"])
        self.assertIsNotNone(ds.variables["u_latitude"])
        self.assertIsNotNone(ds.variables["u_longitude"])
        self.assertIsNotNone(ds.variables["u_satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["u_satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["u_solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["u_solar_zenith_angle"])

    def testCreateTemplateEasy_SSMT2(self):
        ds = FCDRWriter.createTemplateEasy('SSMT2', 722)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(15, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self.verify_SSMT2_specific_variables(ds)

        # easy FCDR variables
        self.assertIsNotNone(ds.variables["u_independent_tb"])
        self.assertIsNotNone(ds.variables["u_structured_tb"])

    def testCreateTemplateFull_SSMT2(self):
        ds = FCDRWriter.createTemplateFull('SSMT2', 722)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(21, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self.verify_SSMT2_specific_variables(ds)

        # easy FCDR variables
        self.assertIsNotNone(ds.variables["u_Temperature_misc_housekeeping"])
        self.assertIsNotNone(ds.variables["u_cold_counts"])
        self.assertIsNotNone(ds.variables["u_counts_to_tb_gain"])
        self.assertIsNotNone(ds.variables["u_counts_to_tb_offset"])
        self.assertIsNotNone(ds.variables["u_gain_control"])
        self.assertIsNotNone(ds.variables["u_tb"])
        self.assertIsNotNone(ds.variables["u_thermal_reference"])
        self.assertIsNotNone(ds.variables["u_warm_counts"])

    def testCreateTemplateEasy_AVHRR(self):
        ds = FCDRWriter.createTemplateEasy('AVHRR', 12198, srf_size=13, corr_dx=14, corr_dy=15)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(42, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self._verify_avhrr_specific_variables(ds)

        # easy FCDR variables
        self.assertIsNotNone(ds.variables["u_independent_Ch1"])
        self.assertIsNotNone(ds.variables["u_structured_Ch1"])
        self.assertIsNotNone(ds.variables["u_common_Ch1"])
        self.assertIsNotNone(ds.variables["u_independent_Ch2"])
        self.assertIsNotNone(ds.variables["u_structured_Ch2"])
        self.assertIsNotNone(ds.variables["u_common_Ch2"])
        self.assertIsNotNone(ds.variables["u_independent_Ch3a"])
        self.assertIsNotNone(ds.variables["u_structured_Ch3a"])
        self.assertIsNotNone(ds.variables["u_common_Ch3a"])
        self.assertIsNotNone(ds.variables["u_independent_Ch3b"])
        self.assertIsNotNone(ds.variables["u_structured_Ch3b"])
        self.assertIsNotNone(ds.variables["u_common_Ch3b"])
        self.assertIsNotNone(ds.variables["u_independent_Ch4"])
        self.assertIsNotNone(ds.variables["u_structured_Ch4"])
        self.assertIsNotNone(ds.variables["u_common_Ch4"])
        self.assertIsNotNone(ds.variables["u_independent_Ch5"])
        self.assertIsNotNone(ds.variables["u_structured_Ch5"])
        self.assertIsNotNone(ds.variables["u_common_Ch5"])

        self.assertIsNotNone(ds.variables["channel_correlation_matrix_independent"])
        self.assertIsNotNone(ds.variables["channel_correlation_matrix_structured"])

        self.assertIsNotNone(ds.variables["cross_element_correlation_coefficients"])
        self.assertIsNotNone(ds.variables["cross_line_correlation_coefficients"])

    def testCreateTemplateFull_AVHRR(self):
        ds = FCDRWriter.createTemplateFull('AVHRR', 13667)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(73, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self._verify_avhrr_specific_variables(ds)

        # variables of full FCDR
        self.assertIsNotNone(ds.variables["u_latitude"])
        self.assertIsNotNone(ds.variables["u_longitude"])
        self.assertIsNotNone(ds.variables["u_time"])
        self.assertIsNotNone(ds.variables["u_satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["u_satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["u_solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["u_solar_zenith_angle"])

        self.assertIsNotNone(ds.variables["PRT_C"])
        self.assertIsNotNone(ds.variables["u_prt"])  # geolocation
        self._verify_geolocation_variables(ds)
        self.assertIsNotNone(ds.variables["R_ICT"])
        self.assertIsNotNone(ds.variables["T_instr"])

        self.assertIsNotNone(ds.variables["Ch1_Csp"])
        self.assertIsNotNone(ds.variables["Ch2_Csp"])
        self.assertIsNotNone(ds.variables["Ch3a_Csp"])
        self.assertIsNotNone(ds.variables["Ch3b_Csp"])
        self.assertIsNotNone(ds.variables["Ch4_Csp"])
        self.assertIsNotNone(ds.variables["Ch5_Csp"])

        self.assertIsNotNone(ds.variables["Ch3b_Cict"])
        self.assertIsNotNone(ds.variables["Ch4_Cict"])
        self.assertIsNotNone(ds.variables["Ch5_Cict"])

        self.assertIsNotNone(ds.variables["Ch1_Ce"])
        self.assertIsNotNone(ds.variables["Ch2_Ce"])
        self.assertIsNotNone(ds.variables["Ch3a_Ce"])
        self.assertIsNotNone(ds.variables["Ch3b_Ce"])
        self.assertIsNotNone(ds.variables["Ch4_Ce"])
        self.assertIsNotNone(ds.variables["Ch5_Ce"])

        self.assertIsNotNone(ds.variables["Ch1_u_Csp"])
        self.assertIsNotNone(ds.variables["Ch2_u_Csp"])
        self.assertIsNotNone(ds.variables["Ch3a_u_Csp"])
        self.assertIsNotNone(ds.variables["Ch3b_u_Csp"])
        self.assertIsNotNone(ds.variables["Ch4_u_Csp"])
        self.assertIsNotNone(ds.variables["Ch5_u_Csp"])

        self.assertIsNotNone(ds.variables["Ch3b_u_Cict"])
        self.assertIsNotNone(ds.variables["Ch4_u_Cict"])
        self.assertIsNotNone(ds.variables["Ch5_u_Cict"])

        self.assertIsNotNone(ds.variables["Ch1_u_Ce"])
        self.assertIsNotNone(ds.variables["Ch2_u_Ce"])
        self.assertIsNotNone(ds.variables["Ch3a_u_Ce"])
        self.assertIsNotNone(ds.variables["Ch3b_u_Ce"])
        self.assertIsNotNone(ds.variables["Ch4_u_Ce"])
        self.assertIsNotNone(ds.variables["Ch5_u_Ce"])

        self.assertIsNotNone(ds.variables["Ch1_u_Refl"])
        self.assertIsNotNone(ds.variables["Ch2_u_Refl"])
        self.assertIsNotNone(ds.variables["Ch3a_u_Refl"])

        self.assertIsNotNone(ds.variables["Ch3b_u_Bt"])
        self.assertIsNotNone(ds.variables["Ch4_u_Bt"])
        self.assertIsNotNone(ds.variables["Ch5_u_Bt"])

        self.assertIsNotNone(ds.variables["Ch3b_ur_Bt"])
        self.assertIsNotNone(ds.variables["Ch4_ur_Bt"])
        self.assertIsNotNone(ds.variables["Ch5_ur_Bt"])

        self.assertIsNotNone(ds.variables["Ch3b_us_Bt"])
        self.assertIsNotNone(ds.variables["Ch4_us_Bt"])
        self.assertIsNotNone(ds.variables["Ch5_us_Bt"])

    def testCreateTemplateEasy_HIRS2(self):
        ds = FCDRWriter.createTemplateEasy('HIRS2', 211)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(19, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self.assertIsNotNone(ds.variables["bt"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["scanline"])
        self.assertIsNotNone(ds.variables["time"])
        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["quality_scanline_bitmask"])
        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_wavelengths"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])
        self.assertIsNotNone(ds.variables["scanline_origl1b"])

        self._assert_HIRS_easy_uncertainties(ds)

    def testCreateTemplateFull_HIRS2(self):
        ds = FCDRWriter.createTemplateFull('HIRS2', 209)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(95, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        self._assert_HIRS_common_variables(ds)

        self.assertIsNotNone(ds.variables["c_earth"])
        self.assertIsNotNone(ds.variables["L_earth"])
        self.assertIsNotNone(ds.variables["u_L_earth_independent"])

        self.assertIsNotNone(ds.variables["navigation_status"])
        self.assertIsNotNone(ds.variables["platform_altitude"])
        self.assertIsNotNone(ds.variables["platform_pitch_angle"])
        self.assertIsNotNone(ds.variables["platform_roll_angle"])
        self.assertIsNotNone(ds.variables["platform_yaw_angle"])
        self.assertIsNotNone(ds.variables["quality_flags"])
        self.assertIsNotNone(ds.variables["scan_angles"])
        self.assertIsNotNone(ds.variables["l1b_scanline_number"])
        self.assertIsNotNone(ds.variables["scanline_position"])
        self.assertIsNotNone(ds.variables["l1b_second_original_calibration_coefficients"])
        self.assertIsNotNone(ds.variables["u_c_earth"])
        self.assertIsNotNone(ds.variables["u_c_earth_chan_corr"])
        self.assertIsNotNone(ds.variables["u_c_space"])
        self.assertIsNotNone(ds.variables["u_c_space_chan_corr"])
        self.assertIsNotNone(ds.variables["u_Earthshine"])
        self.assertIsNotNone(ds.variables["u_O_Re"])
        self.assertIsNotNone(ds.variables["u_O_TIWCT"])
        self.assertIsNotNone(ds.variables["u_O_TPRT"])
        self.assertIsNotNone(ds.variables["u_Rself"])
        self.assertIsNotNone(ds.variables["u_SRF_calib"])
        self.assertIsNotNone(ds.variables["u_d_PRT"])
        self.assertIsNotNone(ds.variables["u_electronics"])
        self.assertIsNotNone(ds.variables["u_nonlinearity"])

        self.assertIsNotNone(ds.variables["temp_corr_slope"])
        self.assertIsNotNone(ds.variables["temp_corr_offset"])
        self.assertIsNotNone(ds.variables["emissivity"])

    def _assert_HIRS_common_variables(self, ds):
        self.assertIsNotNone(ds.variables["bt"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["scanline"])
        self.assertIsNotNone(ds.variables["scnlintime"])
        self.assertIsNotNone(ds.variables["scantype"])
        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["quality_scanline_bitmask"])
        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_wavelengths"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])
        self.assertIsNotNone(ds.variables["scanline_origl1b"])

    def testCreateTemplateEasy_HIRS3(self):
        ds = FCDRWriter.createTemplateEasy('HIRS3', 211)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(22, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self.assertIsNotNone(ds.variables["bt"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["scanline"])
        self.assertIsNotNone(ds.variables["time"])
        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["quality_scanline_bitmask"])
        self.assertIsNotNone(ds.variables["quality_channel_bitmask"])
        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_wavelengths"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])

        self._assert_HIRS_easy_uncertainties(ds)

    def testCreateTemplateFull_HIRS3(self):
        ds = FCDRWriter.createTemplateFull('HIRS3', 209)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(98, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        self._assert_HIRS_common_variables(ds)
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])
        self.assertIsNotNone(ds.variables["scanline_origl1b"])

        self.assertIsNotNone(ds.variables["data_quality_bitmask"])

        self.assertIsNotNone(ds.variables["c_earth"])
        self.assertIsNotNone(ds.variables["L_earth"])
        self.assertIsNotNone(ds.variables["u_L_earth_independent"])

        self.assertIsNotNone(ds.variables["navigation_status"])
        self.assertIsNotNone(ds.variables["platform_altitude"])
        self.assertIsNotNone(ds.variables["platform_pitch_angle"])
        self.assertIsNotNone(ds.variables["platform_roll_angle"])
        self.assertIsNotNone(ds.variables["platform_yaw_angle"])
        self.assertIsNotNone(ds.variables["quality_flags"])
        self.assertIsNotNone(ds.variables["scan_angles"])
        self.assertIsNotNone(ds.variables["l1b_scanline_number"])
        self.assertIsNotNone(ds.variables["scanline_position"])
        self.assertIsNotNone(ds.variables["l1b_second_original_calibration_coefficients"])
        self.assertIsNotNone(ds.variables["u_c_earth"])
        self.assertIsNotNone(ds.variables["u_c_earth_chan_corr"])
        self.assertIsNotNone(ds.variables["u_c_space"])
        self.assertIsNotNone(ds.variables["u_c_space_chan_corr"])
        self.assertIsNotNone(ds.variables["u_Earthshine"])
        self.assertIsNotNone(ds.variables["u_O_Re"])
        self.assertIsNotNone(ds.variables["u_O_TIWCT"])
        self.assertIsNotNone(ds.variables["u_O_TPRT"])
        self.assertIsNotNone(ds.variables["u_d_PRT"])
        self.assertIsNotNone(ds.variables["u_electronics"])
        self.assertIsNotNone(ds.variables["u_nonlinearity"])

        self.assertIsNotNone(ds.variables["temp_corr_slope"])
        self.assertIsNotNone(ds.variables["temp_corr_offset"])
        self.assertIsNotNone(ds.variables["emissivity"])
        self.assertIsNotNone(ds.variables["mnfrqualflags"])

    def testCreateTemplateEasy_HIRS4(self):
        ds = FCDRWriter.createTemplateEasy('HIRS4', 211)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(22, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        # sensor specific
        self.assertIsNotNone(ds.variables["bt"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["scanline"])
        self.assertIsNotNone(ds.variables["time"])
        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["quality_scanline_bitmask"])
        self.assertIsNotNone(ds.variables["quality_channel_bitmask"])
        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_wavelengths"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])
        self.assertIsNotNone(ds.variables["scanline_origl1b"])

        self._assert_HIRS_easy_uncertainties(ds)

    def _assert_HIRS_easy_uncertainties(self, ds):
        self.assertIsNotNone(ds.variables["u_independent"])
        self.assertIsNotNone(ds.variables["u_structured"])
        self.assertIsNotNone(ds.variables["u_common"])
        self.assertIsNotNone(ds.variables["channel_correlation_matrix_independent"])
        self.assertIsNotNone(ds.variables["channel_correlation_matrix_structured"])

    def testCreateTemplateFull_HIRS4(self):
        ds = FCDRWriter.createTemplateFull('HIRS4', 209)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(98, len(ds.data_vars))

        # geolocation + flags
        self._verify_geolocation_variables(ds)
        self._verify_quality_flags(ds)

        self._assert_HIRS_common_variables(ds)
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])

        self.assertIsNotNone(ds.variables["c_earth"])
        self.assertIsNotNone(ds.variables["L_earth"])

        self.assertIsNotNone(ds.variables["u_L_earth_independent"])

        self.assertIsNotNone(ds.variables["navigation_status"])
        self.assertIsNotNone(ds.variables["platform_altitude"])
        self.assertIsNotNone(ds.variables["platform_pitch_angle"])
        self.assertIsNotNone(ds.variables["platform_roll_angle"])
        self.assertIsNotNone(ds.variables["platform_yaw_angle"])
        self.assertIsNotNone(ds.variables["quality_flags"])
        self.assertIsNotNone(ds.variables["scan_angles"])
        self.assertIsNotNone(ds.variables["l1b_scanline_number"])
        self.assertIsNotNone(ds.variables["scanline_position"])
        self.assertIsNotNone(ds.variables["l1b_second_original_calibration_coefficients"])
        self.assertIsNotNone(ds.variables["u_c_earth"])
        self.assertIsNotNone(ds.variables["u_c_earth_chan_corr"])
        self.assertIsNotNone(ds.variables["u_c_space"])
        self.assertIsNotNone(ds.variables["u_c_space_chan_corr"])
        self.assertIsNotNone(ds.variables["u_Earthshine"])
        self.assertIsNotNone(ds.variables["u_O_Re"])
        self.assertIsNotNone(ds.variables["u_O_TIWCT"])
        self.assertIsNotNone(ds.variables["u_O_TPRT"])
        self.assertIsNotNone(ds.variables["u_nonlinearity"])

        self.assertIsNotNone(ds.variables["temp_corr_slope"])
        self.assertIsNotNone(ds.variables["temp_corr_offset"])
        self.assertIsNotNone(ds.variables["emissivity"])
        self.assertIsNotNone(ds.variables["mnfrqualflags"])

    def testCreateTemplateEasy_MVIRI(self):
        ds = FCDRWriter.createTemplateEasy('MVIRI', 5000)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self._verify_quality_flags(ds)

        self.assertEqual(40, len(ds.data_vars))

        # sensor specific
        self.assertIsNotNone(ds.variables["time"])
        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_frequencies"])
        self.assertIsNotNone(ds.variables["covariance_spectral_response_function_vis"])
        self.assertIsNotNone(ds.variables["u_spectral_response_function_ir"])
        self.assertIsNotNone(ds.variables["u_spectral_response_function_wv"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["count_ir"])
        self.assertIsNotNone(ds.variables["count_wv"])
        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["distance_sun_earth"])
        self.assertIsNotNone(ds.variables["solar_irradiance_vis"])
        self.assertIsNotNone(ds.variables["u_solar_irradiance_vis"])
        self.assertIsNotNone(ds.variables["a_ir"])
        self.assertIsNotNone(ds.variables["b_ir"])
        self.assertIsNotNone(ds.variables["u_a_ir"])
        self.assertIsNotNone(ds.variables["u_b_ir"])
        self.assertIsNotNone(ds.variables["a_wv"])
        self.assertIsNotNone(ds.variables["b_wv"])
        self.assertIsNotNone(ds.variables["u_a_wv"])
        self.assertIsNotNone(ds.variables["u_b_wv"])
        self.assertIsNotNone(ds.variables["bt_a_ir"])
        self.assertIsNotNone(ds.variables["bt_b_ir"])
        self.assertIsNotNone(ds.variables["bt_a_wv"])
        self.assertIsNotNone(ds.variables["bt_b_wv"])
        self.assertIsNotNone(ds.variables["years_since_launch"])

        # easy FCDR uncertainties
        self.assertIsNotNone(ds.variables["toa_bidirectional_reflectance_vis"])
        self.assertIsNotNone(ds.variables["u_independent_toa_bidirectional_reflectance"])
        self.assertIsNotNone(ds.variables["u_structured_toa_bidirectional_reflectance"])
        self.assertIsNotNone(ds.variables["u_common_toa_bidirectional_reflectance"])
        self.assertIsNotNone(ds.variables["sub_satellite_latitude_start"])
        self.assertIsNotNone(ds.variables["sub_satellite_longitude_start"])
        self.assertIsNotNone(ds.variables["sub_satellite_latitude_end"])
        self.assertIsNotNone(ds.variables["sub_satellite_longitude_end"])
        self.assertIsNotNone(ds.variables["channel_correlation_matrix_independent"])
        self.assertIsNotNone(ds.variables["channel_correlation_matrix_structured"])

    def testCreateTemplateFull_MVIRI(self):
        ds = FCDRWriter.createTemplateFull('MVIRI', 5000)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)

        self.assertEqual(58, len(ds.data_vars))

        self._verify_quality_flags(ds)

        # sensor specific
        self.assertIsNotNone(ds.variables["count_vis"])
        self.assertIsNotNone(ds.variables["time"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["count_ir"])
        self.assertIsNotNone(ds.variables["count_wv"])
        self.assertIsNotNone(ds.variables["count_vis"])
        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_frequencies"])
        self.assertIsNotNone(ds.variables["covariance_spectral_response_function_vis"])
        self.assertIsNotNone(ds.variables["u_spectral_response_function_ir"])
        self.assertIsNotNone(ds.variables["u_spectral_response_function_wv"])
        self.assertIsNotNone(ds.variables["a0_vis"])
        self.assertIsNotNone(ds.variables["a1_vis"])
        self.assertIsNotNone(ds.variables["a2_vis"])
        self.assertIsNotNone(ds.variables["solar_irradiance_vis"])
        self.assertIsNotNone(ds.variables["u_solar_irradiance_vis"])
        self.assertIsNotNone(ds.variables["distance_sun_earth"])
        self.assertIsNotNone(ds.variables["a_ir"])
        self.assertIsNotNone(ds.variables["b_ir"])
        self.assertIsNotNone(ds.variables["u_a_ir"])
        self.assertIsNotNone(ds.variables["u_b_ir"])
        self.assertIsNotNone(ds.variables["a_wv"])
        self.assertIsNotNone(ds.variables["b_wv"])
        self.assertIsNotNone(ds.variables["u_a_wv"])
        self.assertIsNotNone(ds.variables["u_b_wv"])
        self.assertIsNotNone(ds.variables["bt_a_ir"])
        self.assertIsNotNone(ds.variables["bt_b_ir"])
        self.assertIsNotNone(ds.variables["bt_a_wv"])
        self.assertIsNotNone(ds.variables["bt_b_wv"])
        self.assertIsNotNone(ds.variables["years_since_launch"])

        # full FCDR uncertainties
        self.assertIsNotNone(ds.variables["u_latitude"])
        self.assertIsNotNone(ds.variables["u_longitude"])
        self.assertIsNotNone(ds.variables["u_time"])
        self.assertIsNotNone(ds.variables["u_a0_vis"])
        self.assertIsNotNone(ds.variables["u_a1_vis"])
        self.assertIsNotNone(ds.variables["u_a2_vis"])
        self.assertIsNotNone(ds.variables["u_zero_vis"])
        self.assertIsNotNone(ds.variables["covariance_a_vis"])
        self.assertIsNotNone(ds.variables["u_electronics_counts_vis"])
        self.assertIsNotNone(ds.variables["u_digitization_counts_vis"])
        self.assertIsNotNone(ds.variables["allan_deviation_counts_space_vis"])
        self.assertIsNotNone(ds.variables["u_mean_counts_space_vis"])
        self.assertIsNotNone(ds.variables["u_solar_zenith_angle"])
        self.assertIsNotNone(ds.variables["u_solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["u_satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["u_satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["effect_correlation_matrix"])

        # virtual variables
        self.assertIsNotNone(ds.variables["sensitivity_solar_irradiance_vis"])
        self.assertIsNotNone(ds.variables["sensitivity_count_vis"])
        self.assertIsNotNone(ds.variables["sensitivity_count_space"])
        self.assertIsNotNone(ds.variables["sensitivity_a0_vis"])
        self.assertIsNotNone(ds.variables["sensitivity_a1_vis"])
        self.assertIsNotNone(ds.variables["sensitivity_a2_vis"])

    def testCreate_MVIRI_STATIC(self):
        ds = FCDRWriter.createTemplateEasy('MVIRI_STATIC', 5000)
        self.assertIsNotNone(ds.variables["latitude_vis"])
        self.assertIsNotNone(ds.variables["longitude_vis"])
        self.assertIsNotNone(ds.variables["latitude_ir_wv"])
        self.assertIsNotNone(ds.variables["longitude_ir_wv"])

    def test_create_file_name_FCDR_easy(self):
        start = datetime.datetime(2015, 8, 23, 14, 24, 52)
        end = datetime.datetime(2015, 8, 23, 15, 25, 53)
        self.assertEqual("FIDUCEO_FCDR_L1C_MVIRI_MET7-0.00_20150823142452_20150823152553_EASY_v02.3_fv2.0.0.nc", FCDRWriter.create_file_name_FCDR_easy("MVIRI", "MET7-0.00", start, end, "02.3"))

        start = datetime.datetime(2014, 7, 22, 13, 23, 51)
        end = datetime.datetime(2014, 7, 22, 14, 24, 52)
        self.assertEqual("FIDUCEO_FCDR_L1C_HIRS3_NOAA15_20140722132351_20140722142452_EASY_v03.4_fv2.0.0.nc", FCDRWriter.create_file_name_FCDR_easy("HIRS3", "NOAA15", start, end, "03.4"))

        start = datetime.datetime(2013, 6, 21, 12, 23, 50)
        end = datetime.datetime(2013, 6, 21, 13, 23, 51)
        self.assertEqual("FIDUCEO_FCDR_L1C_HIRS4_METOPA_20130621122350_20130621132351_EASY_v04.5_fv2.0.0.nc", FCDRWriter.create_file_name_FCDR_easy("HIRS4", "METOPA", start, end, "04.5"))

        start = datetime.datetime(2012, 5, 20, 11, 22, 49)
        end = datetime.datetime(2012, 5, 20, 12, 22, 50)
        self.assertEqual("FIDUCEO_FCDR_L1C_AMSUB_NOAA17_20120520112249_20120520122250_EASY_v05.6_fv2.0.0.nc", FCDRWriter.create_file_name_FCDR_easy("AMSUB", "NOAA17", start, end, "05.6"))

    def test_create_file_name_FCDR_full(self):
        start = datetime.datetime(2015, 8, 7, 14, 24, 52)
        end = datetime.datetime(2015, 8, 7, 15, 25, 53)
        self.assertEqual("FIDUCEO_FCDR_L1C_MVIRI_MET7-0.00_20150807142452_20150807152553_FULL_v02.3_fv2.0.0.nc", FCDRWriter.create_file_name_FCDR_full("MVIRI", "MET7-0.00", start, end, "02.3"))

        start = datetime.datetime(2014, 7, 21, 13, 23, 51)
        end = datetime.datetime(2014, 7, 21, 14, 24, 52)
        self.assertEqual("FIDUCEO_FCDR_L1C_HIRS3_NOAA15_20140721132351_20140721142452_FULL_v03.4_fv2.0.0.nc", FCDRWriter.create_file_name_FCDR_full("HIRS3", "NOAA15", start, end, "03.4"))

    def _verify_geolocation_variables(self, ds):
        self.assertIsNotNone(ds.variables["latitude"])
        self.assertIsNotNone(ds.variables["longitude"])

    def _verify_quality_flags(self, ds):
        self.assertIsNotNone(ds.variables["quality_pixel_bitmask"])

    def _verify_geometry_variables(self, ds):
        self.assertIsNotNone(ds.variables["satellite_azimuth_angle"])
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_azimuth_angle"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])

    def _verify_amsub_specific_variables(self, ds):
        self.assertIsNotNone(ds.variables["btemps"])
        self.assertIsNotNone(ds.variables["chanqual"])
        self.assertIsNotNone(ds.variables["instrtemp"])
        self.assertIsNotNone(ds.variables["qualind"])
        self.assertIsNotNone(ds.variables["scanqual"])
        self.assertIsNotNone(ds.variables["scnlin"])
        self.assertIsNotNone(ds.variables["scnlindy"])
        self.assertIsNotNone(ds.variables["scnlintime"])
        self.assertIsNotNone(ds.variables["scnlinyr"])
        # geometry
        self._verify_geometry_variables(ds)

    def verify_SSMT2_specific_variables(self, ds):
        self.assertIsNotNone(ds.variables["Temperature_misc_housekeeping"])
        self.assertIsNotNone(ds.variables["ancil_data"])
        self.assertIsNotNone(ds.variables["channel_quality_flag"])
        self.assertIsNotNone(ds.variables["cold_counts"])
        self.assertIsNotNone(ds.variables["counts_to_tb_gain"])
        self.assertIsNotNone(ds.variables["counts_to_tb_offset"])
        self.assertIsNotNone(ds.variables["gain_control"])
        self.assertIsNotNone(ds.variables["tb"])
        self.assertIsNotNone(ds.variables["thermal_reference"])
        self.assertIsNotNone(ds.variables["warm_counts"])

    def _verify_avhrr_specific_variables(self, ds):
        self.assertIsNotNone(ds.variables["Time"])
        # geometry
        self.assertIsNotNone(ds.variables["satellite_zenith_angle"])
        self.assertIsNotNone(ds.variables["solar_zenith_angle"])

        self.assertIsNotNone(ds.variables["Ch1"])
        self.assertIsNotNone(ds.variables["Ch2"])
        self.assertIsNotNone(ds.variables["Ch3a"])
        self.assertIsNotNone(ds.variables["Ch3b"])
        self.assertIsNotNone(ds.variables["Ch4"])
        self.assertIsNotNone(ds.variables["Ch5"])

        self.assertIsNotNone(ds.variables["data_quality_bitmask"])
        self.assertIsNotNone(ds.variables["quality_scanline_bitmask"])
        self.assertIsNotNone(ds.variables["quality_channel_bitmask"])

        self.assertIsNotNone(ds.variables["SRF_weights"])
        self.assertIsNotNone(ds.variables["SRF_wavelengths"])
        self.assertIsNotNone(ds.variables["scanline_map_to_origl1bfile"])
        self.assertIsNotNone(ds.variables["scanline_origl1b"])
