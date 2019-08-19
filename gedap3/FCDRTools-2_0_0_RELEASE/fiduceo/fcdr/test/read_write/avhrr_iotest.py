import datetime
import os
import tempfile
import unittest

import numpy as np
import xarray as xr

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter

N_PRT = 3
PRODUCT_WIDTH = 409
PRODUCT_HEIGHT = 13198
SRF_SIZE = 409
NUM_CHANNELS = 6

EXPECTED_CHUNKING = (1280, 409)


class AvhrrIoTest(unittest.TestCase):
    temp_dir = None
    target_path = None

    def setUp(self):
        self.temp_dir = tempfile.gettempdir()
        self.target_path = None

    def tearDown(self):
        if self.target_path is not None:
            os.remove(self.target_path)

    def test_write_easy(self):
        avhrr_easy = self.create_easy_dataset()

        start = datetime.datetime(2016, 11, 22, 13, 24, 52)
        end = datetime.datetime(2016, 11, 22, 14, 25, 53)
        file_name = FCDRWriter.create_file_name_FCDR_easy("AVHRR", "NOAA18", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.datetime.now()

        FCDRWriter.write(avhrr_easy, self.target_path)

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        print("AVHRR EASY write time: " + str(elapsed_time.seconds) + "." + str(round(elapsed_time.microseconds / 1000)))

        self.assertTrue(os.path.isfile(self.target_path))

        # target_data = xr.open_dataset(self.target_path, chunks=128)
        target_data = xr.open_dataset(self.target_path)
        try:
            self.assert_geolocation_variables(target_data)
            self.assert_global_flags(target_data)
            self.assert_sensor_variables(target_data)

            variable = target_data["relative_azimuth_angle"]
            self.assertAlmostEqual(0.11, variable.data[11, 11], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_independent_Ch1"]
            self.assertAlmostEqual(0.182, variable.data[14, 14], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_independent_Ch2"]
            self.assertAlmostEqual(0.21, variable.data[15, 15], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_independent_Ch3a"]
            self.assertAlmostEqual(0.24, variable.data[16, 16], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_independent_Ch4"]
            self.assertAlmostEqual(0.272, variable.data[17, 17], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_independent_Ch5"]
            self.assertAlmostEqual(0.306, variable.data[18, 18], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_structured_Ch1"]
            self.assertAlmostEqual(0.34, variable.data[19, 19], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_structured_Ch2"]
            self.assertAlmostEqual(0.38, variable.data[20, 20], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_structured_Ch3a"]
            self.assertAlmostEqual(0.42, variable.data[21, 21], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_structured_Ch4"]
            self.assertAlmostEqual(0.462, variable.data[22, 22], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_structured_Ch5"]
            self.assertAlmostEqual(0.506, variable.data[23, 23], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_common_Ch1"]
            self.assertAlmostEqual(0.072, variable.data[24, 24], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_common_Ch2"]
            self.assertAlmostEqual(0.1, variable.data[25, 25], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_common_Ch3a"]
            self.assertAlmostEqual(0.13, variable.data[26, 26], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_common_Ch4"]
            self.assertAlmostEqual(0.162, variable.data[27, 27], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_common_Ch5"]
            self.assertAlmostEqual(0.196, variable.data[28, 28], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            self.assert_srf_data(target_data)
        finally:
            target_data.close()

    def test_write_full(self):
        avhrr_full = self.create_full_dataset()

        start = datetime.datetime(2016, 11, 22, 14, 25, 53)
        end = datetime.datetime(2016, 11, 22, 15, 26, 54)
        file_name = FCDRWriter.create_file_name_FCDR_full("AVHRR", "NOAA19", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.datetime.now()

        FCDRWriter.write(avhrr_full, self.target_path)

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        print("AVHRR FULL write time: " + str(elapsed_time.seconds) + "." + str(round(elapsed_time.microseconds / 1000)))

        self.assertTrue(os.path.isfile(self.target_path))

        target_data = xr.open_dataset(self.target_path)
        try:
            self.assert_geolocation_variables(target_data)
            self.assert_global_flags(target_data)
            self.assert_sensor_variables(target_data)

            variable = target_data["u_latitude"]
            self.assertAlmostEqual(3.36, variable.data[24, 24], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_longitude"]
            self.assertAlmostEqual(3.75, variable.data[25, 25], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_time"]
            self.assertAlmostEqual(0.16, variable.data[256], 8)

            variable = target_data["u_satellite_azimuth_angle"]
            self.assertAlmostEqual(4.32, variable.data[27, 27], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_satellite_zenith_angle"]
            self.assertAlmostEqual(4.76, variable.data[28, 28], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_solar_azimuth_angle"]
            self.assertAlmostEqual(5.22, variable.data[29, 29], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_solar_zenith_angle"]
            self.assertAlmostEqual(5.7, variable.data[30, 30], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["PRT_C"]
            self.assertEqual(0, variable.data[31, 0])
            self.assertEqual((13198, 3), variable.encoding["chunksizes"])

            variable = target_data["u_prt"]
            self.assertAlmostEqual(0.2, variable.data[32, 1], 8)
            self.assertEqual((13198, 3), variable.encoding["chunksizes"])

            variable = target_data["R_ICT"]
            self.assertAlmostEqual(0.42, variable.data[33, 2], 7)
            self.assertEqual((13198, 3), variable.encoding["chunksizes"])

            variable = target_data["T_instr"]
            self.assertAlmostEqual(0.17, variable.data[257], 8)

            variable = target_data["Ch1_Csp"]
            self.assertEqual(68, variable.data[34, 34])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch2_Csp"]
            self.assertEqual(105, variable.data[35, 35])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3a_Csp"]
            self.assertEqual(144, variable.data[36, 36])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_Csp"]
            self.assertEqual(185, variable.data[37, 37])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_Csp"]
            self.assertEqual(228, variable.data[38, 38])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_Csp"]
            self.assertEqual(273, variable.data[39, 39])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_Cict"]
            self.assertEqual(320, variable.data[40, 40])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_Cict"]
            self.assertEqual(369, variable.data[41, 41])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_Cict"]
            self.assertEqual(420, variable.data[42, 42])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch1_Ce"]
            self.assertEqual(473, variable.data[43, 43])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch2_Ce"]
            self.assertEqual(528, variable.data[44, 44])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3a_Ce"]
            self.assertEqual(585, variable.data[45, 45])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_Ce"]
            self.assertEqual(644, variable.data[46, 46])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_Ce"]
            self.assertEqual(705, variable.data[47, 47])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_Ce"]
            self.assertEqual(768, variable.data[48, 48])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch1_u_Csp"]
            self.assertAlmostEqual(8.33, variable.data[49, 49], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch2_u_Csp"]
            self.assertAlmostEqual(9.0, variable.data[50, 50], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3a_u_Csp"]
            self.assertAlmostEqual(9.69, variable.data[51, 51], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_u_Csp"]
            self.assertAlmostEqual(10.4, variable.data[52, 52], 5)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_u_Csp"]
            self.assertAlmostEqual(11.13, variable.data[53, 53], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_u_Csp"]
            self.assertAlmostEqual(11.88, variable.data[54, 54], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_u_Cict"]
            self.assertAlmostEqual(12.65, variable.data[55, 55], 5)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_u_Cict"]
            self.assertAlmostEqual(13.44, variable.data[56, 56], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_u_Cict"]
            self.assertAlmostEqual(14.25, variable.data[57, 57], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch1_u_Ce"]
            self.assertAlmostEqual(15.08, variable.data[58, 58], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch2_u_Ce"]
            self.assertAlmostEqual(15.93, variable.data[59, 59], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3a_u_Ce"]
            self.assertAlmostEqual(16.8, variable.data[60, 60], 5)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_u_Ce"]
            self.assertAlmostEqual(17.69, variable.data[61, 61], 5)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_u_Ce"]
            self.assertAlmostEqual(18.6, variable.data[62, 62], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_u_Ce"]
            self.assertAlmostEqual(19.53, variable.data[63, 63], 5)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch1_u_Refl"]
            self.assertAlmostEqual(0.16384, variable.data[64, 64], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch2_u_Refl"]
            self.assertAlmostEqual(-0.17688, variable.data[65, 65], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3a_u_Refl"]
            self.assertAlmostEqual(0.15776, variable.data[66, 66], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_u_Bt"]
            self.assertAlmostEqual(23.45, variable.data[67, 67], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_u_Bt"]
            self.assertAlmostEqual(24.48, variable.data[68, 68], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_u_Bt"]
            self.assertAlmostEqual(25.53, variable.data[69, 69], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_ur_Bt"]
            self.assertAlmostEqual(26.6, variable.data[70, 70], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_ur_Bt"]
            self.assertAlmostEqual(27.69, variable.data[71, 71], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_ur_Bt"]
            self.assertAlmostEqual(28.8, variable.data[72, 72], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch3b_us_Bt"]
            self.assertAlmostEqual(29.93, variable.data[73, 73], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch4_us_Bt"]
            self.assertAlmostEqual(31.08, variable.data[74, 74], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["Ch5_us_Bt"]
            self.assertAlmostEqual(32.25, variable.data[75, 75], 6)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])
        finally:
            target_data.close()

    def create_easy_dataset(self):
        avhrr_easy = FCDRWriter.createTemplateEasy("AVHRR", PRODUCT_HEIGHT)
        self.add_global_attributes(avhrr_easy)
        self.add_geolocation_data(avhrr_easy)
        self.add_global_flags(avhrr_easy)
        self.add_sensor_data(avhrr_easy)
        self.add_srf_data(avhrr_easy)

        for x in range(0, PRODUCT_WIDTH):
            avhrr_easy["u_independent_Ch1"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.013
            avhrr_easy["u_independent_Ch2"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.014
            avhrr_easy["u_independent_Ch3a"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.015
            avhrr_easy["u_independent_Ch4"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.016
            avhrr_easy["u_independent_Ch5"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.017
            avhrr_easy["u_structured_Ch1"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.018
            avhrr_easy["u_structured_Ch2"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.019
            avhrr_easy["u_structured_Ch3a"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.020
            avhrr_easy["u_structured_Ch4"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.021
            avhrr_easy["u_structured_Ch5"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.022
            avhrr_easy["u_common_Ch1"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.003
            avhrr_easy["u_common_Ch2"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.004
            avhrr_easy["u_common_Ch3a"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.005
            avhrr_easy["u_common_Ch4"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.006
            avhrr_easy["u_common_Ch5"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x * 0.007

        for x in range(0, 6):
            avhrr_easy["quality_channel_bitmask"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int8) * x

        avhrr_easy["quality_scanline_bitmask"].data[:] = np.ones(PRODUCT_HEIGHT, np.int8)

        return avhrr_easy

    def create_full_dataset(self):
        avhrr_full = FCDRWriter.createTemplateFull("AVHRR", PRODUCT_HEIGHT)
        self.add_global_attributes(avhrr_full)
        self.add_geolocation_data(avhrr_full)
        self.add_global_flags(avhrr_full)
        self.add_sensor_data(avhrr_full)

        for x in range(0, PRODUCT_WIDTH):
            avhrr_full["u_latitude"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.14
            avhrr_full["u_longitude"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.15
            avhrr_full["u_satellite_azimuth_angle"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.16
            avhrr_full["u_satellite_zenith_angle"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.17
            avhrr_full["u_solar_azimuth_angle"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.18
            avhrr_full["u_solar_zenith_angle"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.19
            avhrr_full["Ch1_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 2
            avhrr_full["Ch2_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 3
            avhrr_full["Ch3a_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 4
            avhrr_full["Ch3b_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 5
            avhrr_full["Ch4_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 6
            avhrr_full["Ch5_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 7
            avhrr_full["Ch3b_Cict"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 8
            avhrr_full["Ch4_Cict"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 9
            avhrr_full["Ch5_Cict"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 10
            avhrr_full["Ch1_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 11
            avhrr_full["Ch2_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 12
            avhrr_full["Ch3a_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 13
            avhrr_full["Ch3b_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 14
            avhrr_full["Ch4_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 15
            avhrr_full["Ch5_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int32) * x * 16
            avhrr_full["Ch1_u_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.17
            avhrr_full["Ch2_u_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.18
            avhrr_full["Ch3a_u_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.19
            avhrr_full["Ch3b_u_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.2
            avhrr_full["Ch4_u_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.21
            avhrr_full["Ch5_u_Csp"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.22
            avhrr_full["Ch3b_u_Cict"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.23
            avhrr_full["Ch4_u_Cict"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.24
            avhrr_full["Ch5_u_Cict"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.25
            avhrr_full["Ch1_u_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.26
            avhrr_full["Ch2_u_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.27
            avhrr_full["Ch3a_u_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.28
            avhrr_full["Ch3b_u_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.29
            avhrr_full["Ch4_u_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.3
            avhrr_full["Ch5_u_Ce"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.31
            avhrr_full["Ch1_u_Refl"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.32
            avhrr_full["Ch2_u_Refl"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.33
            avhrr_full["Ch3a_u_Refl"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.34
            avhrr_full["Ch3b_u_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.35
            avhrr_full["Ch4_u_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.36
            avhrr_full["Ch5_u_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.37
            avhrr_full["Ch3b_ur_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.38
            avhrr_full["Ch4_ur_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.39
            avhrr_full["Ch5_ur_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.4
            avhrr_full["Ch3b_us_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.41
            avhrr_full["Ch4_us_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.42
            avhrr_full["Ch5_us_Bt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.43

        for x in range(0, N_PRT):
            avhrr_full["PRT_C"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int16) * x
            avhrr_full["u_prt"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.2
            avhrr_full["R_ICT"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.float32) * x * 0.21

        avhrr_full["u_time"].data[:] = np.ones(PRODUCT_HEIGHT, np.float64) * 0.16
        avhrr_full["T_instr"].data[:] = np.ones(PRODUCT_HEIGHT, np.float32) * 0.17

        for x in range(0, 6):
            avhrr_full["quality_channel_bitmask"].data[:, x] = np.ones(PRODUCT_HEIGHT, np.int8) * x

        avhrr_full["quality_scanline_bitmask"].data[:] = np.ones(PRODUCT_HEIGHT, np.int8)

        return avhrr_full

    def add_global_attributes(self, dataset):
        dataset.attrs["institution"] = "test"
        dataset.attrs["title"] = "sir"
        dataset.attrs["source"] = "invention"
        dataset.attrs["history"] = "new"
        dataset.attrs["references"] = "myself"
        dataset.attrs["comment"] = "should define a test version of this set"

    def add_geolocation_data(self, dataset):
        for x in range(0, PRODUCT_WIDTH):
            dataset["latitude"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.006
            dataset["longitude"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.007

    def add_global_flags(self, dataset):
        for x in range(0, PRODUCT_WIDTH):
            dataset["quality_pixel_bitmask"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int8) * x

    def add_sensor_data(self, dataset):
        for x in range(0, PRODUCT_WIDTH):
            dataset["Ch1"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.001
            dataset["Ch2"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.002
            dataset["Ch3a"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.003
            dataset["Ch4"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.004
            dataset["Ch5"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.005
            dataset["relative_azimuth_angle"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.010
            dataset["satellite_zenith_angle"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.011
            dataset["solar_zenith_angle"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int16) * x * 0.012
            dataset["data_quality_bitmask"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.uint8) * x

        dataset["Time"].data[:] = np.ones((PRODUCT_HEIGHT), np.float64)

    def add_srf_data(self, dataset):
        for x in range(0, SRF_SIZE):
            dataset["SRF_weights"].data[:, x] = np.ones((NUM_CHANNELS), np.float32) * x * 0.03
            dataset["SRF_wavelengths"].data[:, x] = np.ones((NUM_CHANNELS), np.float32) * x * 0.04

    def assert_geolocation_variables(self, target_data):
        variable = target_data["latitude"]
        self.assertAlmostEqual(0.0357066554, variable.data[6, 6], 8)
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])
        variable = target_data["longitude"]
        self.assertAlmostEqual(0.0494399853, variable.data[7, 7], 8)
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

    def assert_sensor_variables(self, target_data):
        variable = target_data["Ch1"]
        self.assertAlmostEqual(0.0, variable.data[0, 0], 8)
        self.assertEqual("toa_reflectance", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["Ch2"]
        self.assertAlmostEqual(0.002, variable.data[1, 1], 8)
        self.assertEqual("toa_reflectance", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["Ch3a"]
        self.assertAlmostEqual(0.006, variable.data[2, 2], 8)
        self.assertEqual("toa_reflectance", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["Ch4"]
        self.assertAlmostEqual(0.01, variable.data[3, 3], 8)
        self.assertEqual("toa_brightness_temperature", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["Ch5"]
        self.assertAlmostEqual(0.02, variable.data[4, 4], 8)
        self.assertEqual("toa_brightness_temperature", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["Time"]
        self.assertAlmostEqual(1.0, variable.data[5], 8)
        self.assertEqual("time", variable.attrs["standard_name"])
        self.assertEqual((13198,), variable.encoding["chunksizes"])

        variable = target_data["satellite_zenith_angle"]
        self.assertAlmostEqual(0.13, variable.data[12, 12], 8)
        self.assertEqual("sensor_zenith_angle", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["solar_zenith_angle"]
        self.assertAlmostEqual(0.16, variable.data[13, 13], 8)
        self.assertEqual("solar_zenith_angle", variable.attrs["standard_name"])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["data_quality_bitmask"]
        self.assertEqual(2, variable.data[9, 2])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["quality_channel_bitmask"]
        self.assertEqual(1, variable.data[8, 1])
        self.assertEqual((13198, 6), variable.encoding["chunksizes"])

        variable = target_data["quality_scanline_bitmask"]
        self.assertEqual(1, variable.data[10])
        self.assertEqual((13198,), variable.encoding["chunksizes"])

    def assert_global_flags(self, target_data):
        variable = target_data["quality_pixel_bitmask"]
        self.assertEqual(11, variable.data[9, 9])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

    def assert_srf_data(self, target_data):
        variable = target_data["SRF_weights"]
        self.assertAlmostEqual(0.300003, variable.data[0, 10], 8)

        variable = target_data["SRF_wavelengths"]
        self.assertAlmostEqual(0.44, variable.data[1, 11], 8)
