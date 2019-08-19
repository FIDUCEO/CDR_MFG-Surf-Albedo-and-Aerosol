import os
import tempfile
import unittest
from datetime import datetime

import numpy as np
import xarray as xr

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter

EXPECTED_CHUNKING = (500, 500)

PRODUCT_WIDTH = 5000
PRODUCT_HEIGHT = 5000


class MviriIoTest(unittest.TestCase):
    temp_dir = None
    target_path = None

    def setUp(self):
        self.temp_dir = tempfile.gettempdir()
        self.target_path = None

    def tearDown(self):
        if self.target_path is not None:
            os.remove(self.target_path)

    def test_write_easy(self):
        mviri_easy = self.create_easy_dataset()

        start = datetime(2011, 9, 12, 13, 24, 52)
        end = datetime(2011, 9, 12, 13, 27, 51)
        file_name = FCDRWriter.create_file_name_FCDR_easy("MVIRI", "Meteosat8", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.now()

        FCDRWriter.write(mviri_easy, self.target_path)

        end_time = datetime.now()
        elapsed_time = end_time - start_time
        print("MVIRI EASY write time: " + str(elapsed_time.seconds) + "." + str(round(elapsed_time.microseconds / 1000)))

        self.assertTrue(os.path.isfile(self.target_path))

        target_data = xr.open_dataset(self.target_path)
        try:
            self.assert_global_flags(target_data)
            self.assert_sensor_variables(target_data)

            variable = target_data["toa_bidirectional_reflectance_vis"]
            self.assertAlmostEqual(0.75, variable.data[25, 25], 5)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_independent_toa_bidirectional_reflectance"]
            self.assertAlmostEqual(1.0400092904, variable.data[26, 26], 4)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_structured_toa_bidirectional_reflectance"]
            self.assertAlmostEqual(1.3500070712000001, variable.data[27, 27], 4)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["sub_satellite_latitude_start"]
            self.assertAlmostEqual(25.25, variable.data, 8)

            variable = target_data["sub_satellite_longitude_start"]
            self.assertAlmostEqual(26.26, variable.data, 8)

            variable = target_data["sub_satellite_latitude_end"]
            self.assertAlmostEqual(27.27, variable.data, 8)

            variable = target_data["sub_satellite_longitude_end"]
            self.assertAlmostEqual(28.28, variable.data, 8)
        finally:
            target_data.close()

    def test_write_full(self):
        mviri_full = self.create_full_dataset()

        start = datetime(2010, 8, 11, 13, 24, 52)
        end = datetime(2010, 8, 11, 13, 27, 51)
        file_name = FCDRWriter.create_file_name_FCDR_full("MVIRI", "Meteosat7", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.now()

        FCDRWriter.write(mviri_full, self.target_path)

        end_time = datetime.now()
        elapsed_time = end_time - start_time
        print("MVIRI FULL write time: " + str(elapsed_time.seconds) + "." + str(round(elapsed_time.microseconds / 1000)))

        self.assertTrue(os.path.isfile(self.target_path))

        target_data = xr.open_dataset(self.target_path)
        try:
            self.assert_global_flags(target_data)
            self.assert_sensor_variables(target_data)

            variable = target_data["count_vis"]
            self.assertEqual(11, variable.data[11, 11])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_latitude"]
            self.assertAlmostEqual(0.21696, variable.data[12, 12], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_longitude"]
            self.assertAlmostEqual(0.633915, variable.data[13, 13], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_time"]
            self.assertAlmostEqual(0.402832012, variable.data[14], 8)
            self.assertEqual((2500,), variable.encoding["chunksizes"])

            variable = target_data["u_satellite_zenith_angle"]
            self.assertAlmostEqual(4.4999668098, variable.data[15, 15], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_satellite_azimuth_angle"]
            self.assertAlmostEqual(1.399993065, variable.data[16, 16], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_solar_zenith_angle"]
            self.assertAlmostEqual(3.4999826625, variable.data[17, 17], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["u_solar_azimuth_angle"]
            self.assertAlmostEqual(0.8000178354, variable.data[18, 18], 8)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["a0_vis"]
            self.assertAlmostEqual(7.7, variable.data, 8)

            variable = target_data["a1_vis"]
            self.assertAlmostEqual(8.8, variable.data, 8)

            variable = target_data["mean_count_space_vis"]
            self.assertAlmostEqual(9.9, variable.data, 8)

            variable = target_data["u_a0_vis"]
            self.assertAlmostEqual(10.1, variable.data, 8)

            variable = target_data["u_a1_vis"]
            self.assertAlmostEqual(11.11, variable.data, 8)

            variable = target_data["u_a2_vis"]
            self.assertAlmostEqual(12.12, variable.data, 8)

            variable = target_data["u_zero_vis"]
            self.assertAlmostEqual(13.13, variable.data, 8)

            variable = target_data["covariance_a_vis"]
            self.assertEqual((3,3), variable.shape)
            self.assertAlmostEqual(0.56, variable.data[1, 2], 8)

            variable = target_data["u_electronics_counts_vis"]
            self.assertAlmostEqual(13.13, variable.data, 8)

            variable = target_data["u_digitization_counts_vis"]
            self.assertAlmostEqual(14.14, variable.data, 8)

            variable = target_data["allan_deviation_counts_space_vis"]
            self.assertAlmostEqual(15.15, variable.data, 8)
        finally:
            target_data.close()

    def create_easy_dataset(self):
        mviri_easy = FCDRWriter.createTemplateEasy("MVIRI", 5000)
        self.add_global_attributes(mviri_easy)
        self.add_global_flags(mviri_easy)
        self.add_sensor_data(mviri_easy)

        for x in range(0, 5000):
            mviri_easy["toa_bidirectional_reflectance_vis"].data[:, x] = np.ones((5000), np.uint16) * x * 0.03
            mviri_easy["u_independent_toa_bidirectional_reflectance"].data[:, x] = np.ones((5000), np.uint16) * x * 0.04
            mviri_easy["u_structured_toa_bidirectional_reflectance"].data[:, x] = np.ones((5000), np.uint16) * x * 0.05

        mviri_easy["sub_satellite_latitude_start"].data = 25.25
        mviri_easy["sub_satellite_longitude_start"].data = 26.26
        mviri_easy["sub_satellite_latitude_end"].data = 27.27
        mviri_easy["sub_satellite_longitude_end"].data = 28.28

        return mviri_easy

    def create_full_dataset(self):
        mviri_full = FCDRWriter.createTemplateFull("MVIRI", 5000)
        self.add_global_attributes(mviri_full)
        self.add_global_flags(mviri_full)
        self.add_sensor_data(mviri_full)

        for x in range(0, 5000):
            mviri_full["count_vis"].data[:, x] = np.ones((5000), np.uint8) * x

        for x in range(0, 500):
            mviri_full["u_latitude"].data[:, x] = np.ones((500), np.float32) * x * 0.1
            mviri_full["u_longitude"].data[:, x] = np.ones((500), np.float32) * x * 0.2
            mviri_full["u_satellite_zenith_angle"].data[:, x] = np.ones((500), np.float32) * x * 0.3
            mviri_full["u_satellite_azimuth_angle"].data[:, x] = np.ones((500), np.float32) * x * 0.4
            mviri_full["u_solar_zenith_angle"].data[:, x] = np.ones((500), np.float32) * x * 0.5
            mviri_full["u_solar_azimuth_angle"].data[:, x] = np.ones((500), np.float32) * x * 0.6

        mviri_full["u_time"].data[:] = np.ones((2500), np.float32) * 0.4

        mviri_full["covariance_a_vis"].data = np.ones((3, 3), np.float32) * 0.56

        mviri_full["a0_vis"].data = 7.7
        mviri_full["a1_vis"].data = 8.8
        mviri_full["mean_count_space_vis"].data = 9.9
        mviri_full["u_a0_vis"].data = 10.1
        mviri_full["u_a1_vis"].data = 11.11
        mviri_full["u_a2_vis"].data = 12.12
        mviri_full["u_zero_vis"].data = 13.13
        mviri_full["u_electronics_counts_vis"].data = 13.13
        mviri_full["u_digitization_counts_vis"].data = 14.14
        mviri_full["allan_deviation_counts_space_vis"].data = 15.15

        return mviri_full

    def add_global_attributes(self, dataset):
        dataset.attrs["institution"] = "test"
        dataset.attrs["title"] = "sir"
        dataset.attrs["source"] = "invention"
        dataset.attrs["history"] = "new"
        dataset.attrs["references"] = "myself"
        dataset.attrs["comment"] = "should define a test version of this set"

    def add_global_flags(self, dataset):
        for x in range(0, PRODUCT_WIDTH):
            dataset["quality_pixel_bitmask"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int8) * x
            dataset["data_quality_bitmask"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int8) * x + 1

    def assert_global_flags(self, target_data):
        variable = target_data["quality_pixel_bitmask"]
        self.assertEqual(11, variable.data[10, 10])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

    def assert_sensor_variables(self, target_data):
        variable = target_data["time"]
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["solar_azimuth_angle"]
        self.assertAlmostEqual(0.010986328, variable.data[1, 1], 8)
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["solar_zenith_angle"]
        self.assertAlmostEqual(0.038452736, variable.data[2, 2], 8)
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["count_ir"]
        self.assertEqual(3, variable.data[3, 3])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["count_wv"]
        self.assertEqual(5, variable.data[4, 4])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["data_quality_bitmask"]
        self.assertEqual(6, variable.data[5, 5])
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["distance_sun_earth"]
        self.assertAlmostEqual(3.3, variable.data, 8)

        variable = target_data["solar_irradiance_vis"]
        self.assertAlmostEqual(4.4, variable.data, 8)

        variable = target_data["u_solar_irradiance_vis"]
        self.assertAlmostEqual(5.5, variable.data, 8)

        variable = target_data["SRF_weights"]
        self.assertAlmostEqual(0.159984, variable.data[0, 4], 8)
        self.assertEqual((3, 1011), variable.encoding["chunksizes"])

        variable = target_data["SRF_frequencies"]
        self.assertAlmostEqual(0.25, variable.data[1, 5], 8)
        self.assertEqual((3, 1011), variable.encoding["chunksizes"])

        variable = target_data["covariance_spectral_response_function_vis"]
        self.assertAlmostEqual(40.4, variable.data[5, 5], 5)
        self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        variable = target_data["u_spectral_response_function_ir"]
        self.assertAlmostEqual(0.05, variable.data[7], 8)
        self.assertEqual((1011,), variable.encoding["chunksizes"])

        variable = target_data["u_spectral_response_function_wv"]
        self.assertAlmostEqual(0.07, variable.data[9], 8)
        self.assertEqual((1011,), variable.encoding["chunksizes"])

        variable = target_data["a_ir"]
        self.assertAlmostEqual(8.8, variable.data, 8)

        variable = target_data["b_ir"]
        self.assertAlmostEqual(9.9, variable.data, 8)

        variable = target_data["u_a_ir"]
        self.assertAlmostEqual(10.1, variable.data, 8)

        variable = target_data["u_b_ir"]
        self.assertAlmostEqual(11.11, variable.data, 8)

        variable = target_data["a_wv"]
        self.assertAlmostEqual(12.12, variable.data, 8)

        variable = target_data["b_wv"]
        self.assertAlmostEqual(13.13, variable.data, 8)

        variable = target_data["u_a_wv"]
        self.assertAlmostEqual(14.14, variable.data, 8)

        variable = target_data["u_b_wv"]
        self.assertAlmostEqual(15.15, variable.data, 8)

        variable = target_data["bt_a_ir"]
        self.assertAlmostEqual(20.2, variable.data, 8)

        variable = target_data["bt_b_ir"]
        self.assertAlmostEqual(21.21, variable.data, 8)

        variable = target_data["bt_a_wv"]
        self.assertAlmostEqual(22.22, variable.data, 8)

        variable = target_data["bt_b_wv"]
        self.assertAlmostEqual(23.23, variable.data, 8)

        variable = target_data["years_since_launch"]
        self.assertAlmostEqual(24.24, variable.data, 8)

    def add_sensor_data(self, dataset):
        for x in range(0, 2500):
            dataset["time"].data[:, x] = np.ones((2500), np.uint32) * x
            dataset["count_ir"].data[:, x] = np.ones((2500), np.uint8) * x
            dataset["count_wv"].data[:, x] = np.ones((2500), np.uint8) * x + 1

        for x in range(0, 500):
            dataset["solar_azimuth_angle"].data[:, x] = np.ones((500), np.uint32) * x * 0.01
            dataset["solar_zenith_angle"].data[:, x] = np.ones((500), np.uint32) * x * 0.02

        for x in range(0, 1011):
            dataset["covariance_spectral_response_function_vis"].data[:] = np.ones((1011), np.float32) * x * 0.04

        for freq in range(0, 1011):
            dataset["SRF_weights"].data[:, freq] = np.ones((3), np.float32) * freq * 0.04
            dataset["SRF_frequencies"].data[:, freq] = np.ones((3), np.float32) * freq * 0.05

        dataset["u_spectral_response_function_ir"].data[:] = np.ones((1011), np.float32) * 0.05
        dataset["u_spectral_response_function_wv"].data[:] = np.ones((1011), np.float32) * 0.07

        dataset["distance_sun_earth"].data = 3.3
        dataset["solar_irradiance_vis"].data = 4.4
        dataset["u_solar_irradiance_vis"].data = 5.5
        dataset["a_ir"].data = 8.8
        dataset["b_ir"].data = 9.9
        dataset["u_a_ir"].data = 10.1
        dataset["u_b_ir"].data = 11.11
        dataset["a_wv"].data = 12.12
        dataset["b_wv"].data = 13.13
        dataset["u_a_wv"].data = 14.14
        dataset["u_b_wv"].data = 15.15
        dataset["bt_a_ir"].data = 20.2
        dataset["bt_b_ir"].data = 21.21
        dataset["bt_a_wv"].data = 22.22
        dataset["bt_b_wv"].data = 23.23
        dataset["years_since_launch"].data = 24.24
