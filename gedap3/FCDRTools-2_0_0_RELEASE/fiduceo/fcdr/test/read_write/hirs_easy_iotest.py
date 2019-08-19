import os
import tempfile
import unittest
from datetime import datetime

import numpy as np
import xarray as xr

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter

EXPECTED_CHUNKING_2D = (512, 56)
EXPECTED_CHUNKING_3D = (10, 512, 56)

SRF_SIZE = 103
NUM_CHANNELS = 19


class HirsEASYIoTest(unittest.TestCase):
    temp_dir = None
    target_path = None

    def setUp(self):
        self.temp_dir = tempfile.gettempdir()
        self.target_path = None

    def tearDown(self):
        if self.target_path is not None:
            os.remove(self.target_path)

    def test_write_HIRS2(self):
        hirs_easy = self.create_easy_dataset("HIRS2")

        start = datetime(2015, 10, 21, 13, 24, 52)
        end = datetime(2015, 10, 21, 14, 25, 53)
        file_name = FCDRWriter.create_file_name_FCDR_easy("HIRS2", "NOAA12", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.now()

        FCDRWriter.write(hirs_easy, self.target_path)

        end_time = datetime.now()
        elapsed_time = end_time - start_time
        print("HIRS/2 EASY write time: " + "%d.%03d" % (elapsed_time.seconds, int(round(elapsed_time.microseconds / 1000))))

        self.assertTrue(os.path.isfile(self.target_path))
        target_data = xr.open_dataset(self.target_path)
        try:
            variable = target_data["bt"]
            self.assertAlmostEqual(0.0, variable.data[0, 0, 0], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["latitude"]
            self.assertAlmostEqual(0.0192266606, variable.data[1, 1], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["longitude"]
            self.assertAlmostEqual(0.0604266487, variable.data[2, 2], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["quality_pixel_bitmask"]
            self.assertEqual(3, variable.data[3, 3])
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["data_quality_bitmask"]
            self.assertEqual(2, variable.data[9, 2])
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["quality_scanline_bitmask"]
            self.assertEqual(1, variable.data[4])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["satellite_zenith_angle"]
            self.assertEqual(2, variable.data[5])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["scanline"]
            self.assertEqual(3, variable.data[6])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["solar_azimuth_angle"]
            self.assertAlmostEqual(0.35, variable.data[7, 7], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["time"]
            self.assertEqual(4, variable.data[8])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["u_independent"]
            self.assertAlmostEqual(0.54, variable.data[9, 9, 9], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["u_structured"]
            self.assertAlmostEqual(0.7, variable.data[10, 10, 10], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["u_common"]
            self.assertAlmostEqual(0.88, variable.data[11, 11, 11], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            self.assert_srf(target_data)
        finally:
            target_data.close()

    def test_write_HIRS3(self):
        hirs_easy = self.create_easy_dataset("HIRS3")

        start = datetime(2014, 9, 20, 13, 24, 52)
        end = datetime(2014, 9, 20, 14, 25, 53)
        file_name = FCDRWriter.create_file_name_FCDR_easy("HIRS3", "NOAA15", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.now()

        FCDRWriter.write(hirs_easy, self.target_path)

        end_time = datetime.now()
        elapsed_time = end_time - start_time
        print("HIRS/3 EASY write time: " + "%d.%03d" % (elapsed_time.seconds, int(round(elapsed_time.microseconds / 1000))))

        self.assertTrue(os.path.isfile(self.target_path))
        target_data = xr.open_dataset(self.target_path)
        try:
            variable = target_data["bt"]
            self.assertAlmostEqual(0.01, variable.data[1, 1, 1], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["latitude"]
            self.assertAlmostEqual(0.041199987, variable.data[2, 2], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["longitude"]
            self.assertAlmostEqual(0.0878933072, variable.data[3, 3], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["quality_pixel_bitmask"]
            self.assertEqual(6, variable.data[4, 4])
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["data_quality_bitmask"]
            self.assertEqual(3, variable.data[10, 3])
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["quality_scanline_bitmask"]
            self.assertEqual(1, variable.data[5])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["quality_channel_bitmask"]
            self.assertEqual(18, variable.data[6, 6])
            self.assertEqual((944, 19), variable.encoding["chunksizes"])

            variable = target_data["satellite_zenith_angle"]
            self.assertAlmostEqual(0.24, variable.data[6, 6], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["satellite_azimuth_angle"]
            self.assertAlmostEqual(0.35, variable.data[7, 7], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["scanline"]
            self.assertEqual(3, variable.data[7])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["solar_azimuth_angle"]
            self.assertAlmostEqual(0.4, variable.data[8, 8], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["solar_zenith_angle"]
            self.assertAlmostEqual(0.54, variable.data[9, 9], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["time"]
            self.assertEqual(4, variable.data[9])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["u_independent"]
            self.assertAlmostEqual(0.6, variable.data[10, 10, 10], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["u_structured"]
            self.assertAlmostEqual(0.77, variable.data[11, 11, 11], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["u_common"]
            self.assertAlmostEqual(0.96, variable.data[12, 12, 12], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            self.assert_srf(target_data)
        finally:
            target_data.close()

    def test_write_HIRS4(self):
        hirs_easy = self.create_easy_dataset("HIRS4")

        start = datetime(2013, 8, 19, 13, 24, 52)
        end = datetime(2013, 8, 19, 14, 25, 53)
        file_name = FCDRWriter.create_file_name_FCDR_easy("HIRS4", "NOAA18", start, end, "1.0")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.now()

        FCDRWriter.write(hirs_easy, self.target_path)

        end_time = datetime.now()
        elapsed_time = end_time - start_time
        print("HIRS/4 EASY write time: " + "%d.%03d" % (elapsed_time.seconds, int(round(elapsed_time.microseconds / 1000))))

        self.assertTrue(os.path.isfile(self.target_path))
        target_data = xr.open_dataset(self.target_path)
        try:
            variable = target_data["bt"]
            self.assertAlmostEqual(0.02, variable.data[2, 2, 2], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["latitude"]
            self.assertAlmostEqual(0.0604266476, variable.data[3, 3], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["longitude"]
            self.assertAlmostEqual(0.1208532974, variable.data[4, 4], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["quality_pixel_bitmask"]
            self.assertEqual(7, variable.data[5, 5])
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["data_quality_bitmask"]
            self.assertEqual(4, variable.data[11, 4])
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["quality_scanline_bitmask"]
            self.assertEqual(1, variable.data[6])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["quality_channel_bitmask"]
            self.assertEqual(21, variable.data[7, 7])
            self.assertEqual((944, 19), variable.encoding["chunksizes"])

            variable = target_data["satellite_zenith_angle"]
            self.assertAlmostEqual(0.28, variable.data[7, 7], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["satellite_azimuth_angle"]
            self.assertAlmostEqual(0.4, variable.data[8, 8], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["scanline"]
            self.assertEqual(3, variable.data[8])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["solar_azimuth_angle"]
            self.assertAlmostEqual(0.45, variable.data[9, 9], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["solar_zenith_angle"]
            self.assertAlmostEqual(0.6, variable.data[10, 10], 8)
            self.assertEqual(EXPECTED_CHUNKING_2D, variable.encoding["chunksizes"])

            variable = target_data["time"]
            self.assertEqual(4, variable.data[10])
            self.assertEqual((944,), variable.encoding["chunksizes"])

            variable = target_data["u_independent"]
            self.assertAlmostEqual(0.66, variable.data[11, 11, 11], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["u_structured"]
            self.assertAlmostEqual(0.84, variable.data[12, 12, 12], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            variable = target_data["u_common"]
            self.assertAlmostEqual(1.04, variable.data[13, 13, 13], 8)
            self.assertEqual(EXPECTED_CHUNKING_3D, variable.encoding["chunksizes"])

            self.assert_srf(target_data)
        finally:
            target_data.close()

    def create_easy_dataset(self, type):
        hirs_easy = FCDRWriter.createTemplateEasy(type, 944, SRF_SIZE)
        hirs_easy.attrs["institution"] = "test"
        hirs_easy.attrs["title"] = "sir"
        hirs_easy.attrs["source"] = "invention"
        hirs_easy.attrs["history"] = "new"
        hirs_easy.attrs["references"] = "myself"
        hirs_easy.attrs["comment"] = "should define a test version of this set"

        for x in range(0, 56):
            hirs_easy["bt"].data[:, :, x] = np.ones((944), np.int16) * x * 0.01
            hirs_easy["latitude"].data[:, x] = np.ones((944), np.int16) * x * 0.02
            hirs_easy["longitude"].data[:, x] = np.ones((944), np.int16) * x * 0.03
            hirs_easy["data_quality_bitmask"].data[:, x] = np.ones((944), np.int8) * x
            hirs_easy["quality_pixel_bitmask"].data[:, x] = np.ones((944), np.int8) * x
            if type != "HIRS2":
                hirs_easy["satellite_zenith_angle"].data[:, x] = np.ones((944), np.int16) * x * 0.04
                hirs_easy["satellite_azimuth_angle"].data[:, x] = np.ones((944), np.int16) * x * 0.05
                hirs_easy["solar_zenith_angle"].data[:, x] = np.ones((944), np.int16) * x * 0.06

            hirs_easy["solar_azimuth_angle"].data[:, x] = np.ones((944), np.int16) * x * 0.05
            hirs_easy["u_independent"].data[:, :, x] = np.ones((944), np.int16) * x * 0.06
            hirs_easy["u_structured"].data[:, :, x] = np.ones((944), np.int16) * x * 0.07
            hirs_easy["u_common"].data[:, :, x] = np.ones((944), np.int16) * x * 0.08

        if type != "HIRS2":
            for x in range(0, 19):
                hirs_easy["quality_channel_bitmask"].data[:, x] = np.ones((944), np.int32) * x * 3

        hirs_easy["quality_scanline_bitmask"].data[:] = np.ones((944), np.int8)
        hirs_easy["scanline"].data[:] = np.ones((944), np.int8) * 3
        hirs_easy["time"].data[:] = np.ones((944), np.int32) * 4

        for x in range(0, SRF_SIZE):
            hirs_easy["SRF_weights"].data[:, x] = np.ones((NUM_CHANNELS), np.float32) * x * 0.04
            hirs_easy["SRF_wavelengths"].data[:, x] = np.ones((NUM_CHANNELS), np.float32) * x * 0.05

        if type == "HIRS2":
            hirs_easy["satellite_zenith_angle"].data[:] = np.ones((944), np.int8) * 2

        return hirs_easy

    def assert_srf(self, dataset):
        variable = dataset["SRF_weights"]
        self.assertAlmostEqual(0.439989, variable.data[1, 11], 8)

        variable = dataset["SRF_wavelengths"]
        self.assertAlmostEqual(0.6, variable.data[2, 12], 8)
