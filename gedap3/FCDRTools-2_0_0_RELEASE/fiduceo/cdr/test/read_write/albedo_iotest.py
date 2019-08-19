import datetime
import os
import tempfile
import unittest

import numpy as np
import xarray as xr

from fiduceo.cdr.writer.cdr_writer import CDRWriter

PRODUCT_WIDTH = 5000
PRODUCT_HEIGHT = 5000

EXPECTED_CHUNKING = (500, 500)


class AlbedoIoTest(unittest.TestCase):
    temp_dir = None
    target_path = None

    def setUp(self):
        self.temp_dir = tempfile.gettempdir()
        self.target_path = None

    def tearDown(self):
        if self.target_path is not None:
            os.remove(self.target_path)

    def test_write_albedo(self):
        template = self.create_albedo_dataset()

        start = datetime.datetime(2016, 11, 22, 13, 24, 52)
        end = datetime.datetime(2016, 11, 22, 14, 25, 53)
        file_name = CDRWriter.create_file_name_CDR("ALBEDO", "MVIRI", "METEOSAT-7", start, end, "L2", "1.1")
        self.target_path = os.path.join(self.temp_dir, file_name)

        start_time = datetime.datetime.now()

        CDRWriter.write(template, self.target_path)

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        print("ALBEDO write time: " + str(elapsed_time.seconds) + "." + str(round(elapsed_time.microseconds / 1000)))

        self.assertTrue(os.path.isfile(self.target_path))

        target_data = xr.open_dataset(self.target_path)
        try:
            variable = target_data["quality_pixel_bitmask"]
            self.assertEqual(8, variable.data[1, 1])
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

            variable = target_data["time"]
            self.assertEqual(8, variable.data[2])

            variable = target_data["surface_albedo"]
            self.assertAlmostEqual(2.7, variable.data[3, 3], 7)
            self.assertEqual(EXPECTED_CHUNKING, variable.encoding["chunksizes"])

        finally:
            target_data.close()

    def create_albedo_dataset(self):
        dataset = CDRWriter.createTemplate("ALBEDO", PRODUCT_WIDTH, PRODUCT_HEIGHT)

        dataset.attrs["institution"] = "Brockmann Consult GmbH"
        dataset.attrs["title"] = "integration test dataset"
        dataset.attrs["source"] = "made-up data"
        dataset.attrs["history"] = "initial version"
        dataset.attrs["references"] = "no"
        dataset.attrs["comment"] = "just testing"

        dataset.attrs["auxiliary_data"] = "nope"
        dataset.attrs["configuration"] = "test"
        dataset.attrs["time_coverage_start"] = "now"
        dataset.attrs["time_coverage_end"] = "later"
        dataset.attrs["time_coverage_duration"] = "P8S"
        dataset.attrs["time_coverage_resolution"] = "P1S"

        for x in range(0, PRODUCT_WIDTH):
            dataset["quality_pixel_bitmask"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int8) * x + 7
            dataset["surface_albedo"].data[:, x] = np.ones((PRODUCT_HEIGHT), np.int8) * x * 0.9

        dataset["time"].data[:] = np.ones((PRODUCT_HEIGHT), np.int32) * 8

        return dataset
