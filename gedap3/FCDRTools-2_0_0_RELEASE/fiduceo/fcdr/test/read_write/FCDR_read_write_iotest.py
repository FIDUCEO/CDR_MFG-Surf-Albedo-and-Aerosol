import os
import tempfile
import unittest

import numpy as np
import xarray as xr
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter


class FCDRReadWriteTests(unittest.TestCase):

    def setUp(self):
        self.dataset = xr.Dataset()
        self.dataset.attrs["template_key"] = "HIRS2"

        default_array = DefaultData.create_default_array(5, 5, np.uint16, fill_value=0)
        variable = Variable(["y", "x"], default_array)
        self.dataset["data_quality_bitmask"] = variable
        self.dataset["quality_pixel_bitmask"] = variable

        default_array = DefaultData.create_default_vector(5, np.int32, fill_value=0)
        variable = Variable(["y"], default_array)
        self.dataset["quality_scanline_bitmask"] = variable

        default_array = DefaultData.create_default_array(19, 5, np.uint8, fill_value=0)
        variable = Variable(["y", "channel"], default_array)
        self.dataset["quality_channel_bitmask"] = variable

        tempDir = tempfile.gettempdir()
        self.testDir = os.path.join(tempDir, 'fcdr')
        os.mkdir(self.testDir)

    def tearDown(self):
        if os.path.isdir(self.testDir):
            for i in os.listdir(self.testDir):
                os.remove(os.path.join(self.testDir, i))
            os.rmdir(self.testDir)

    def test_write_overwrite_true(self):
        testFile = os.path.join(self.testDir, 'delete_me.nc')

        FCDRWriter.write(self.dataset, testFile, overwrite=True)
        self.assertTrue(os.path.isfile(testFile))

        FCDRWriter.write(self.dataset, testFile, overwrite=True)
        self.assertTrue(os.path.isfile(testFile))

    def test_write_overwrite_false(self):
        testFile = os.path.join(self.testDir, 'delete_me.nc')

        FCDRWriter.write(self.dataset, testFile)
        self.assertTrue(os.path.isfile(testFile))

        try:
            FCDRWriter.write(self.dataset, testFile, overwrite=False)
            self.fail("IOError expected")
        except IOError:
            pass
