import os
import tempfile
import unittest

import xarray as xr

from fiduceo.cdr.writer.cdr_writer import CDRWriter


class CDRReadWriteTests(unittest.TestCase):

    def setUp(self):
        self.dataset = xr.Dataset()

        tempDir = tempfile.gettempdir()
        self.testDir = os.path.join(tempDir, 'cdr')
        os.mkdir(self.testDir)

    def tearDown(self):
        if os.path.isdir(self.testDir):
            for i in os.listdir(self.testDir):
                os.remove(os.path.join(self.testDir, i))
            os.rmdir(self.testDir)

    def test_write_overwrite_true(self):
        testFile = os.path.join(self.testDir, 'delete_me.nc')

        CDRWriter.write(self.dataset, testFile, overwrite=True)
        self.assertTrue(os.path.isfile(testFile))

        CDRWriter.write(self.dataset, testFile, overwrite=True)
        self.assertTrue(os.path.isfile(testFile))

    def test_write_overwrite_false(self):
        testFile = os.path.join(self.testDir, 'delete_me.nc')

        CDRWriter.write(self.dataset, testFile)
        self.assertTrue(os.path.isfile(testFile))

        try:
            CDRWriter.write(self.dataset, testFile, overwrite=False)
            self.fail("IOError expected")
        except IOError:
            pass
