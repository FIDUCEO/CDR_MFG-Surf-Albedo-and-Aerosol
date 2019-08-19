import unittest

import datetime

from fiduceo.cdr.writer.cdr_writer import CDRWriter
from fiduceo.common.test.assertions import Assertions


class CDRWriterTest(unittest.TestCase):

    def testTemplateEasy_ALBEDO(self):
        ds = CDRWriter.createTemplate('ALBEDO', 5000, 5000)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)
        Assertions.assert_cdr_global_attributes(self, ds.attrs)

    def testTemplateEasy_AOT(self):
        ds = CDRWriter.createTemplate('AOT', 409, 12941)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)
        Assertions.assert_cdr_global_attributes(self, ds.attrs)

    def testTemplateEasy_SST(self):
        ds = CDRWriter.createTemplate('SST', 409, 11732)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)
        Assertions.assert_cdr_global_attributes(self, ds.attrs)

    def testTemplateEasy_SST_ENSEMBLE(self):
        ds = CDRWriter.createTemplate('SST_ENSEMBLE', 409, 11732, 20)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)
        Assertions.assert_cdr_global_attributes(self, ds.attrs)

    def testTemplateEasy_UTH(self):
        ds = CDRWriter.createTemplate('UTH', 360, 120)
        self.assertIsNotNone(ds)

        Assertions.assert_global_attributes(self, ds.attrs)
        Assertions.assert_cdr_global_attributes(self, ds.attrs)
        Assertions.assert_gridded_global_attributes(self, ds.attrs)

    def test_create_file_name_CDR(self):
        start = datetime.datetime(2015, 8, 23, 14, 24, 52)
        end = datetime.datetime(2015, 8, 23, 15, 25, 53)
        self.assertEqual("FIDUCEO_CDR_AOT_MVIRI_MET7-0.00_20150823142452_20150823152553_L2_v02.3_fv2.0.0.nc",
                         CDRWriter.create_file_name_CDR("AOT", "MVIRI", "MET7-0.00", start, end, "L2", "02.3"))

        start = datetime.datetime(2018, 1, 9, 23, 24, 52)
        end = datetime.datetime(2018, 1, 10, 1, 25, 53)
        self.assertEqual("FIDUCEO_CDR_SST_AVHRR_NOAA18_20180109232452_20180110012553_ENSEMBLE_v03.3_fv2.0.0.nc",
                         CDRWriter.create_file_name_CDR("SST", "AVHRR", "NOAA18", start, end, "ENSEMBLE", "03.3"))