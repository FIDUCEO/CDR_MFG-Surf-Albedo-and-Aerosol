import unittest

import numpy as np
import xarray as xr

from fiduceo.cdr.writer.templates.sst_ensemble import SST_ENSEMBLE
from fiduceo.common.test.assertions import Assertions


class SSTEnsembleTest(unittest.TestCase):

    CHUNKING = (1280, 409)

    def test_add_variables(self):
        ds = xr.Dataset()

        SST_ENSEMBLE.add_variables(ds, 409, 12623, 11)

        Assertions.assert_geolocation_variables(self, ds, 409, 12623, chunking=self.CHUNKING)
        Assertions.assert_quality_flags(self, ds, 409, 12623, chunking=self.CHUNKING)

        time = ds.variables["time"]
        self.assertEqual((12623,), time.shape)
        self.assertEqual(-1, time.values[166])
        self.assertEqual(4294967295, time.attrs["_FillValue"])
        self.assertEqual("time", time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", time.attrs["long_name"])
        self.assertEqual("s", time.attrs["units"])

        sst = ds.variables["sst"]
        self.assertEqual((11, 12623, 409), sst.shape)
        self.assertTrue(np.isnan(sst.values[6, 198, 199]))
        self.assertTrue(np.isnan(sst.attrs["_FillValue"]))
        self.assertEqual("sea_surface_temperature", sst.attrs["standard_name"])
        self.assertEqual("K", sst.attrs["units"])
        self.assertEqual("longitude latitude", sst.attrs["coordinates"])