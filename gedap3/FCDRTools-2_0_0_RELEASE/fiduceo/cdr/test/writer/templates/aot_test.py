import unittest

import numpy as np
import xarray as xr

from fiduceo.cdr.writer.templates.aot import AOT
from fiduceo.common.test.assertions import Assertions


class AOTTest(unittest.TestCase):

    CHUNKING = (1280, 409)

    def test_add_variables(self):
        ds = xr.Dataset()

        AOT.add_variables(ds, 409, 12876)

        Assertions.assert_geolocation_variables(self, ds, 409, 12876, chunking=self.CHUNKING)
        Assertions.assert_quality_flags(self, ds, 409, 12876, chunking=self.CHUNKING)

        time = ds.variables["time"]
        self.assertEqual((12876,), time.shape)
        self.assertEqual(-1, time.values[165])
        self.assertEqual(4294967295, time.attrs["_FillValue"])
        self.assertEqual("time", time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", time.attrs["long_name"])
        self.assertEqual("s", time.attrs["units"])

        aot = ds.variables["aot"]
        self.assertEqual((12876, 409), aot.shape)
        self.assertTrue(np.isnan(aot.values[176, 177]))
        self.assertTrue(np.isnan(aot.attrs["_FillValue"]))
        self.assertEqual("longitude latitude", aot.attrs["coordinates"])

        u_ind_aot = ds.variables["u_independent_aot"]
        self.assertEqual((12876, 409), u_ind_aot.shape)
        self.assertTrue(np.isnan(u_ind_aot.values[178, 179]))
        self.assertTrue(np.isnan(u_ind_aot.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of aot due to independent effects", u_ind_aot.attrs["description"])
        self.assertEqual("longitude latitude", u_ind_aot.attrs["coordinates"])

        u_str_aot = ds.variables["u_structured_aot"]
        self.assertEqual((12876, 409), u_str_aot.shape)
        self.assertTrue(np.isnan(u_str_aot.values[180, 181]))
        self.assertTrue(np.isnan(u_str_aot.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of aot due to structured effects", u_str_aot.attrs["description"])
        self.assertEqual("longitude latitude", u_str_aot.attrs["coordinates"])

        u_com_aot = ds.variables["u_common_aot"]
        self.assertEqual((12876, 409), u_com_aot.shape)
        self.assertTrue(np.isnan(u_com_aot.values[182, 183]))
        self.assertTrue(np.isnan(u_com_aot.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of aot due to common effects", u_com_aot.attrs["description"])
        self.assertEqual("longitude latitude", u_com_aot.attrs["coordinates"])
