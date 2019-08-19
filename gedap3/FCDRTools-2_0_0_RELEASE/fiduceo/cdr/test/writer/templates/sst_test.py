import unittest

import numpy as np
import xarray as xr

from fiduceo.cdr.writer.templates.sst import SST
from fiduceo.common.test.assertions import Assertions


class SSTTest(unittest.TestCase):
    CHUNKING = (1280, 409)

    def test_add_variables(self):
        ds = xr.Dataset()

        SST.add_variables(ds, 409, 12877)

        Assertions.assert_geolocation_variables(self, ds, 409, 12877, chunking=self.CHUNKING)
        Assertions.assert_quality_flags(self, ds, 409, 12877, chunking=self.CHUNKING)

        time = ds.variables["time"]
        self.assertEqual((12877,), time.shape)
        self.assertEqual(-1, time.values[166])
        self.assertEqual(4294967295, time.attrs["_FillValue"])
        self.assertEqual("time", time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", time.attrs["long_name"])
        self.assertEqual("s", time.attrs["units"])

        sst = ds.variables["sst"]
        self.assertEqual((12877, 409), sst.shape)
        self.assertTrue(np.isnan(sst.values[196, 197]))
        self.assertTrue(np.isnan(sst.attrs["_FillValue"]))
        self.assertEqual("sea_surface_temperature", sst.attrs["standard_name"])
        self.assertEqual("K", sst.attrs["units"])
        self.assertEqual("longitude latitude", sst.attrs["coordinates"])

        u_ind_sst = ds.variables["u_independent_sst"]
        self.assertEqual((12877, 409), u_ind_sst.shape)
        self.assertTrue(np.isnan(u_ind_sst.values[198, 199]))
        self.assertTrue(np.isnan(u_ind_sst.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of sst due to independent effects", u_ind_sst.attrs["description"])
        self.assertEqual("longitude latitude", u_ind_sst.attrs["coordinates"])

        u_str_sst = ds.variables["u_structured_sst"]
        self.assertEqual((12877, 409), u_str_sst.shape)
        self.assertTrue(np.isnan(u_str_sst.values[200, 201]))
        self.assertTrue(np.isnan(u_str_sst.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of sst due to structured effects", u_str_sst.attrs["description"])
        self.assertEqual("longitude latitude", u_str_sst.attrs["coordinates"])

        u_com_sst = ds.variables["u_common_sst"]
        self.assertEqual((12877, 409), u_com_sst.shape)
        self.assertTrue(np.isnan(u_com_sst.values[170, 171]))
        self.assertTrue(np.isnan(u_com_sst.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of sst due to common effects", u_com_sst.attrs["description"])
        self.assertEqual("longitude latitude", u_com_sst.attrs["coordinates"])
