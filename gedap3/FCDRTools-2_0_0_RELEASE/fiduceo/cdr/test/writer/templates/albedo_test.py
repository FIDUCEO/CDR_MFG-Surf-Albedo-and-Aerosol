import unittest

import numpy as np
import xarray as xr

from fiduceo.cdr.writer.templates.albedo import Albedo
from fiduceo.common.test.assertions import Assertions


class AlbedoTest(unittest.TestCase):
    CHUNKING = (500, 500)

    def test_add_variables(self):
        ds = xr.Dataset()
        Albedo.add_variables(ds, 5000, 5000)

        Assertions.assert_quality_flags(self, ds, 5000, 5000, chunking=self.CHUNKING)

        time = ds.variables["time"]
        self.assertEqual((5000,), time.shape)
        self.assertEqual(-1, time.values[165])
        self.assertEqual(-1, time.attrs["_FillValue"])
        self.assertEqual("time", time.attrs["standard_name"])
        self.assertEqual("Acquisition time in seconds since 1970-01-01 00:00:00", time.attrs["long_name"])
        self.assertEqual("s", time.attrs["units"])

        albedo = ds.variables["surface_albedo"]
        self.assertEqual((5000, 5000), albedo.shape)
        self.assertTrue(np.isnan(albedo.values[166, 167]))
        self.assertTrue(np.isnan(albedo.attrs["_FillValue"]))
        self.assertEqual("surface_albedo", albedo.attrs["standard_name"])
        self.assertEqual("longitude latitude", albedo.attrs["coordinates"])

        u_ind_albedo = ds.variables["u_independent_surface_albedo"]
        self.assertEqual((5000, 5000), u_ind_albedo.shape)
        self.assertTrue(np.isnan(u_ind_albedo.values[168, 169]))
        self.assertTrue(np.isnan(u_ind_albedo.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of surface_albedo due to independent effects", u_ind_albedo.attrs["description"])
        self.assertEqual("longitude latitude", u_ind_albedo.attrs["coordinates"])

        u_str_albedo = ds.variables["u_structured_surface_albedo"]
        self.assertEqual((5000, 5000), u_str_albedo.shape)
        self.assertTrue(np.isnan(u_str_albedo.values[170, 171]))
        self.assertTrue(np.isnan(u_str_albedo.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of surface_albedo due to structured effects", u_str_albedo.attrs["description"])
        self.assertEqual("longitude latitude", u_str_albedo.attrs["coordinates"])

        u_com_albedo = ds.variables["u_common_surface_albedo"]
        self.assertEqual((5000, 5000), u_com_albedo.shape)
        self.assertTrue(np.isnan(u_com_albedo.values[170, 171]))
        self.assertTrue(np.isnan(u_com_albedo.attrs["_FillValue"]))
        self.assertEqual("Uncertainty of surface_albedo due to common effects", u_com_albedo.attrs["description"])
        self.assertEqual("longitude latitude", u_com_albedo.attrs["coordinates"])
