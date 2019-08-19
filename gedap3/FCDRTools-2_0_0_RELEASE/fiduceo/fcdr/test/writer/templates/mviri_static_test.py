import unittest

import numpy as np
import xarray as xr

from fiduceo.fcdr.writer.templates.mviri_static import MVIRI_STATIC


class MVIRI_STATICTest(unittest.TestCase):
    def test_add_original_variables(self):
        ds = xr.Dataset()
        MVIRI_STATIC.add_original_variables(ds, 5)

        lat_vis = ds.variables["latitude_vis"]
        self.assertEqual((5000, 5000), lat_vis.shape)
        self.assertTrue(np.isnan(lat_vis.data[4, 5]))
        self.assertEqual("latitude", lat_vis.attrs["standard_name"])
        self.assertEqual("degrees_north", lat_vis.attrs["units"])
        self.assertEqual(np.int16, lat_vis.encoding['dtype'])
        self.assertEqual(-32768, lat_vis.encoding['_FillValue'])
        self.assertEqual(0.0027466658, lat_vis.encoding['scale_factor'])
        self.assertEqual(0.0, lat_vis.encoding['add_offset'])

        lon_vis = ds.variables["longitude_vis"]
        self.assertEqual((5000, 5000), lon_vis.shape)
        self.assertTrue(np.isnan(lon_vis.data[5, 5]))
        self.assertEqual("longitude", lon_vis.attrs["standard_name"])
        self.assertEqual("degrees_east", lon_vis.attrs["units"])
        self.assertEqual(np.int16, lon_vis.encoding['dtype'])
        self.assertEqual(-32768, lon_vis.encoding['_FillValue'])
        self.assertEqual(0.0054933317, lon_vis.encoding['scale_factor'])
        self.assertEqual(0.0, lon_vis.encoding['add_offset'])

        lat_ir_wv = ds.variables["latitude_ir_wv"]
        self.assertEqual((2500, 2500), lat_ir_wv.shape)
        self.assertTrue(np.isnan(lat_ir_wv.data[5, 6]))
        self.assertEqual("latitude", lat_ir_wv.attrs["standard_name"])
        self.assertEqual("degrees_north", lat_ir_wv.attrs["units"])
        self.assertEqual(np.int16, lat_ir_wv.encoding['dtype'])
        self.assertEqual(-32768, lat_ir_wv.encoding['_FillValue'])
        self.assertEqual(0.0027466658, lat_ir_wv.encoding['scale_factor'])
        self.assertEqual(0.0, lat_ir_wv.encoding['add_offset'])

        lon_ir_wv = ds.variables["longitude_ir_wv"]
        self.assertEqual((2500, 2500), lon_ir_wv.shape)
        self.assertTrue(np.isnan(lon_ir_wv.data[6, 6]))
        self.assertEqual("longitude", lon_ir_wv.attrs["standard_name"])
        self.assertEqual("degrees_east", lon_ir_wv.attrs["units"])
        self.assertEqual(np.int16, lon_ir_wv.encoding['dtype'])
        self.assertEqual(-32768, lon_ir_wv.encoding['_FillValue'])
        self.assertEqual(0.0054933317, lon_ir_wv.encoding['scale_factor'])
        self.assertEqual(0.0, lon_ir_wv.encoding['add_offset'])

    def test_get_swath_width(self):
        self.assertEqual(5000, MVIRI_STATIC.get_swath_width())

    def test_add_template_key(self):
        ds = xr.Dataset()

        MVIRI_STATIC.add_template_key(ds)

        self.assertEqual("MVIRI_STATIC", ds.attrs["template_key"])