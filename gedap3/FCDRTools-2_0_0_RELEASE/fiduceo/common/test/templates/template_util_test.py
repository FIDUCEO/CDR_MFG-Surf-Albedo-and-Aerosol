import unittest

import numpy as np
import xarray as xr
from xarray import Variable

from fiduceo.common.writer.templates.templateutil import TemplateUtil
from fiduceo.common.writer.default_data import DefaultData


class TemplateUtilTest(unittest.TestCase):

    def test_add_geolocation_variables(self):
        ds = xr.Dataset()
        TemplateUtil.add_geolocation_variables(ds, 8, 10)

        latitude = ds.variables["latitude"]
        self.assertEqual((10, 8), latitude.shape)
        self.assertTrue(np.isnan(latitude.data[4, 4]))
        self.assertEqual("latitude", latitude.attrs["standard_name"])
        self.assertEqual("degrees_north", latitude.attrs["units"])
        self.assertEqual(np.int16, latitude.encoding['dtype'])
        self.assertEqual(-32768, latitude.encoding['_FillValue'])
        self.assertEqual(0.0027466658, latitude.encoding['scale_factor'])
        self.assertEqual(0.0, latitude.encoding['add_offset'])

        longitude = ds.variables["longitude"]
        self.assertEqual((10, 8), longitude.shape)
        self.assertTrue(np.isnan(longitude.data[5, 6]))
        self.assertEqual("longitude", longitude.attrs["standard_name"])
        self.assertEqual("degrees_east", longitude.attrs["units"])
        self.assertEqual(np.int16, longitude.encoding['dtype'])
        self.assertEqual(-32768, longitude.encoding['_FillValue'])
        self.assertEqual(0.0054933317, longitude.encoding['scale_factor'])
        self.assertEqual(0.0, longitude.encoding['add_offset'])

    def test_add_gridded_geolocation_variables(self):
        ds = xr.Dataset()
        TemplateUtil.add_gridded_geolocation_variables(ds, 9, 11)

        lat = ds.variables["lat"]
        self.assertEqual((11,), lat.shape)
        self.assertTrue(np.isnan(lat.data[5]))
        self.assertTrue(np.isnan(lat.attrs['_FillValue']))
        self.assertEqual("latitude", lat.attrs["standard_name"])
        self.assertEqual("latitude", lat.attrs["long_name"])
        self.assertEqual("degrees_north", lat.attrs["units"])
        self.assertEqual("lat_bnds", lat.attrs["bounds"])

        lat_bnds = ds.variables["lat_bnds"]
        self.assertEqual((11, 2), lat_bnds.shape)
        self.assertTrue(np.isnan(lat_bnds.data[6, 0]))
        self.assertTrue(np.isnan(lat_bnds.attrs['_FillValue']))
        self.assertEqual("latitude cell boundaries", lat_bnds.attrs["long_name"])
        self.assertEqual("degrees_north", lat_bnds.attrs["units"])

        lon = ds.variables["lon"]
        self.assertEqual((9,), lon.shape)
        self.assertTrue(np.isnan(lon.data[6]))
        self.assertTrue(np.isnan(lon.attrs['_FillValue']))
        self.assertEqual("longitude", lon.attrs["standard_name"])
        self.assertEqual("longitude", lon.attrs["long_name"])
        self.assertEqual("degrees_east", lon.attrs["units"])
        self.assertEqual("lon_bnds", lon.attrs["bounds"])

        lon_bnds = ds.variables["lon_bnds"]
        self.assertEqual((9, 2), lon_bnds.shape)
        self.assertTrue(np.isnan(lon_bnds.data[7, 1]))
        self.assertTrue(np.isnan(lon_bnds.attrs['_FillValue']))
        self.assertEqual("longitude cell boundaries", lon_bnds.attrs["long_name"])
        self.assertEqual("degrees_east", lon_bnds.attrs["units"])

    def test_add_quality_flags(self):
        ds = xr.Dataset()
        TemplateUtil.add_quality_flags(ds, 9, 11)

        quality = ds.variables["quality_pixel_bitmask"]
        self.assertEqual((11, 9), quality.shape)
        self.assertEqual(0, quality.data[5, 5])
        self.assertEqual(np.uint8, quality.dtype)
        self.assertEqual("status_flag", quality.attrs["standard_name"])
        self.assertEqual("1, 2, 4, 8, 16, 32, 64, 128", quality.attrs["flag_masks"])
        self.assertEqual("invalid use_with_caution invalid_input invalid_geoloc invalid_time sensor_error padded_data incomplete_channel_data", quality.attrs["flag_meanings"])

    def test_add_coordinates(self):
        ds = xr.Dataset()
        TemplateUtil.add_geolocation_variables(ds, 13, 108)

        TemplateUtil.add_coordinates(ds)

        x = ds["x"]
        self.assertEqual((13,), x.shape)
        self.assertEqual(np.uint16, x.dtype)

        y = ds["y"]
        self.assertEqual((108,), y.shape)
        self.assertEqual(np.uint16, y.dtype)

    def test_add_coordinates_with_channel(self):
        ds = xr.Dataset()
        TemplateUtil.add_geolocation_variables(ds, 8, 11)
        default_array = DefaultData.create_default_array_3d(8, 11, 4, np.float32, np.NaN)
        ds["three_d"] = Variable(["channel", "y", "x"], default_array)

        TemplateUtil.add_coordinates(ds)

        x = ds["x"]
        self.assertEqual((8,), x.shape)
        self.assertEqual(np.uint16, x.dtype)

        y = ds["y"]
        self.assertEqual((11,), y.shape)
        self.assertEqual(np.uint16, y.dtype)

        channel = ds["channel"]
        self.assertEqual((4,), channel.shape)
        self.assertEqual(np.uint16, channel.dtype)

    def test_add_geolocation_attribute(self):
        default_array = DefaultData.create_default_array_3d(8, 6, 4, np.float32, np.NaN)
        variable = Variable(["channel", "y", "x"], default_array)

        TemplateUtil.add_geolocation_attribute(variable)
        self.assertEqual("longitude latitude", variable.attrs["coordinates"])
