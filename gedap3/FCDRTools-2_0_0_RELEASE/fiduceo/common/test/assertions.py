import numpy as np


class Assertions:

    @staticmethod
    def assert_geolocation_variables(test_case, dataset, width, height, chunking=None):
        latitude = dataset.variables["latitude"]
        test_case.assertEqual((height, width), latitude.shape)
        test_case.assertTrue(np.isnan(latitude.data[1, 1]))
        test_case.assertEqual("latitude", latitude.attrs["standard_name"])
        test_case.assertEqual("degrees_north", latitude.attrs["units"])
        test_case.assertEqual(np.int16, latitude.encoding['dtype'])
        test_case.assertEqual(-32768, latitude.encoding['_FillValue'])
        test_case.assertEqual(0.0027466658, latitude.encoding['scale_factor'])
        test_case.assertEqual(0.0, latitude.encoding['add_offset'])
        if chunking is not None:
            test_case.assertEqual(chunking, latitude.encoding['chunksizes'])

        longitude = dataset.variables["longitude"]
        test_case.assertEqual((height, width), longitude.shape)
        test_case.assertTrue(np.isnan(longitude.data[2, 2]))
        test_case.assertEqual("longitude", longitude.attrs["standard_name"])
        test_case.assertEqual("degrees_east", longitude.attrs["units"])
        test_case.assertEqual(np.int16, longitude.encoding['dtype'])
        test_case.assertEqual(-32768, longitude.encoding['_FillValue'])
        test_case.assertEqual(0.0054933317, longitude.encoding['scale_factor'])
        test_case.assertEqual(0.0, longitude.encoding['add_offset'])
        if chunking is not None:
            test_case.assertEqual(chunking, longitude.encoding['chunksizes'])

    @staticmethod
    def assert_gridded_geolocation_variables(test_case, dataset, width, height):
        lat = dataset.variables["lat"]
        test_case.assertEqual((height,), lat.shape)
        test_case.assertTrue(np.isnan(lat.data[1]))
        test_case.assertTrue(np.isnan(lat.attrs['_FillValue']))
        test_case.assertEqual("latitude", lat.attrs["standard_name"])
        test_case.assertEqual("latitude", lat.attrs["long_name"])
        test_case.assertEqual("degrees_north", lat.attrs["units"])
        test_case.assertEqual("lat_bnds", lat.attrs["bounds"])

        lat_bnds = dataset.variables["lat_bnds"]
        test_case.assertEqual((height, 2), lat_bnds.shape)
        test_case.assertTrue(np.isnan(lat_bnds.data[2, 0]))
        test_case.assertTrue(np.isnan(lat_bnds.attrs['_FillValue']))
        test_case.assertEqual("latitude cell boundaries", lat_bnds.attrs["long_name"])
        test_case.assertEqual("degrees_north", lat_bnds.attrs["units"])

        lon = dataset.variables["lon"]
        test_case.assertEqual((width,), lon.shape)
        test_case.assertTrue(np.isnan(lon.data[3]))
        test_case.assertTrue(np.isnan(lon.attrs['_FillValue']))
        test_case.assertEqual("longitude", lon.attrs["standard_name"])
        test_case.assertEqual("longitude", lon.attrs["long_name"])
        test_case.assertEqual("degrees_east", lon.attrs["units"])
        test_case.assertEqual("lon_bnds", lon.attrs["bounds"])

        lon_bnds = dataset.variables["lon_bnds"]
        test_case.assertEqual((width, 2), lon_bnds.shape)
        test_case.assertTrue(np.isnan(lon_bnds.data[3, 1]))
        test_case.assertTrue(np.isnan(lon_bnds.attrs['_FillValue']))
        test_case.assertEqual("longitude cell boundaries", lon_bnds.attrs["long_name"])
        test_case.assertEqual("degrees_east", lon_bnds.attrs["units"])

    @staticmethod
    def assert_quality_flags(test_case, dataset, width, height, chunking=None, masks_append=None, meanings_append=None):
        quality = dataset.variables["quality_pixel_bitmask"]
        test_case.assertEqual((height, width), quality.shape)
        test_case.assertEqual(0, quality.data[1, 1])
        test_case.assertEqual("status_flag", quality.attrs["standard_name"])
        test_case.assertEqual("longitude latitude", quality.attrs["coordinates"])

        masks = "1, 2, 4, 8, 16, 32, 64, 128"
        if masks_append is not None:
            masks = masks + masks_append
        test_case.assertEqual(masks, quality.attrs["flag_masks"])

        meanings = "invalid use_with_caution invalid_input invalid_geoloc invalid_time sensor_error padded_data incomplete_channel_data"
        if meanings_append is not None:
            meanings = meanings + meanings_append
        test_case.assertEqual(meanings, quality.attrs["flag_meanings"])

        if chunking is not None:
            test_case.assertEqual(chunking, quality.encoding['chunksizes'])

    @staticmethod
    def assert_correlation_matrices(test_case, ds, num_channels):
        ch_corr_indep = ds.variables["channel_correlation_matrix_independent"]
        test_case.assertEqual((num_channels, num_channels), ch_corr_indep.shape)
        test_case.assertAlmostEqual(1.0, ch_corr_indep[0, 0])
        test_case.assertAlmostEqual(0.0, ch_corr_indep[1, 0])
        test_case.assertEqual(np.int16, ch_corr_indep.encoding['dtype'])
        test_case.assertEqual(-32768, ch_corr_indep.encoding['_FillValue'])
        test_case.assertEqual(0.0001, ch_corr_indep.encoding['scale_factor'])
        test_case.assertEqual(0.0, ch_corr_indep.encoding['add_offset'])
        test_case.assertEqual("Channel_correlation_matrix_independent_effects", ch_corr_indep.attrs["long_name"])
        test_case.assertEqual("1", ch_corr_indep.attrs["units"])
        test_case.assertEqual("-10000", ch_corr_indep.attrs["valid_min"])
        test_case.assertEqual("10000", ch_corr_indep.attrs["valid_max"])
        test_case.assertEqual("Channel error correlation matrix for independent effects", ch_corr_indep.attrs["description"])

        ch_corr_stuct = ds.variables["channel_correlation_matrix_structured"]
        test_case.assertEqual((num_channels, num_channels), ch_corr_stuct.shape)
        test_case.assertAlmostEqual(1.0, ch_corr_stuct[1, 1])
        test_case.assertAlmostEqual(0.0, ch_corr_stuct[0, 2])
        test_case.assertEqual(np.int16, ch_corr_indep.encoding['dtype'])
        test_case.assertEqual(-32768, ch_corr_indep.encoding['_FillValue'])
        test_case.assertEqual(0.0001, ch_corr_indep.encoding['scale_factor'])
        test_case.assertEqual(0.0, ch_corr_indep.encoding['add_offset'])
        test_case.assertEqual("Channel_correlation_matrix_structured_effects", ch_corr_stuct.attrs["long_name"])
        test_case.assertEqual("1", ch_corr_stuct.attrs["units"])
        test_case.assertEqual("-10000", ch_corr_indep.attrs["valid_min"])
        test_case.assertEqual("10000", ch_corr_indep.attrs["valid_max"])
        test_case.assertEqual("Channel error correlation matrix for structured effects", ch_corr_stuct.attrs["description"])

    @staticmethod
    def assert_lookup_tables(test_case, ds, num_channels, lut_size):
        lut_bt = ds.variables["lookup_table_BT"]
        test_case.assertEqual((lut_size, num_channels), lut_bt.shape)
        test_case.assertTrue(np.isnan(lut_bt[1, 2]))
        test_case.assertTrue(np.isnan(lut_bt.attrs["_FillValue"]))
        test_case.assertEqual("Lookup table to convert radiance to brightness temperatures", lut_bt.attrs["description"])

        lut_rad = ds.variables["lookup_table_radiance"]
        test_case.assertEqual((lut_size, num_channels), lut_rad.shape)
        test_case.assertTrue(np.isnan(lut_rad[2, 0]))
        test_case.assertTrue(np.isnan(lut_rad.attrs["_FillValue"]))
        test_case.assertEqual("Lookup table to convert brightness temperatures to radiance", lut_rad.attrs["description"])

    @staticmethod
    def assert_correlation_coefficients(test_case, ds, num_channels, delta_x, delta_y):
        x_element_corr = ds.variables["cross_element_correlation_coefficients"]
        test_case.assertEqual((delta_x, num_channels), x_element_corr.shape)
        test_case.assertTrue(np.isnan(x_element_corr[2, 2]))
        test_case.assertTrue(np.isnan(x_element_corr.attrs["_FillValue"]))
        test_case.assertEqual("cross_element_correlation_coefficients", x_element_corr.attrs["long_name"])
        test_case.assertEqual("Correlation coefficients per channel for scanline correlation", x_element_corr.attrs["description"])

        x_line_corr = ds.variables["cross_line_correlation_coefficients"]
        test_case.assertEqual((delta_y, num_channels), x_line_corr.shape)
        test_case.assertTrue(np.isnan(x_line_corr[3, 0]))
        test_case.assertTrue(np.isnan(x_line_corr.attrs["_FillValue"]))
        test_case.assertEqual("cross_line_correlation_coefficients", x_line_corr.attrs["long_name"])
        test_case.assertEqual("Correlation coefficients per channel for inter scanline correlation", x_line_corr.attrs["description"])

    @staticmethod
    def assert_global_attributes(test_case, attributes):
        test_case.assertIsNotNone(attributes)
        test_case.assertEqual("CF-1.6", attributes["Conventions"])
        test_case.assertEqual("This dataset is released for use under CC-BY licence (https://creativecommons.org/licenses/by/4.0/) and was developed in the EC "
                              "FIDUCEO project \"Fidelity and Uncertainty in Climate Data Records from Earth "
                              "Observations\". Grant Agreement: 638822.", attributes["licence"])
        test_case.assertEqual("2.0.0", attributes["writer_version"])

        test_case.assertIsNone(attributes["institution"])
        test_case.assertIsNone(attributes["source"])
        test_case.assertIsNone(attributes["title"])
        test_case.assertIsNone(attributes["history"])
        test_case.assertIsNone(attributes["references"])  # @todo 3 tb/tb add "id" attribute when we have the DOI issue resolved 2018-06-27  # test_case.assertIsNone(attributes["id"])

    @staticmethod
    def assert_cdr_global_attributes(test_case, attributes):
        test_case.assertIsNone(attributes["source"])
        test_case.assertIsNone(attributes["auxiliary_data"])
        test_case.assertIsNone(attributes["configuration"])
        test_case.assertIsNone(attributes["time_coverage_start"])
        test_case.assertIsNone(attributes["time_coverage_end"])
        test_case.assertIsNone(attributes["time_coverage_duration"])
        test_case.assertIsNone(attributes["time_coverage_resolution"])

    @staticmethod
    def assert_gridded_global_attributes(test_case, attributes):
        test_case.assertIsNone(attributes["geospatial_lat_units"])
        test_case.assertIsNone(attributes["geospatial_lon_units"])
        test_case.assertIsNone(attributes["geospatial_lat_resolution"])
        test_case.assertIsNone(attributes["geospatial_lon_resolution"])
