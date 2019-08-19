import unittest

import xarray as xr

from fiduceo.common.test.assertions import Assertions
from fiduceo.common.writer.writer_utils import WriterUtils


class WriterUtilsTest(unittest.TestCase):

    def test_add_standard_global_attributes(self):
        dataset = xr.Dataset()

        WriterUtils.add_standard_global_attributes(dataset)
        Assertions.assert_global_attributes(self, dataset.attrs)

    def test_add_cdr_global_attributes(self):
        dataset = xr.Dataset()

        WriterUtils.add_cdr_global_attributes(dataset)
        Assertions.assert_cdr_global_attributes(self, dataset.attrs)

    def test_add_gridded_global_attributes(self):
        dataset = xr.Dataset()

        WriterUtils.add_gridded_global_attributes(dataset)
        Assertions.assert_gridded_global_attributes(self, dataset.attrs)
