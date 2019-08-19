import os

import xarray as xr

from fiduceo.common.version import __version__
from fiduceo.common.writer.writer_utils import WriterUtils
from fiduceo.fcdr.writer.templates.template_factory import TemplateFactory

DATE_PATTERN = "%Y%m%d%H%M%S"


class FCDRWriter:

    @staticmethod
    def write(ds, file, compression_level=None, overwrite=False):
        """
        Save a dataset to NetCDF file.
        :param ds: The dataset
        :param file: File path
        :param compression_level: the file compression level, 0 - 9, default is 5
        :param overwrite: set true to overwrite existing files
         """
        if os.path.isfile(file):
            if overwrite is True:
                os.remove(file)
            else:
                raise IOError("The file already exists: " + file)

        # trigger mapping of sensor specific flags to the global flag variable
        template_factory = TemplateFactory()
        flag_mapper = template_factory.get_flag_mapper(ds.attrs["template_key"])
        flag_mapper.map_global_flags(ds)

        # set up compression parameter for ALL variables. Unfortunately, xarray does not allow
        # one set of compression params per file, only per variable. tb 2017-01-25
        if compression_level is None:
            compression_level = 5

        comp = dict(zlib=True, complevel=compression_level)
        encoding = dict()
        for var_name in ds.data_vars:
            var_encoding = dict(comp)
            var_encoding.update(ds[var_name].encoding)
            encoding.update({var_name: var_encoding})

        ds.to_netcdf(file, format='netCDF4', engine='netcdf4', encoding=encoding)

    @staticmethod
    def createTemplateEasy(sensorType, height, srf_size=None, corr_dx=None, corr_dy=None, lut_size=None):
        """
        Create a template dataset in EASY FCDR format for the sensor given as argument.
        :param sensorType: the sensor type to create the template for
        :param height: the height in pixels of the data product
        :param srf_size: if set, the length of the spectral response function in frequency steps
        :param corr_dx: correlation length across track
        :param corr_dy: correlation length along track
        :param lut_size: size of a BT/radiance conversion lookup table
        :return the template dataset
         """
        dataset = xr.Dataset()
        WriterUtils.add_standard_global_attributes(dataset)

        template_factory = TemplateFactory()

        sensor_template = template_factory.get_sensor_template(sensorType)
        sensor_template.add_original_variables(dataset, height, srf_size)
        sensor_template.add_easy_fcdr_variables(dataset, height, corr_dx, corr_dy, lut_size)
        sensor_template.add_template_key(dataset)

        return dataset

    @staticmethod
    def createTemplateFull(sensorType, height):
        """
        Create a template dataset in FULL FCDR format for the sensor given as argument.
        :param sensorType: the sensor type to create the template for
        :param height the hheight in pixels of the data product
        :return the template dataset
         """
        dataset = xr.Dataset()

        WriterUtils.add_standard_global_attributes(dataset)

        template_factory = TemplateFactory()

        sensor_template = template_factory.get_sensor_template(sensorType)
        sensor_template.add_original_variables(dataset, height)
        sensor_template.add_full_fcdr_variables(dataset, height)
        sensor_template.add_template_key(dataset)

        return dataset

    @staticmethod
    def create_file_name_FCDR_easy(sensor, platform, start, end, version):
        """
        Create a file name for EASY FCDR format .
        :param sensor: the sensor name
        :param platform: the name of the satellite platform
        :param start: the acquisition start date and time, type datetime
        :param end: the acquisition end date and time, type datetime
        :param version: the processor version string, format "xx.x"
        :return a valid file name
         """
        return FCDRWriter._create_file_name(end, platform, sensor, start, "EASY", version)

    @staticmethod
    def create_file_name_FCDR_full(sensor, platform, start, end, version):
        """
        Create a file name for FULL FCDR format .
        :param sensor: the sensor name
        :param platform: the name of the satellite platform
        :param start: the acquisition start date and time, type datetime
        :param end: the acquisition end date and time, type datetime
        :param version: the processor version string, format "xx.x"
        :return a valid file name
         """
        return FCDRWriter._create_file_name(end, platform, sensor, start, "FULL", version)

    @staticmethod
    def _create_file_name(end, platform, sensor, start, type, version):
        start_string = start.strftime(DATE_PATTERN)
        end_string = end.strftime(DATE_PATTERN)
        return "FIDUCEO_FCDR_L1C_" + sensor + "_" + platform + "_" + start_string + "_" + end_string + "_" + type + "_v" + version + "_fv" + __version__ + ".nc"
