import numpy as np
from xarray import Variable

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu
from fiduceo.common.writer.writer_utils import WriterUtils


class UTH:

    @staticmethod
    def add_variables(dataset, width, height):
        WriterUtils.add_gridded_global_attributes(dataset)

        tu.add_gridded_geolocation_variables(dataset, width, height)
        tu.add_quality_flags(dataset, width, height)

        dataset["time_ranges_ascend"] = UTH._create_time_ranges_variable(height, width, "Minimum and maximum seconds of day pixel contribution time, ascending nodes")
        dataset["time_ranges_descend"] = UTH._create_time_ranges_variable(height, width, "Minimum and maximum seconds of day pixel contribution time, descending nodes")

        dataset["observation_count_ascend"] = UTH._create_observation_counts_variable(height, width, "Number of UTH/brightness temperature observations in a grid box for ascending passes")
        dataset["observation_count_descend"] = UTH._create_observation_counts_variable(height, width, "Number of UTH/brightness temperature observations in a grid box for descending passes")

        dataset["overpass_count_ascend"] = UTH._create_overpass_counts_variable(height, width, "Number of satellite overpasses in a grid box for ascending passes")
        dataset["overpass_count_descend"] = UTH._create_overpass_counts_variable(height, width, "Number of satellite overpasses in a grid box for descending passes")

        dataset["uth_ascend"] = UTH._create_uth_variable(width, height, description="Monthly average of all UTH retrievals in a grid box for ascending passes (calculated from daily averages)",)
        dataset["uth_descend"] = UTH._create_uth_variable(width, height, description="Monthly average of all UTH retrievals in a grid box for descending passes (calculated from daily averages)")

        dataset["u_independent_uth_ascend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of UTH due to independent effects for ascending passes", coordinates="lon lat")
        dataset["u_independent_uth_descend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of UTH due to independent effects for descending passes", coordinates="lon lat")
        dataset["u_structured_uth_ascend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of UTH due to structured effects for ascending passes", coordinates="lon lat")
        dataset["u_structured_uth_descend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of UTH due to structured effects for descending passes", coordinates="lon lat")
        dataset["u_common_uth_ascend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of UTH due to common effects for ascending passes", coordinates="lon lat")
        dataset["u_common_uth_descend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of UTH due to common effects for descending passes", coordinates="lon lat")

        dataset["uth_inhomogeneity_ascend"] = tu.create_CDR_uncertainty(width, height, "Standard deviation of all daily UTH averages which were used to calculate the monthly UTH average in a grid box for ascending passes", coordinates="lon lat")
        dataset["uth_inhomogeneity_descend"] = tu.create_CDR_uncertainty(width, height, "Standard deviation of all daily UTH averages which were used to calculate the monthly UTH average in a grid box for descending passes", coordinates="lon lat")

        dataset["BT_ascend"] = UTH._create_bt_variable(width, height, description="Monthly average of all brightness temperatures which were used to retrieve UTH in a grid box for ascending passes (calculated from daily averages)")
        dataset["BT_descend"] = UTH._create_bt_variable(width, height, description="Monthly average of all brightness temperatures which were used to retrieve UTH in a grid box for descending passes (calculated from daily averages)")

        dataset["u_independent_BT_ascend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of brightness temperature due to independent effects for ascending passes", coordinates="lon lat",
                                                                       units="K")
        dataset["u_independent_BT_descend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of brightness temperature due to independent effects for descending passes", coordinates="lon lat",
                                                                        units="K")

        dataset["u_structured_BT_ascend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of brightness temperature due to structured effects for ascending passes", coordinates="lon lat",
                                                                      units="K")
        dataset["u_structured_BT_descend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of brightness temperature due to structured effects for descending passes", coordinates="lon lat",
                                                                       units="K")

        dataset["u_common_BT_ascend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of brightness temperature due to common effects for ascending passes", coordinates="lon lat",
                                                                  units="K")
        dataset["u_common_BT_descend"] = tu.create_CDR_uncertainty(width, height, "Uncertainty of brightness temperature due to common effects for descending passes", coordinates="lon lat",
                                                                   units="K")

        dataset["BT_inhomogeneity_ascend"] = tu.create_CDR_uncertainty(width, height,
                                                                       "Standard deviation of all daily brightness temperature averages which were used to calculate the monthly brightness temperature average for ascending passes",
                                                                       coordinates="lon lat", units="K")
        dataset["BT_inhomogeneity_descend"] = tu.create_CDR_uncertainty(width, height,
                                                                       "Standard deviation of all daily brightness temperature averages which were used to calculate the monthly brightness temperature average for descending passes",
                                                                       coordinates="lon lat", units="K")

        dataset["observation_count_all_ascend"] = UTH._create_observation_counts_variable(height, width, "Number of all observations in a grid box for ascending passes - no filtering done")
        dataset["observation_count_all_descend"] = UTH._create_observation_counts_variable(height, width, "Number of all observations in a grid box for descending passes - no filtering done")

    @staticmethod
    def _create_uth_variable(width, height, description=None):
        default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["coordinates"] = "lon lat"
        variable.attrs["long_name"] = "upper_tropospheric_humidity"
        tu.add_units(variable, "%")
        if description is not None:
            variable.attrs["description"] = description
        return variable

    @staticmethod
    def _create_bt_variable(width, height, description=None):
        default_array = DefaultData.create_default_array(width, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["coordinates"] = "lon lat"
        variable.attrs["standard_name"] = "toa_brightness_temperature"
        tu.add_units(variable, "K")

        if description is not None:
            variable.attrs["description"] = description

        return variable

    @staticmethod
    def _create_observation_counts_variable(height, width, description):
        fill_value = DefaultData.get_default_fill_value(np.int16)
        default_array = DefaultData.create_default_array(width, height, np.int16, fill_value=fill_value)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, fill_value)
        variable.attrs["description"] = description
        variable.attrs["coordinates"] = "lon lat"
        return variable

    @staticmethod
    def _create_overpass_counts_variable(height, width, description):
        fill_value = DefaultData.get_default_fill_value(np.uint8)
        default_array = DefaultData.create_default_array(width, height, np.uint8, fill_value=fill_value)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, fill_value)
        variable.attrs["description"] = description
        variable.attrs["coordinates"] = "lon lat"
        return variable

    @staticmethod
    def _create_time_ranges_variable(height, width, description):
        default_array = DefaultData.create_default_array_3d(width, height, 2, np.int32, fill_value=4294967295)
        variable = Variable(["bounds", "y", "x"], default_array)
        tu.add_fill_value(variable, 4294967295)
        tu.add_units(variable, "s")
        variable.attrs["description"] = description
        variable.attrs["coordinates"] = "lon lat"
        return variable
