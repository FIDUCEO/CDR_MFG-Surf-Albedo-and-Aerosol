import numpy as np
from xarray import Variable

from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu
from fiduceo.fcdr.writer.correlation import Correlation as corr
from fiduceo.common.writer.default_data import DefaultData

SWATH_WIDTH = 409
PRT_WIDTH = 3
N_CHANS = 6
MAX_SRF_SIZE = 5902
CHUNKS_2D = (1280, 409)

COUNT_CORRELATION_ATTRIBUTES = {corr.PIX_CORR_FORM: corr.RECT_ABS, corr.PIX_CORR_UNIT: corr.PIXEL, corr.PIX_CORR_SCALE: [-np.inf, np.inf], corr.SCAN_CORR_FORM: corr.TRI_REL,
                                corr.SCAN_CORR_UNIT: corr.LINE, corr.SCAN_CORR_SCALE: [-25, 25], "pdf_shape": "digitised_gaussian"}


class AVHRR:
    @staticmethod
    def add_original_variables(dataset, height, srf_size=None):
        tu.add_geolocation_variables(dataset, SWATH_WIDTH, height, chunksizes=CHUNKS_2D)
        tu.add_quality_flags(dataset, SWATH_WIDTH, height, chunksizes=CHUNKS_2D)

        # Time
        default_array = DefaultData.create_default_vector(height, np.float64, fill_value=np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        tu.add_units(variable, "s")
        variable.attrs["standard_name"] = "time"
        variable.attrs["long_name"] = "Acquisition time in seconds since 1970-01-01 00:00:00"
        dataset["Time"] = variable

        # relative_azimuth_angle
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "relative_azimuth_angle"
        tu.add_units(variable, "degree")
        tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.01, chunksizes=CHUNKS_2D)
        variable.attrs["valid_max"] = 18000
        variable.attrs["valid_min"] = -18000
        tu.add_geolocation_attribute(variable)
        dataset["relative_azimuth_angle"] = variable

        # satellite_zenith_angle
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "sensor_zenith_angle"
        tu.add_units(variable, "degree")
        tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.01, chunksizes=CHUNKS_2D)
        variable.attrs["valid_max"] = 9000
        variable.attrs["valid_min"] = 0
        tu.add_geolocation_attribute(variable)
        dataset["satellite_zenith_angle"] = variable

        # solar_zenith_angle
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "solar_zenith_angle"
        tu.add_units(variable, "degree")
        tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.01, chunksizes=CHUNKS_2D)
        variable.attrs["valid_max"] = 18000
        variable.attrs["valid_min"] = 0
        tu.add_geolocation_attribute(variable)
        dataset["solar_zenith_angle"] = variable

        dataset["Ch1"] = AVHRR._create_channel_refl_variable(height, "Channel 1 Reflectance")
        dataset["Ch2"] = AVHRR._create_channel_refl_variable(height, "Channel 2 Reflectance")
        dataset["Ch3a"] = AVHRR._create_channel_refl_variable(height, "Channel 3a Reflectance")
        dataset["Ch3b"] = AVHRR._create_channel_bt_variable(height, "Channel 3b Brightness Temperature")
        dataset["Ch4"] = AVHRR._create_channel_bt_variable(height, "Channel 4 Brightness Temperature")
        dataset["Ch5"] = AVHRR._create_channel_bt_variable(height, "Channel 5 Brightness Temperature")

        # data_quality_bitmask
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.uint8, fill_value=0)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = 'status_flag'
        variable.attrs["long_name"] = 'bitmask for quality per pixel'
        variable.attrs["flag_masks"] = '1,2'
        variable.attrs['flag_meanings'] = 'bad_geolocation_timing_err bad_calibration_radiometer_err'
        tu.add_chunking(variable, CHUNKS_2D)
        tu.add_geolocation_attribute(variable)
        dataset['data_quality_bitmask'] = variable

        default_array = DefaultData.create_default_vector(height, np.uint8, fill_value=0)
        variable = Variable(["y"], default_array)
        variable.attrs["long_name"] = 'bitmask for quality per scanline'
        variable.attrs["standard_name"] = 'status_flag'
        variable.attrs["flag_masks"] = '1,2,4,8,16,32,64'
        variable.attrs['flag_meanings'] = 'do_not_use bad_time bad_navigation bad_calibration channel3a_present solar_contamination_failure solar_contamination'
        dataset['quality_scanline_bitmask'] = variable

        default_array = DefaultData.create_default_array(N_CHANS, height, np.uint8, fill_value=0)
        variable = Variable(["y", "channel"], default_array)
        variable.attrs["long_name"] = 'bitmask for quality per channel'
        variable.attrs["standard_name"] = 'status_flag'
        variable.attrs["flag_masks"] = '1,2'
        variable.attrs['flag_meanings'] = 'bad_channel some_pixels_not_detected_2sigma'
        dataset['quality_channel_bitmask'] = variable

        if srf_size is None:
            srf_size = MAX_SRF_SIZE

        default_array = DefaultData.create_default_array(srf_size, N_CHANS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "n_frequencies"], default_array)
        variable.attrs["long_name"] = 'Spectral Response Function weights'
        variable.attrs["description"] = 'Per channel: weights for the relative spectral response function'
        tu.add_encoding(variable, np.int16, -32768, 0.000033)
        dataset['SRF_weights'] = variable

        default_array = DefaultData.create_default_array(srf_size, N_CHANS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "n_frequencies"], default_array)
        variable.attrs["long_name"] = 'Spectral Response Function wavelengths'
        variable.attrs["description"] = 'Per channel: wavelengths for the relative spectral response function'
        tu.add_encoding(variable, np.int32, -2147483648, 0.0001)
        tu.add_units(variable, "um")
        dataset['SRF_wavelengths'] = variable

        default_vector = DefaultData.create_default_vector(height, np.uint8, fill_value=255)
        variable = Variable(["y"], default_vector)
        tu.add_fill_value(variable, 255)
        variable.attrs["long_name"] = 'Indicator of original file'
        variable.attrs["description"] = "Indicator for mapping each line to its corresponding original level 1b file. See global attribute 'source' for the filenames. 0 corresponds to 1st listed file, 1 to 2nd file."
        dataset["scanline_map_to_origl1bfile"] = variable

        default_vector = DefaultData.create_default_vector(height, np.int16, fill_value=DefaultData.get_default_fill_value(np.int16))
        variable = Variable(["y"], default_vector)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int16))
        variable.attrs["long_name"] = 'Original_Scan_line_number'
        variable.attrs["description"] = 'Original scan line numbers from corresponding l1b records'
        dataset["scanline_origl1b"] = variable

        tu.add_coordinates(dataset, ["Ch1", "Ch2", "Ch3a", "Ch3b", "Ch4", "Ch5"])

    @staticmethod
    def get_swath_width():
        return SWATH_WIDTH

    @staticmethod
    def add_easy_fcdr_variables(dataset, height, corr_dx=None, corr_dy=None, lut_size=None):
        # u_independent_Ch1-3a
        long_names = ["independent uncertainty per pixel for channel 1", "independent uncertainty per pixel for channel 2", "independent uncertainty per pixel for channel 3a"]
        names = ["u_independent_Ch1", "u_independent_Ch2", "u_independent_Ch3a"]
        AVHRR._add_refl_uncertainties_variables(dataset, height, names, long_names)

        # u_structured_Ch1-3a
        long_names = ["structured uncertainty per pixel for channel 1", "structured uncertainty per pixel for channel 2", "structured uncertainty per pixel for channel 3a"]
        names = ["u_structured_Ch1", "u_structured_Ch2", "u_structured_Ch3a"]
        AVHRR._add_refl_uncertainties_variables(dataset, height, names, long_names, structured=True)

        # u_common_Ch1-3a
        long_names = ["common uncertainty per pixel for channel 1", "common uncertainty per pixel for channel 2", "common uncertainty per pixel for channel 3a"]
        names = ["u_common_Ch1", "u_common_Ch2", "u_common_Ch3a"]
        AVHRR._add_refl_uncertainties_variables(dataset, height, names, long_names)

        # u_independent_Ch3b-5
        long_names = ["independent uncertainty per pixel for channel 3b", "independent uncertainty per pixel for channel 4", "independent uncertainty per pixel for channel 5"]
        names = ["u_independent_Ch3b", "u_independent_Ch4", "u_independent_Ch5"]
        AVHRR._add_bt_uncertainties_variables(dataset, height, names, long_names)

        # u_structured_Ch3b-5
        long_names = ["structured uncertainty per pixel for channel 3b", "structured uncertainty per pixel for channel 4", "structured uncertainty per pixel for channel 5"]
        names = ["u_structured_Ch3b", "u_structured_Ch4", "u_structured_Ch5"]
        AVHRR._add_bt_uncertainties_variables(dataset, height, names, long_names)

        # u_common_Ch3b-5
        long_names = ["common uncertainty per pixel for channel 3b", "common uncertainty per pixel for channel 4", "common uncertainty per pixel for channel 5"]
        names = ["u_common_Ch3b", "u_common_Ch4", "u_common_Ch5"]
        AVHRR._add_bt_uncertainties_variables(dataset, height, names, long_names)

        tu.add_correlation_matrices(dataset, N_CHANS)

        if lut_size is not None:
            tu.add_lookup_tables(dataset, N_CHANS, lut_size)

        if corr_dx is not None and corr_dy is not None:
            tu.add_correlation_coefficients(dataset, N_CHANS, corr_dx, corr_dy)

    @staticmethod
    def add_full_fcdr_variables(dataset, height):
        # u_latitude
        variable = AVHRR._create_angle_uncertainty_variable("latitude", height)
        dataset["u_latitude"] = variable

        # u_longitude
        variable = AVHRR._create_angle_uncertainty_variable("longitude", height)
        dataset["u_longitude"] = variable

        # u_time
        default_array = DefaultData.create_default_vector(height, np.float64, fill_value=np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        tu.add_units(variable, "s")
        variable.attrs["long_name"] = "uncertainty of acquisition time"
        dataset["u_time"] = variable

        # u_satellite_azimuth_angle
        variable = AVHRR._create_angle_uncertainty_variable("satellite azimuth angle", height)
        dataset["u_satellite_azimuth_angle"] = variable

        # u_satellite_zenith_angle
        variable = AVHRR._create_angle_uncertainty_variable("satellite zenith angle", height)
        dataset["u_satellite_zenith_angle"] = variable

        # u_solar_azimuth_angle
        variable = AVHRR._create_angle_uncertainty_variable("solar azimuth angle", height)
        dataset["u_solar_azimuth_angle"] = variable

        # u_solar_zenith_angle
        variable = AVHRR._create_angle_uncertainty_variable("solar zenith angle", height)
        dataset["u_solar_zenith_angle"] = variable

        # PRT_C
        default_array = DefaultData.create_default_array(PRT_WIDTH, height, np.int16)
        variable = Variable(["y", "n_prt"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int16))
        variable.attrs["long_name"] = "Prt counts"
        tu.add_units(variable, "count")
        dataset["PRT_C"] = variable

        # u_prt
        default_array = DefaultData.create_default_array(PRT_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "n_prt"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Uncertainty on the PRT counts"
        tu.add_units(variable, "count")
        variable.attrs[corr.PIX_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.PIX_CORR_UNIT] = corr.PIXEL
        variable.attrs[corr.PIX_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs[corr.SCAN_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.SCAN_CORR_UNIT] = corr.LINE
        variable.attrs[corr.SCAN_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs["pdf_shape"] = "rectangle"
        variable.attrs["pdf_parameter"] = 0.1
        dataset["u_prt"] = variable

        # R_ICT
        default_array = DefaultData.create_default_array(PRT_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "n_prt"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Radiance of the PRT"
        tu.add_units(variable, "mW m^-2 sr^-1 cm")
        dataset["R_ICT"] = variable

        # T_instr
        default_array = DefaultData.create_default_vector(height, np.float32, fill_value=np.NaN)
        variable = Variable(["y"], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Instrument temperature"
        tu.add_units(variable, "K")
        dataset["T_instr"] = variable

        # Chx_Csp
        standard_names = ["Ch1 Space counts", "Ch2 Space counts", "Ch3a Space counts", "Ch3b Space counts", "Ch4 Space counts", "Ch5 Space counts"]
        names = ["Ch1_Csp", "Ch2_Csp", "Ch3a_Csp", "Ch3b_Csp", "Ch4_Csp", "Ch5_Csp"]
        AVHRR._add_counts_variables(dataset, height, names, standard_names)

        # Chx_Cict
        standard_names = ["Ch3b ICT counts", "Ch4 ICT counts", "Ch5 ICT counts"]
        names = ["Ch3b_Cict", "Ch4_Cict", "Ch5_Cict"]
        AVHRR._add_counts_variables(dataset, height, names, standard_names)

        # Chx_Ce
        standard_names = ["Ch1 Earth counts", "Ch2 Earth counts", "Ch3a Earth counts", "Ch3b Earth counts", "Ch4 Earth counts", "Ch5 Earth counts"]
        names = ["Ch1_Ce", "Ch2_Ce", "Ch3a_Ce", "Ch3b_Ce", "Ch4_Ce", "Ch5_Ce"]
        AVHRR._add_counts_variables(dataset, height, names, standard_names)

        # Chx_u_Csp
        standard_names = ["Ch1 Uncertainty on space counts", "Ch2 Uncertainty on space counts", "Ch3a Uncertainty on space counts", "Ch3b Uncertainty on space counts",
                          "Ch4 Uncertainty on space counts", "Ch5 Uncertainty on space counts"]
        names = ["Ch1_u_Csp", "Ch2_u_Csp", "Ch3a_u_Csp", "Ch3b_u_Csp", "Ch4_u_Csp", "Ch5_u_Csp"]
        AVHRR._add_counts_uncertainties_variables(dataset, height, names, standard_names, COUNT_CORRELATION_ATTRIBUTES)

        # Chx_Cict
        standard_names = ["Ch3b Uncertainty on ICT counts", "Ch4 Uncertainty on ICT counts", "Ch5 Uncertainty on ICT counts"]
        names = ["Ch3b_u_Cict", "Ch4_u_Cict", "Ch5_u_Cict"]
        AVHRR._add_counts_uncertainties_variables(dataset, height, names, standard_names, COUNT_CORRELATION_ATTRIBUTES)

        # Chx_u_Ce
        standard_names = ["Ch1 Uncertainty on earth counts", "Ch2 Uncertainty on earth counts", "Ch3a Uncertainty on earth counts", "Ch3b Uncertainty on earth counts",
                          "Ch4 Uncertainty on earth counts", "Ch5 Uncertainty on earth counts"]
        names = ["Ch1_u_Ce", "Ch2_u_Ce", "Ch3a_u_Ce", "Ch3b_u_Ce", "Ch4_u_Ce", "Ch5_u_Ce"]
        attributes = {"pdf_shape": "digitised_gaussian"}
        AVHRR._add_counts_uncertainties_variables(dataset, height, names, standard_names, attributes)

        # Chx_u_Refl
        long_names = ["Ch1 Total uncertainty on toa reflectance", "Ch2 Total uncertainty on toa reflectance", "Ch3a Total uncertainty on toa reflectance"]
        names = ["Ch1_u_Refl", "Ch2_u_Refl", "Ch3a_u_Refl"]
        AVHRR._add_refl_uncertainties_variables(dataset, height, names, long_names)

        # Chx_u_Bt
        standard_names = ["Ch3b Total uncertainty on brightness temperature", "Ch4 Total uncertainty on brightness temperature", "Ch5 Total uncertainty on brightness temperature"]
        names = ["Ch3b_u_Bt", "Ch4_u_Bt", "Ch5_u_Bt"]
        AVHRR._add_bt_uncertainties_variables(dataset, height, names, standard_names)

        # Chx_ur_Bt
        standard_names = ["Ch3b Random uncertainty on brightness temperature", "Ch4 Random uncertainty on brightness temperature", "Ch5 Random uncertainty on brightness temperature"]
        names = ["Ch3b_ur_Bt", "Ch4_ur_Bt", "Ch5_ur_Bt"]
        AVHRR._add_bt_uncertainties_variables(dataset, height, names, standard_names)

        # Chx_us_Bt
        standard_names = ["Ch3b Systematic uncertainty on brightness temperature", "Ch4 Systematic uncertainty on brightness temperature", "Ch5 Systematic uncertainty on brightness temperature"]
        names = ["Ch3b_us_Bt", "Ch4_us_Bt", "Ch5_us_Bt"]
        AVHRR._add_bt_uncertainties_variables(dataset, height, names, standard_names)

    @staticmethod
    def add_template_key(dataset):
        dataset.attrs["template_key"] = "AVHRR"

    @staticmethod
    def _add_counts_uncertainties_variables(dataset, height, names, long_names, attributes=None):
        for i, name in enumerate(names):
            variable = AVHRR._create_counts_uncertainty_variable(height, long_names[i])
            if attributes is not None:
                for key, value in attributes.items():
                    variable.attrs[key] = value
            dataset[name] = variable

    @staticmethod
    def _add_bt_uncertainties_variables(dataset, height, names, long_names):
        for i, name in enumerate(names):
            variable = AVHRR._create_bt_uncertainty_variable(height, long_name=long_names[i])
            dataset[name] = variable

    @staticmethod
    def _add_refl_uncertainties_variables(dataset, height, names, long_names, structured=False):
        for i, name in enumerate(names):
            variable = AVHRR._create_refl_uncertainty_variable(height, long_name=long_names[i], structured=structured)
            dataset[name] = variable

    @staticmethod
    def _add_counts_variables(dataset, height, names, standard_names):
        for i, name in enumerate(names):
            variable = AVHRR._create_counts_variable(height, standard_names[i])
            dataset[name] = variable

    @staticmethod
    def _create_counts_uncertainty_variable(height, long_name):
        variable = tu.create_float_variable(SWATH_WIDTH, height, long_name=long_name, fill_value=np.NaN)
        tu.add_units(variable, "count")
        tu.add_geolocation_attribute(variable)
        tu.add_chunking(variable, CHUNKS_2D)
        return variable

    @staticmethod
    def _create_counts_variable(height, long_name):
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.int32)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.int32))
        variable.attrs["long_name"] = long_name
        tu.add_units(variable, "count")
        tu.add_geolocation_attribute(variable)
        tu.add_chunking(variable, CHUNKS_2D)
        return variable

    @staticmethod
    def _create_angle_uncertainty_variable(angle_name, height):
        variable = tu.create_float_variable(SWATH_WIDTH, height, long_name="uncertainty of " + angle_name, fill_value=np.NaN)
        tu.add_units(variable, "degree")
        tu.add_geolocation_attribute(variable)
        tu.add_chunking(variable, CHUNKS_2D)
        return variable

    @staticmethod
    def _create_refl_uncertainty_variable(height, long_name=None, structured=False):
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)

        tu.add_units(variable, "percent")
        tu.add_geolocation_attribute(variable)
        variable.attrs["long_name"] = long_name

        if structured:
            tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.01, chunksizes=CHUNKS_2D)
            variable.attrs["valid_min"] = 3
            variable.attrs["valid_max"] = 5
        else:
            tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.00001, chunksizes=CHUNKS_2D)
            variable.attrs["valid_max"] = 1000
            variable.attrs["valid_min"] = 10
        return variable

    @staticmethod
    def _create_bt_uncertainty_variable(height, long_name):
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        tu.add_units(variable, "K")
        tu.add_geolocation_attribute(variable)
        tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.001, chunksizes=CHUNKS_2D)
        variable.attrs["valid_max"] = 15000
        variable.attrs["valid_min"] = 1
        variable.attrs["long_name"] = long_name
        return variable

    @staticmethod
    def _create_channel_refl_variable(height, long_name):
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "toa_reflectance"
        variable.attrs["long_name"] = long_name
        tu.add_units(variable, "1")
        tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.0001, chunksizes=CHUNKS_2D)
        variable.attrs["valid_max"] = 15000
        variable.attrs["valid_min"] = 0
        tu.add_geolocation_attribute(variable)
        return variable

    @staticmethod
    def _create_channel_bt_variable(height, long_name):
        default_array = DefaultData.create_default_array(SWATH_WIDTH, height, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "toa_brightness_temperature"
        variable.attrs["long_name"] = long_name
        tu.add_units(variable, "K")
        variable.attrs["valid_max"] = 10000
        variable.attrs["valid_min"] = -20000
        tu.add_geolocation_attribute(variable)
        tu.add_encoding(variable, np.int16, DefaultData.get_default_fill_value(np.int16), 0.01, 273.15, chunksizes=CHUNKS_2D)
        return variable
