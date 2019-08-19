import numpy as np
from xarray import Variable, Coordinate

from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu
from fiduceo.fcdr.writer.correlation import Correlation as corr

CHUNKSIZES = (500, 500)

SRF_VIS_DIMENSION = "srf_size"
SRF_IR_WV_DIMENSION = "srf_size_ir_wv"
IR_X_DIMENSION = "x_ir_wv"
IR_Y_DIMENSION = "y_ir_wv"

FULL_SIZE = 5000
IR_SIZE = 2500
TIE_SIZE = 500
COV_SIZE = 3

SRF_SIZE = 1011
SOL_IRR_SIZE = 24

NUM_CHANNELS = 3

TIME_FILL_VALUE = -32768


class MVIRI:

    @staticmethod
    def add_original_variables(dataset, height, srf_size=None):
        # height is ignored - supplied just for interface compatibility tb 2017-02-05

        tu.add_quality_flags(dataset, FULL_SIZE, FULL_SIZE, chunksizes=CHUNKSIZES)

        # time
        default_array = DefaultData.create_default_array(IR_SIZE, IR_SIZE, np.uint32)
        variable = Variable([IR_Y_DIMENSION, IR_X_DIMENSION], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint32))
        variable.attrs["standard_name"] = "time"
        variable.attrs["long_name"] = "Acquisition time of pixel"
        tu.add_units(variable, "seconds since 1970-01-01 00:00:00")
        tu.add_offset(variable, TIME_FILL_VALUE)
        tu.add_chunking(variable, CHUNKSIZES)
        dataset["time"] = variable

        dataset["solar_azimuth_angle"] = MVIRI._create_angle_variable_int(0.005493164, standard_name="solar_azimuth_angle", unsigned=True)
        dataset["solar_zenith_angle"] = MVIRI._create_angle_variable_int(0.005493248, standard_name="solar_zenith_angle")
        dataset["satellite_azimuth_angle"] = MVIRI._create_angle_variable_int(0.01, standard_name="sensor_azimuth_angle", long_name="sensor_azimuth_angle", unsigned=True)
        dataset["satellite_zenith_angle"] = MVIRI._create_angle_variable_int(0.01, standard_name="platform_zenith_angle", unsigned=True)

        # count_ir
        default_array = DefaultData.create_default_array(IR_SIZE, IR_SIZE, np.uint8)
        variable = Variable([IR_Y_DIMENSION, IR_X_DIMENSION], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint8))
        variable.attrs["long_name"] = "Infrared Image Counts"
        tu.add_units(variable, "count")
        tu.add_chunking(variable, CHUNKSIZES)
        dataset["count_ir"] = variable

        # count_wv
        default_array = DefaultData.create_default_array(IR_SIZE, IR_SIZE, np.uint8)
        variable = Variable([IR_Y_DIMENSION, IR_X_DIMENSION], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint8))
        variable.attrs["long_name"] = "WV Image Counts"
        tu.add_units(variable, "count")
        tu.add_chunking(variable, CHUNKSIZES)
        dataset["count_wv"] = variable

        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.uint8, fill_value=0)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["flag_masks"] = "1, 2, 4, 8, 16, 32"
        variable.attrs["flag_meanings"] = "uncertainty_suspicious uncertainty_too_large space_view_suspicious not_on_earth suspect_time suspect_geo"
        variable.attrs["standard_name"] = "status_flag"
        tu.add_chunking(variable, CHUNKSIZES)
        dataset["data_quality_bitmask"] = variable

        # distance_sun_earth
        dataset["distance_sun_earth"] = tu.create_scalar_float_variable(long_name="Sun-Earth distance", units="au")

        # solar_irradiance_vis
        dataset["solar_irradiance_vis"] = tu.create_scalar_float_variable(standard_name="solar_irradiance_vis", long_name="Solar effective Irradiance", units="W*m-2")

        # u_solar_irradiance_vis
        default_array = np.full([], np.NaN, np.float32)
        variable = Variable([], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Uncertainty in Solar effective Irradiance"
        tu.add_units(variable, "Wm^-2")
        variable.attrs[corr.PIX_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.PIX_CORR_UNIT] = corr.PIXEL
        variable.attrs[corr.PIX_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs[corr.SCAN_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.SCAN_CORR_UNIT] = corr.LINE
        variable.attrs[corr.SCAN_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs[corr.IMG_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.IMG_CORR_UNIT] = corr.DAYS
        variable.attrs[corr.IMG_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs["pdf_shape"] = "rectangle"
        dataset["u_solar_irradiance_vis"] = variable

        if srf_size is None:
            srf_size = SRF_SIZE

        default_array = DefaultData.create_default_array(srf_size, NUM_CHANNELS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "n_frequencies"], default_array)
        variable.attrs["long_name"] = 'Spectral Response Function weights'
        variable.attrs["description"] = 'Per channel: weights for the relative spectral response function'
        tu.add_encoding(variable, np.int16, -32768, 0.000033)
        dataset['SRF_weights'] = variable

        default_array = DefaultData.create_default_array(srf_size, NUM_CHANNELS, np.float32, fill_value=np.NaN)
        variable = Variable(["channel", "n_frequencies"], default_array)
        variable.attrs["long_name"] = 'Spectral Response Function frequencies'
        variable.attrs["description"] = 'Per channel: frequencies for the relative spectral response function'
        tu.add_encoding(variable, np.int32, -2147483648, 0.0001)
        tu.add_units(variable, "nm")
        variable.attrs["source"] = "Filename of SRF"
        variable.attrs["Valid(YYYYDDD)"] = "datestring"
        dataset['SRF_frequencies'] = variable

        # srf covariance_
        default_array = DefaultData.create_default_array(srf_size, srf_size, np.float32, fill_value=np.NaN)
        variable = Variable([SRF_VIS_DIMENSION, SRF_VIS_DIMENSION], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Covariance of the Visible Band Spectral Response Function"
        tu.add_chunking(variable, CHUNKSIZES)
        dataset["covariance_spectral_response_function_vis"] = variable

        # u_srf_ir
        default_array = DefaultData.create_default_vector(srf_size, np.float32, fill_value=np.NaN)
        variable = Variable([SRF_IR_WV_DIMENSION], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Uncertainty in Spectral Response Function for IR channel"
        dataset["u_spectral_response_function_ir"] = variable

        # u_srf_wv
        default_array = DefaultData.create_default_vector(srf_size, np.float32, fill_value=np.NaN)
        variable = Variable([SRF_IR_WV_DIMENSION], default_array)
        tu.add_fill_value(variable, np.NaN)
        variable.attrs["long_name"] = "Uncertainty in Spectral Response Function for WV channel"
        dataset["u_spectral_response_function_wv"] = variable

        dataset["a_ir"] = tu.create_scalar_float_variable(long_name="Calibration parameter a for IR Band", units="mWm^-2sr^-1cm^-1")
        dataset["b_ir"] = tu.create_scalar_float_variable(long_name="Calibration parameter b for IR Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["u_a_ir"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter a for IR Band", units="mWm^-2sr^-1cm^-1")
        dataset["u_b_ir"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter b for IR Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["a_wv"] = tu.create_scalar_float_variable(long_name="Calibration parameter a for WV Band", units="mWm^-2sr^-1cm^-1")
        dataset["b_wv"] = tu.create_scalar_float_variable(long_name="Calibration parameter b for WV Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["u_a_wv"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter a for WV Band", units="mWm^-2sr^-1cm^-1")
        dataset["u_b_wv"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter b for WV Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["bt_a_ir"] = tu.create_scalar_float_variable(long_name="IR Band BT conversion parameter A", units="1")
        dataset["bt_b_ir"] = tu.create_scalar_float_variable(long_name="IR Band BT conversion parameter B", units="1")
        dataset["bt_a_wv"] = tu.create_scalar_float_variable(long_name="WV Band BT conversion parameter A", units="1")
        dataset["bt_b_wv"] = tu.create_scalar_float_variable(long_name="WV Band BT conversion parameter B", units="1")
        dataset["years_since_launch"] = tu.create_scalar_float_variable(long_name="Fractional year since launch of satellite", units="years")

        x_ir_wv_dim = dataset.dims["x_ir_wv"]
        dataset["x_ir_wv"] = Coordinate("x_ir_wv", np.arange(x_ir_wv_dim, dtype=np.uint16))

        y_ir_wv_dim = dataset.dims["y_ir_wv"]
        dataset["y_ir_wv"] = Coordinate("y_ir_wv", np.arange(y_ir_wv_dim, dtype=np.uint16))

        srf_size_dim = dataset.dims["srf_size"]
        dataset["srf_size"] = Coordinate("srf_size", np.arange(srf_size_dim, dtype=np.uint16))

    @staticmethod
    def get_swath_width():
        return FULL_SIZE

    @staticmethod
    def add_easy_fcdr_variables(dataset, height, corr_dx=None, corr_dy=None, lut_size=None):
        # height is ignored - supplied just for interface compatibility tb 2017-02-05

        # reflectance
        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["standard_name"] = "toa_bidirectional_reflectance_vis"
        variable.attrs["long_name"] = "top of atmosphere bidirectional reflectance factor per pixel of the visible band with central wavelength 0.7"
        tu.add_units(variable, "1")
        tu.add_encoding(variable, np.uint16, DefaultData.get_default_fill_value(np.uint16), 3.05176E-05, chunksizes=CHUNKSIZES)
        dataset["toa_bidirectional_reflectance_vis"] = variable

        # u_independent
        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["long_name"] = "independent uncertainty per pixel"
        tu.add_units(variable, "1")
        tu.add_encoding(variable, np.uint16, DefaultData.get_default_fill_value(np.uint16), 3.05176E-05, chunksizes=CHUNKSIZES)
        dataset["u_independent_toa_bidirectional_reflectance"] = variable

        # u_structured
        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable(["y", "x"], default_array)
        variable.attrs["long_name"] = "structured uncertainty per pixel"
        tu.add_units(variable, "1")
        tu.add_encoding(variable, np.uint16, DefaultData.get_default_fill_value(np.uint16), 3.05176E-05, chunksizes=CHUNKSIZES)
        dataset["u_structured_toa_bidirectional_reflectance"] = variable

        # u_common
        dataset["u_common_toa_bidirectional_reflectance"] = tu.create_scalar_float_variable(long_name="common uncertainty per slot", units="1")

        dataset["sub_satellite_latitude_start"] = tu.create_scalar_float_variable(long_name="Latitude of the sub satellite point at image start", units="degrees_north")
        dataset["sub_satellite_longitude_start"] = tu.create_scalar_float_variable(long_name="Longitude of the sub satellite point at image start", units="degrees_east")
        dataset["sub_satellite_latitude_end"] = tu.create_scalar_float_variable(long_name="Latitude of the sub satellite point at image end", units="degrees_north")
        dataset["sub_satellite_longitude_end"] = tu.create_scalar_float_variable(long_name="Longitude of the sub satellite point at image end", units="degrees_east")

        tu.add_correlation_matrices(dataset, NUM_CHANNELS)

        if lut_size is not None:
            tu.add_lookup_tables(dataset, NUM_CHANNELS, lut_size=lut_size)

        if corr_dx is not None and corr_dy is not None:
            tu.add_correlation_coefficients(dataset, NUM_CHANNELS, corr_dx, corr_dy)

        tu.add_coordinates(dataset, ["vis", "wv", "ir"])

    @staticmethod
    def add_full_fcdr_variables(dataset, height):
        # height is ignored - supplied just for interface compatibility tb 2017-02-05

        # count_vis
        default_array = DefaultData.create_default_array(FULL_SIZE, FULL_SIZE, np.uint8)
        variable = Variable(["y", "x"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint8))
        variable.attrs["long_name"] = "Image counts"
        tu.add_units(variable, "count")
        tu.add_chunking(variable, CHUNKSIZES)
        dataset["count_vis"] = variable

        dataset["u_latitude"] = MVIRI._create_angle_variable_int(1.5E-05, long_name="Uncertainty in Latitude", unsigned=True)
        MVIRI._add_geo_correlation_attributes(dataset["u_latitude"])

        dataset["u_longitude"] = MVIRI._create_angle_variable_int(1.5E-05, long_name="Uncertainty in Longitude", unsigned=True)
        MVIRI._add_geo_correlation_attributes(dataset["u_longitude"])

        # u_time
        default_array = DefaultData.create_default_vector(IR_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable([IR_Y_DIMENSION], default_array)
        variable.attrs["standard_name"] = "Uncertainty in Time"
        tu.add_units(variable, "s")
        tu.add_encoding(variable, np.uint16, DefaultData.get_default_fill_value(np.uint16), 0.009155273)
        variable.attrs["pdf_shape"] = "rectangle"
        dataset["u_time"] = variable

        dataset["u_satellite_zenith_angle"] = MVIRI._create_angle_variable_int(7.62939E-05, long_name="Uncertainty in Satellite Zenith Angle", unsigned=True)
        dataset["u_satellite_azimuth_angle"] = MVIRI._create_angle_variable_int(7.62939E-05, long_name="Uncertainty in Satellite Azimuth Angle", unsigned=True)
        dataset["u_solar_zenith_angle"] = MVIRI._create_angle_variable_int(7.62939E-05, long_name="Uncertainty in Solar Zenith Angle", unsigned=True)
        dataset["u_solar_azimuth_angle"] = MVIRI._create_angle_variable_int(7.62939E-05, long_name="Uncertainty in Solar Azimuth Angle", unsigned=True)

        dataset["a0_vis"] = tu.create_scalar_float_variable("Calibration Coefficient at Launch", units="Wm^-2sr^-1/count")
        dataset["a1_vis"] = tu.create_scalar_float_variable("Time variation of a0", units="Wm^-2sr^-1/count day^-1 10^5")
        dataset["a2_vis"] = tu.create_scalar_float_variable("Time variation of a0, quadratic term", units="Wm^-2sr^-1/count year^-2")
        dataset["mean_count_space_vis"] = tu.create_scalar_float_variable("Space count", units="count")

        # u_a0_vis
        variable = tu.create_scalar_float_variable("Uncertainty in a0", units="Wm^-2sr^-1/count")
        MVIRI._add_calibration_coeff_correlation_attributes(variable)
        dataset["u_a0_vis"] = variable

        # u_a1_vis
        variable = tu.create_scalar_float_variable("Uncertainty in a1", units="Wm^-2sr^-1/count day^-1 10^5")
        MVIRI._add_calibration_coeff_correlation_attributes(variable)
        dataset["u_a1_vis"] = variable

        # u_a2_vis
        variable = tu.create_scalar_float_variable("Uncertainty in a2", units="Wm^-2sr^-1/count year^-2")
        MVIRI._add_calibration_coeff_correlation_attributes(variable)
        dataset["u_a2_vis"] = variable

        # u_zero_vis
        variable = tu.create_scalar_float_variable("Uncertainty zero term", units="Wm^-2sr^-1/count")
        MVIRI._add_calibration_coeff_correlation_attributes(variable, image_correlation_scale=[-np.inf, np.inf])
        dataset["u_zero_vis"] = variable

        # covariance_a_vis
        variable = tu.create_float_variable(COV_SIZE, COV_SIZE, long_name="Covariance of calibration coefficients from fit to calibration runs", dim_names=["cov_size", "cov_size"], fill_value=np.NaN)
        tu.add_fill_value(variable, np.NaN)
        tu.add_units(variable, "Wm^-2sr^-1/count")
        MVIRI._add_calibration_coeff_correlation_attributes(variable, image_correlation_scale=[-np.inf, np.inf])
        dataset["covariance_a_vis"] = variable

        dataset["u_electronics_counts_vis"] = tu.create_scalar_float_variable("Uncertainty due to Electronics noise", units="count")
        dataset["u_digitization_counts_vis"] = tu.create_scalar_float_variable("Uncertainty due to digitization", units="count")

        # allan_deviation_counts_space_vis
        variable = tu.create_scalar_float_variable("Uncertainty of space count", units="count")
        variable.attrs[corr.SCAN_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.SCAN_CORR_UNIT] = corr.LINE
        variable.attrs[corr.SCAN_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs["pdf_shape"] = "digitised_gaussian"
        dataset["allan_deviation_counts_space_vis"] = variable

        # u_mean_counts_space_vis
        variable = tu.create_scalar_float_variable("Uncertainty of space count", units="count")
        variable.attrs[corr.PIX_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.PIX_CORR_UNIT] = corr.PIXEL
        variable.attrs[corr.PIX_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs[corr.SCAN_CORR_FORM] = corr.RECT_ABS
        variable.attrs[corr.SCAN_CORR_UNIT] = corr.LINE
        variable.attrs[corr.SCAN_CORR_SCALE] = [-np.inf, np.inf]
        variable.attrs["pdf_shape"] = "digitised_gaussian"
        dataset["u_mean_counts_space_vis"] = variable

        # sensitivity_solar_irradiance_vis
        variable = tu.create_scalar_float_variable()
        variable.attrs["virtual"] = "true"
        variable.attrs["dimension"] = "y, x"
        variable.attrs[
            "expression"] = "distance_sun_earth * distance_sun_earth * PI * (count_vis - mean_count_space_vis) * (a2_vis * years_since_launch * years_since_launch + a1_vis * years_since_launch + a0_vis) / (cos(solar_zenith_angle * PI / 180.0) * solar_irradiance_vis * solar_irradiance_vis)"
        dataset["sensitivity_solar_irradiance_vis"] = variable

        # sensitivity_count_vis
        variable = tu.create_scalar_float_variable()
        variable.attrs["virtual"] = "true"
        variable.attrs["dimension"] = "y, x"
        variable.attrs[
            "expression"] = "distance_sun_earth * distance_sun_earth * PI * (a2_vis * years_since_launch * years_since_launch + a1_vis * years_since_launch + a0_vis) / (cos(solar_zenith_angle * PI / 180.0) * solar_irradiance_vis)"
        dataset["sensitivity_count_vis"] = variable

        # sensitivity_count_space
        variable = tu.create_scalar_float_variable()
        variable.attrs["virtual"] = "true"
        variable.attrs["dimension"] = "y, x"
        variable.attrs[
            "expression"] = "-1.0 * distance_sun_earth * distance_sun_earth * PI * (a2_vis * years_since_launch * years_since_launch + a1_vis * years_since_launch + a0_vis) / (cos(solar_zenith_angle * PI / 180.0) * solar_irradiance_vis)"
        dataset["sensitivity_count_space"] = variable

        # sensitivity_a0_vis
        variable = tu.create_scalar_float_variable()
        variable.attrs["virtual"] = "true"
        variable.attrs["dimension"] = "y, x"
        variable.attrs["expression"] = "distance_sun_earth * distance_sun_earth * PI * (count_vis - mean_count_space_vis) / (cos(solar_zenith_angle * PI / 180.0) * solar_irradiance_vis)"
        dataset["sensitivity_a0_vis"] = variable

        # sensitivity_a1_vis
        variable = tu.create_scalar_float_variable()
        variable.attrs["virtual"] = "true"
        variable.attrs["dimension"] = "y, x"
        variable.attrs[
            "expression"] = "distance_sun_earth * distance_sun_earth * PI * (count_vis - mean_count_space_vis) * years_since_launch / (cos(solar_zenith_angle * PI / 180.0) * solar_irradiance_vis)"
        dataset["sensitivity_a1_vis"] = variable

        # sensitivity_a2_vis
        variable = tu.create_scalar_float_variable()
        variable.attrs["virtual"] = "true"
        variable.attrs["dimension"] = "y, x"
        variable.attrs[
            "expression"] = "distance_sun_earth * distance_sun_earth * PI * (count_vis - mean_count_space_vis) * years_since_launch*years_since_launch / (cos(solar_zenith_angle * PI / 180.0) * solar_irradiance_vis)"
        dataset["sensitivity_a2_vis"] = variable

        effect_names = ["u_solar_irradiance_vis", "u_a0_vis", "u_a1_vis", "u_a2_vis", "u_zero_vis", "u_solar_zenith_angle", "u_mean_count_space_vis"]
        dataset["Ne"] = Coordinate("Ne", effect_names)

        num_effects = len(effect_names)
        default_array = DefaultData.create_default_array(num_effects, num_effects, np.float32, fill_value=np.NaN)
        variable = Variable(["Ne", "Ne"], default_array)
        tu.add_encoding(variable, np.int16, -32768, 3.05176E-05)
        variable.attrs["valid_min"] = -1
        variable.attrs["valid_max"] = 1
        variable.attrs["long_name"] = "Channel error correlation matrix for structured effects."
        variable.attrs["description"] = "Matrix_describing correlations between errors of the uncertainty_effects due to spectral response function errors (determined using Monte Carlo approach)"
        dataset["effect_correlation_matrix"] = variable

    @staticmethod
    def add_template_key(dataset):
        dataset.attrs["template_key"] = "MVIRI"

    @staticmethod
    def _add_geo_correlation_attributes(geo_variable):
        geo_variable.attrs[corr.PIX_CORR_FORM] = corr.TRI_REL
        geo_variable.attrs[corr.PIX_CORR_UNIT] = corr.PIXEL
        geo_variable.attrs[corr.PIX_CORR_SCALE] = [-250, 250]
        geo_variable.attrs[corr.SCAN_CORR_FORM] = corr.TRI_REL
        geo_variable.attrs[corr.SCAN_CORR_UNIT] = corr.LINE
        geo_variable.attrs[corr.SCAN_CORR_SCALE] = [-250, 250]
        geo_variable.attrs[corr.IMG_CORR_FORM] = corr.TRI_REL
        geo_variable.attrs[corr.IMG_CORR_UNIT] = corr.IMG
        geo_variable.attrs[corr.IMG_CORR_SCALE] = [-12, 0]
        geo_variable.attrs["pdf_shape"] = "gaussian"

    @staticmethod
    def _add_calibration_coeff_correlation_attributes(coeff_variable, image_correlation_scale=None):
        coeff_variable.attrs[corr.PIX_CORR_FORM] = corr.RECT_ABS
        coeff_variable.attrs[corr.PIX_CORR_UNIT] = corr.PIXEL
        coeff_variable.attrs[corr.PIX_CORR_SCALE] = [-np.inf, np.inf]
        coeff_variable.attrs[corr.SCAN_CORR_FORM] = corr.RECT_ABS
        coeff_variable.attrs[corr.SCAN_CORR_UNIT] = corr.LINE
        coeff_variable.attrs[corr.SCAN_CORR_SCALE] = [-np.inf, np.inf]
        coeff_variable.attrs[corr.IMG_CORR_FORM] = corr.TRI_REL
        coeff_variable.attrs[corr.IMG_CORR_UNIT] = corr.MONTHS
        if image_correlation_scale is None:
            coeff_variable.attrs[corr.IMG_CORR_SCALE] = [-1.5, 1.5]
        else:
            coeff_variable.attrs[corr.IMG_CORR_SCALE] = image_correlation_scale

        coeff_variable.attrs["pdf_shape"] = "gaussian"

    @staticmethod
    def _create_angle_variable_int(scale_factor, standard_name=None, long_name=None, unsigned=False, fill_value=None):
        default_array = DefaultData.create_default_array(TIE_SIZE, TIE_SIZE, np.float32, fill_value=np.NaN)
        variable = Variable(["y_tie", "x_tie"], default_array)

        if unsigned is True:
            data_type = np.uint16
        else:
            data_type = np.int16

        if fill_value is None:
            fill_value = DefaultData.get_default_fill_value(data_type)

        if standard_name is not None:
            variable.attrs["standard_name"] = standard_name

        if long_name is not None:
            variable.attrs["long_name"] = long_name

        tu.add_units(variable, "degree")
        variable.attrs["tie_points"] = "true"
        tu.add_encoding(variable, data_type, fill_value, scale_factor, chunksizes=CHUNKSIZES)
        return variable
