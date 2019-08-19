import numpy as np

from fiduceo.fcdr.writer.global_flags import GlobalFlags as gf
from fiduceo.fcdr.writer.templates.default_flag_mapper import DefaultFlagMapper


class AVHRR_FlagMapper(DefaultFlagMapper):
    BAD_GEOLOCATION_TIMING_ERR = np.uint8(1)
    BAD_CALIBRATION_RADIOMETER_ERR = np.uint8(2)

    source_masks = [BAD_GEOLOCATION_TIMING_ERR, BAD_CALIBRATION_RADIOMETER_ERR]
    target_masks = [gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION]

    def map_global_flags(self, dataset):
        global_flag_data = dataset["quality_pixel_bitmask"].data
        avhrr_flag_data = dataset["data_quality_bitmask"].data

        global_flag_data = self.evaluate_masks_uint8(avhrr_flag_data, global_flag_data, self.source_masks, self.target_masks)

        dataset["quality_pixel_bitmask"].data = global_flag_data