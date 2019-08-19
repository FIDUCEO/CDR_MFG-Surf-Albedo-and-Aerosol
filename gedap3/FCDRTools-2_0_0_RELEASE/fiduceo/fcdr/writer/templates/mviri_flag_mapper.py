import numpy as np

from fiduceo.fcdr.writer.global_flags import GlobalFlags as gf
from fiduceo.fcdr.writer.templates.default_flag_mapper import DefaultFlagMapper


class MVIRI_FlagMapper(DefaultFlagMapper):
    UNCERTAINTY_SUSPICIOUS = np.uint8(1)
    UNCERTAINTY_TOO_LARGE = np.uint8(2)
    SPACE_VIEW_SUSPICIOUS = np.uint8(4)
    NOT_ON_EARTH = np.uint8(8)
    SUSPECT_TIME = np.uint8(16)
    SUSPECT_GEOLOCATION = np.uint8(32)

    source_masks = [UNCERTAINTY_SUSPICIOUS, UNCERTAINTY_TOO_LARGE, SPACE_VIEW_SUSPICIOUS, NOT_ON_EARTH, SUSPECT_TIME, SUSPECT_GEOLOCATION]
    target_masks = [gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION, gf.INVALID, np.bitwise_or(gf.USE_WITH_CAUTION, gf.INVALID_TIME), np.bitwise_or(gf.USE_WITH_CAUTION, gf.INVALID_GEOLOC)]

    def map_global_flags(self, dataset):
        global_flag_data = dataset["quality_pixel_bitmask"].data
        hirs_flag_data = dataset["data_quality_bitmask"].data

        global_flag_data = self.evaluate_masks_uint8(hirs_flag_data, global_flag_data, self.source_masks, self.target_masks)

        dataset["quality_pixel_bitmask"].data = global_flag_data
