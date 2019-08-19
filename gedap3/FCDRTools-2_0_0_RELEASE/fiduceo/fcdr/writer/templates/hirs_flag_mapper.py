import numpy as np

from fiduceo.fcdr.writer.global_flags import GlobalFlags as gf
from fiduceo.fcdr.writer.templates.default_flag_mapper import DefaultFlagMapper


class HIRS_FlagMapper(DefaultFlagMapper):
    # pixel_quality
    SUSPECT_MIRROR = np.uint8(1)
    SUSPECT_GEO = np.uint8(2)
    SUSPECT_TIME = np.uint8(4)
    OUTLIER_NOS = np.uint8(8)
    UNCERTAINTY_TOO_LARGE = np.uint8(16)

    # scanline_quality
    REDUCED_CONTEXT = np.int32(536870912)
    BAD_TEMP_NO_RSELF = np.int32(1073741824)

    # channel_quality
    DO_NOT_USE = np.uint8(1)
    UNCERTAINTY_SUSPICIOUS = np.uint8(2)
    SELF_EMISSION_FAILS = np.uint8(4)
    CALIBRATION_IMPOSSIBLE = np.uint8(8)
    CALIBRATION_SUSPECT = np.uint8(16)

    source_masks = [SUSPECT_MIRROR, SUSPECT_GEO, SUSPECT_TIME, OUTLIER_NOS, UNCERTAINTY_TOO_LARGE]
    target_masks = [gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION, gf.USE_WITH_CAUTION]

    source_scanline_masks = [REDUCED_CONTEXT, BAD_TEMP_NO_RSELF]
    target_scanline_masks = [gf.USE_WITH_CAUTION, gf.INVALID]

    source_channel_masks_dual = [DO_NOT_USE, SELF_EMISSION_FAILS, CALIBRATION_IMPOSSIBLE]
    source_channel_masks = [UNCERTAINTY_SUSPICIOUS, CALIBRATION_SUSPECT]

    def map_global_flags(self, dataset):
        global_flag_data = dataset["quality_pixel_bitmask"].data
        hirs_flag_data = dataset["data_quality_bitmask"].data

        global_flag_data = self.evaluate_masks_uint8(hirs_flag_data, global_flag_data, self.source_masks, self.target_masks)

        global_flag_data = self.apply_scanline_flags(dataset, global_flag_data)

        global_flag_data = self.apply_channel_flags(dataset, global_flag_data)

        dataset["quality_pixel_bitmask"].data = global_flag_data

    def apply_channel_flags(self, dataset, global_flag_data):
        if not "quality_channel_bitmask" in dataset.data_vars:  # special case for HIRS2, does not contain this variable tb 2018-02-19
            return global_flag_data

        channel_flag_data = dataset["quality_channel_bitmask"].data
        num_channels = channel_flag_data.shape[1]
        for line in range(0, channel_flag_data.shape[0]):
            for source_mask in self.source_channel_masks_dual:
                flag_count = 0
                flag_set = False
                for channel in range(0, num_channels):
                    channel_flag_set = np.bitwise_and(source_mask, channel_flag_data[line, channel]) > 0
                    flag_set |= channel_flag_set
                    flag_count += channel_flag_set.astype(np.uint8)

                if flag_count == num_channels:
                    global_flag_data[line, :] = np.bitwise_or(global_flag_data[line, :], gf.INVALID)
                elif flag_set:
                    global_flag_data[line, :] = np.bitwise_or(global_flag_data[line, :], gf.USE_WITH_CAUTION)

            for source_mask in self.source_channel_masks:
                flag_set = False
                for channel in range(0, num_channels):
                    channel_flag_set = np.bitwise_and(source_mask, channel_flag_data[line, channel]) > 0
                    flag_set |= channel_flag_set

                if flag_set:
                    global_flag_data[line, :] = np.bitwise_or(global_flag_data[line, :], gf.USE_WITH_CAUTION)

        return global_flag_data

    def apply_scanline_flags(self, dataset, global_flag_data):
        scanline_flag_data = dataset["quality_scanline_bitmask"].data
        for line in range(0, scanline_flag_data.shape[0]):
            for source_mask, target_mask in zip(self.source_scanline_masks, self.target_scanline_masks):
                flag_set = np.bitwise_and(source_mask, scanline_flag_data[line]) > 0
                if flag_set:
                    global_flag_data[line, :] = np.bitwise_or(global_flag_data[line, :], target_mask)

        return global_flag_data
