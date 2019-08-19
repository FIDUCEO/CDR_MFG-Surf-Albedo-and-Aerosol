import numpy as np

class DefaultFlagMapper:

    def map_global_flags(self, dataset):
        pass

    def evaluate_masks_uint8(self, avhrr_flag_data, global_flag_data, source_masks, target_masks):
        for source_mask, target_mask in zip(source_masks, target_masks):
            intermediate = np.bitwise_and(avhrr_flag_data, source_mask) > 0
            intermediate = intermediate.astype(np.uint8) * target_mask
            global_flag_data = np.bitwise_or(global_flag_data, intermediate)
        return global_flag_data