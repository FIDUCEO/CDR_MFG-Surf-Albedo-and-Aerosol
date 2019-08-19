import numpy as np

class GlobalFlags:

    INVALID = np.uint8(1)
    USE_WITH_CAUTION = np.uint8(2)
    INVALID_INPUT = np.uint8(4)
    INVALID_GEOLOC = np.uint8(8)
    INVALID_TIME = np.uint8(16)
    SENSOR_ERROR = np.uint8(32)
    PADDED_DATA = np.uint8(64)
    INCOMPLETE_CHANNEL_DATA = np.uint8(128)