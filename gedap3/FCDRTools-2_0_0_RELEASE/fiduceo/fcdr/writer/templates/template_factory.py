from fiduceo.fcdr.writer.templates.amsub_mhs import AMSUB_MHS
from fiduceo.fcdr.writer.templates.avhrr import AVHRR
from fiduceo.fcdr.writer.templates.avhrr_flag_mapper import AVHRR_FlagMapper
from fiduceo.fcdr.writer.templates.default_flag_mapper import DefaultFlagMapper
from fiduceo.fcdr.writer.templates.hirs_2 import HIRS2
from fiduceo.fcdr.writer.templates.hirs_3 import HIRS3
from fiduceo.fcdr.writer.templates.hirs_4 import HIRS4
from fiduceo.fcdr.writer.templates.hirs_flag_mapper import HIRS_FlagMapper
from fiduceo.fcdr.writer.templates.mviri import MVIRI
from fiduceo.fcdr.writer.templates.mviri_flag_mapper import MVIRI_FlagMapper
from fiduceo.fcdr.writer.templates.mviri_static import MVIRI_STATIC
from fiduceo.fcdr.writer.templates.ssmt2 import SSMT2


class TemplateFactory:
    def __init__(self):
        self.templates = dict(
            [("AMSUB", AMSUB_MHS), ("MHS", AMSUB_MHS), ("SSMT2", SSMT2), ("AVHRR", AVHRR), ("HIRS2", HIRS2), ("HIRS3", HIRS3), ("HIRS4", HIRS4), ("MVIRI", MVIRI), ("MVIRI_STATIC", MVIRI_STATIC)])

        self.flag_mapper = dict(
            [("AMSUB", DefaultFlagMapper()), ("MHS", DefaultFlagMapper()), ("SSMT2", DefaultFlagMapper()), ("AVHRR", AVHRR_FlagMapper()), ("HIRS2", HIRS_FlagMapper()), ("HIRS3", HIRS_FlagMapper()),
             ("HIRS4", HIRS_FlagMapper()), ("MVIRI", MVIRI_FlagMapper())])

    def get_sensor_template(self, name):
        return self.templates[name]

    def get_flag_mapper(self, name):
        return self.flag_mapper[name]
