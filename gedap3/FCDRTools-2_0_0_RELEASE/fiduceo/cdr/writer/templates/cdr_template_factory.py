from fiduceo.cdr.writer.templates.albedo import Albedo
from fiduceo.cdr.writer.templates.aot import AOT
from fiduceo.cdr.writer.templates.sst import SST
from fiduceo.cdr.writer.templates.sst_ensemble import SST_ENSEMBLE
from fiduceo.cdr.writer.templates.uth import UTH


class CDR_TemplateFactory:

    def __init__(self):
        self.templates = dict([("ALBEDO", Albedo), ("AOT", AOT), ("SST", SST), ("SST_ENSEMBLE", SST_ENSEMBLE), ("UTH", UTH)])

    def get_cdr_template(self, name):
        return self.templates[name]
