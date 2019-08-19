import unittest

from fiduceo.cdr.writer.templates.cdr_template_factory import CDR_TemplateFactory


class TemplateFactoryTest(unittest.TestCase):

    def setUp(self):
        self.factory = CDR_TemplateFactory()

    def test_get_cdr_template(self):
        albedo = self.factory.get_cdr_template("ALBEDO")
        self.assertIsNotNone(albedo)

        aot = self.factory.get_cdr_template("AOT")
        self.assertIsNotNone(aot)

        sst = self.factory.get_cdr_template("SST")
        self.assertIsNotNone(sst)

        sst_ensemble = self.factory.get_cdr_template("SST_ENSEMBLE")
        self.assertIsNotNone(sst_ensemble)

        uth = self.factory.get_cdr_template("UTH")
        self.assertIsNotNone(uth)
