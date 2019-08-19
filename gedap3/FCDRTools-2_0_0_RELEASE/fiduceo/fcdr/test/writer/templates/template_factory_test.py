import unittest

from fiduceo.fcdr.writer.templates.avhrr_flag_mapper import AVHRR_FlagMapper
from fiduceo.fcdr.writer.templates.default_flag_mapper import DefaultFlagMapper
from fiduceo.fcdr.writer.templates.hirs_flag_mapper import HIRS_FlagMapper
from fiduceo.fcdr.writer.templates.mviri_flag_mapper import MVIRI_FlagMapper
from fiduceo.fcdr.writer.templates.template_factory import TemplateFactory


class TemplateFactoryTest(unittest.TestCase):

    def setUp(self):
        self.factory = TemplateFactory()

    def test_get_sensor_template(self):
        amsub = self.factory.get_sensor_template("AMSUB")
        self.assertIsNotNone(amsub)

        mhs = self.factory.get_sensor_template("MHS")
        self.assertIsNotNone(mhs)

        ssmts = self.factory.get_sensor_template("SSMT2")
        self.assertIsNotNone(ssmts)

        avhrr = self.factory.get_sensor_template("AVHRR")
        self.assertIsNotNone(avhrr)

        hirs_2 = self.factory.get_sensor_template("HIRS2")
        self.assertIsNotNone(hirs_2)

        hirs_3 = self.factory.get_sensor_template("HIRS3")
        self.assertIsNotNone(hirs_3)

        hirs_4 = self.factory.get_sensor_template("HIRS4")
        self.assertIsNotNone(hirs_4)

        mviri = self.factory.get_sensor_template("MVIRI")
        self.assertIsNotNone(mviri)

        mviri_static = self.factory.get_sensor_template("MVIRI_STATIC")
        self.assertIsNotNone(mviri_static)

    def test_get_flag_mapper(self):
        amsub = self.factory.get_flag_mapper("AMSUB")
        self.assertIsNotNone(amsub)
        self.assertIsInstance(amsub, DefaultFlagMapper)

        mhs = self.factory.get_flag_mapper("MHS")
        self.assertIsNotNone(mhs)
        self.assertIsInstance(mhs, DefaultFlagMapper)

        ssmt2 = self.factory.get_flag_mapper("SSMT2")
        self.assertIsNotNone(ssmt2)
        self.assertIsInstance(ssmt2, DefaultFlagMapper)

        avhrr = self.factory.get_flag_mapper("AVHRR")
        self.assertIsNotNone(avhrr)
        self.assertIsInstance(avhrr, AVHRR_FlagMapper)

        hirs2 = self.factory.get_flag_mapper("HIRS2")
        self.assertIsNotNone(hirs2)
        self.assertIsInstance(hirs2, HIRS_FlagMapper)

        hirs3 = self.factory.get_flag_mapper("HIRS3")
        self.assertIsNotNone(hirs3)
        self.assertIsInstance(hirs3, HIRS_FlagMapper)

        hirs4 = self.factory.get_flag_mapper("HIRS4")
        self.assertIsNotNone(hirs4)
        self.assertIsInstance(hirs4, HIRS_FlagMapper)

        mviri = self.factory.get_flag_mapper("MVIRI")
        self.assertIsNotNone(mviri)
        self.assertIsInstance(mviri, MVIRI_FlagMapper)
