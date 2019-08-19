import unittest

from fiduceo.fcdr.reader.fcdr_reader import FCDRReader
from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter


class FCDRReaderMviriVirtualVariablesTest(unittest.TestCase):

    def setUp(self):
        self.fcdr_reader = FCDRReader()

        fcdr_writer = FCDRWriter()
        self.dataset = fcdr_writer.createTemplateFull("MVIRI", 5000)

    def testCalculate_sensitivitiy_solar_irradiance(self):
        self.dataset["distance_sun_earth"].data = 1.0166579484939575
        self.dataset["count_vis"].data[:, :] = 63
        self.dataset["mean_count_space_vis"] = 4.961684375
        self.dataset["a0_vis"].data = 0.9800095200636486
        self.dataset["a1_vis"].data = 0.01179638707394702
        self.dataset["a2_vis"].data = 0.001179638707394702
        self.dataset["years_since_launch"].data = 8.830136986301369
        self.dataset["solar_zenith_angle"].data[:, :] = 49.3436
        self.dataset["solar_irradiance_vis"].data = 688.144781045

        self.fcdr_reader._load_virtual_variable(self.dataset, "sensitivity_solar_irradiance_vis")

        virtual_variable = self.dataset["sensitivity_solar_irradiance_vis"]
        self.assertIsNotNone(virtual_variable)
        self.assertEqual((5000, 5000), virtual_variable.shape)

        self.assertAlmostEqual(0.00071843801878964449, virtual_variable.data[0, 0])

    def testCalculate_sensitivitiy_count_vis(self):
        self.dataset["distance_sun_earth"].data = 1.0166579484939575
        self.dataset["a0_vis"].data = 0.9800095200636486
        self.dataset["a1_vis"].data = 0.01179638707394702
        self.dataset["a2_vis"].data = 0.02179638707394702
        self.dataset["years_since_launch"].data = 8.830136986301369
        self.dataset["solar_zenith_angle"].data[:, :] = 35.4024
        self.dataset["solar_irradiance_vis"].data = 688.144781045

        self.fcdr_reader._load_virtual_variable(self.dataset, "sensitivity_count_vis")

        virtual_variable = self.dataset["sensitivity_count_vis"]
        self.assertIsNotNone(virtual_variable)
        self.assertEqual((500, 500), virtual_variable.shape)

        self.assertAlmostEqual(0.016114772327457307, virtual_variable.data[1, 1])

    def testCalculate_sensitivitiy_count_space(self):
        self.dataset["distance_sun_earth"].data = 1.0166579484939575
        self.dataset["a0_vis"].data = 0.9800095200636486
        self.dataset["a1_vis"].data = 0.01179638707394702
        self.dataset["a2_vis"].data = 0.01379638707394702
        self.dataset["years_since_launch"].data = 8.830136986301369
        self.dataset["solar_zenith_angle"].data[:, :] = 26.9212
        self.dataset["solar_irradiance_vis"].data = 688.144781045

        self.fcdr_reader._load_virtual_variable(self.dataset, "sensitivity_count_space")

        virtual_variable = self.dataset["sensitivity_count_space"]
        self.assertIsNotNone(virtual_variable)
        self.assertEqual((500, 500), virtual_variable.shape)

        self.assertAlmostEqual(-0.011430569161934807, virtual_variable.data[2, 2])

    def testCalculate_sensitivitiy_a0_vis(self):
        self.dataset["distance_sun_earth"].data = 1.0166579484939575
        self.dataset["count_vis"].data[:, :] = 14
        self.dataset["mean_count_space_vis"] = 4.961684375
        self.dataset["solar_zenith_angle"].data[:, :] = 18.2038
        self.dataset["solar_irradiance_vis"].data = 688.144781045

        self.fcdr_reader._load_virtual_variable(self.dataset, "sensitivity_a0_vis")

        virtual_variable = self.dataset["sensitivity_a0_vis"]
        self.assertIsNotNone(virtual_variable)
        self.assertEqual((5000, 5000), virtual_variable.shape)

        self.assertAlmostEqual(0.044895821172659549, virtual_variable.data[3, 3])

    def testCalculate_sensitivitiy_a1_vis(self):
        self.dataset["distance_sun_earth"].data = 1.0166579484939575
        self.dataset["years_since_launch"].data = 8.830136986301369
        self.dataset["count_vis"].data[:, :] = 24
        self.dataset["mean_count_space_vis"] = 4.961684375
        self.dataset["solar_zenith_angle"].data[:, :] = 22.1907
        self.dataset["solar_irradiance_vis"].data = 688.144781045

        self.fcdr_reader._load_virtual_variable(self.dataset, "sensitivity_a1_vis")

        virtual_variable = self.dataset["sensitivity_a1_vis"]
        self.assertIsNotNone(virtual_variable)
        self.assertEqual((5000, 5000), virtual_variable.shape)

        self.assertAlmostEqual(0.85671563074783341, virtual_variable.data[4, 4])

    def testCalculate_sensitivitiy_a2_vis(self):
        self.dataset["distance_sun_earth"].data = 1.0166579484939575
        self.dataset["years_since_launch"].data = 8.830136986301369
        self.dataset["count_vis"].data[:, :] = 24
        self.dataset["mean_count_space_vis"] = 4.961684375
        self.dataset["solar_zenith_angle"].data[:, :] = 22.1907
        self.dataset["solar_irradiance_vis"].data = 688.144781045

        self.fcdr_reader._load_virtual_variable(self.dataset, "sensitivity_a2_vis")

        virtual_variable = self.dataset["sensitivity_a2_vis"]
        self.assertIsNotNone(virtual_variable)
        self.assertEqual((5000, 5000), virtual_variable.shape)

        self.assertAlmostEqual(7.5649163778089505, virtual_variable.data[4, 4])
