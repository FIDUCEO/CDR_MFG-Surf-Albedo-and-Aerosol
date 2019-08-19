import datetime

from data_utility import DataUtility
from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter


def main():
    writer = FCDRWriter()

    # get a template for sensor name in FULL format, supply product height
    # The scan-width is set automatically
    # ---------------------------------------------------------------------
    dataset = writer.createTemplateFull("AVHRR", 128)

    # set some mandatory global attributes (CF standards). Writing will fail if not all of them are filled
    # automatically set: CF version and FIDUCEO license
    # ----------------------------------------------------------------------------------------------------
    dataset.attrs["institution"] = "Brockmann Consult GmbH"
    dataset.attrs["title"] = "FIDUCEO test dataset"
    dataset.attrs["source"] = "arbitray stuff"
    dataset.attrs["history"] = "none"
    dataset.attrs["references"] = "CDR_FCDR sensor reference documentation"
    dataset.attrs["comment"] = "just to show how things are intended to be used"

    # write real data to the variables. All variables initially contain "_FillValue".
    # Not writing to the whole array is completely OK
    # -------------------------------------------------------------------------------
    Time = dataset.variables["Time"]
    Time.data[44] = 0.456
    Time.data[45] = 0.457

    raa = dataset.variables["relative_azimuth_angle"]
    raa.data[3, 0] = 0.567
    raa.data[3, 1] = 0.568

    # ensure not to generate over/underflows
    # --------------------------------------
    DataUtility.check_scaling_ranges(raa)

    # create a standardized file name
    # -------------------------------
    start = datetime.datetime(2006, 8, 23, 14, 24, 52)
    end = datetime.datetime(2006, 8, 23, 15, 25, 53)
    file_name = writer.create_file_name_FCDR_full("AVHRR", "NOAA12", start, end, "01.2")

    # dump it to disk, netcdf4, medium compression
    # overwrite existing file
    # --------------------------------------------
    writer.write(dataset, "D:\\Satellite\\DELETE\\" + file_name, overwrite=True)


if __name__ == "__main__":
    main()
