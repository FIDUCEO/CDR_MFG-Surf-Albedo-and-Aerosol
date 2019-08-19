from fiduceo.common.version import __version__

class WriterUtils:

    @staticmethod
    def add_standard_global_attributes(dataset):
        dataset.attrs["Conventions"] = "CF-1.6"
        dataset.attrs[
            "licence"] = "This dataset is released for use under CC-BY licence (https://creativecommons.org/licenses/by/4.0/) and was developed in the EC " \
                         "FIDUCEO project \"Fidelity and Uncertainty in Climate Data Records from Earth " \
                         "Observations\". Grant Agreement: 638822."
        dataset.attrs["writer_version"] = __version__

        # The following dictionary entries have to be supplied by the data generators
        dataset.attrs["institution"] = None
        dataset.attrs["title"] = None
        dataset.attrs["source"] = None
        dataset.attrs["history"] = None
        dataset.attrs["references"] = None
        dataset.attrs["comment"] = None

    @staticmethod
    def add_cdr_global_attributes(dataset):
        # The following dictionary entries have to be supplied by the data generators
        dataset.attrs["source"] = None
        dataset.attrs["auxiliary_data"] = None
        dataset.attrs["configuration"] = None
        dataset.attrs["time_coverage_start"] = None
        dataset.attrs["time_coverage_end"] = None
        dataset.attrs["time_coverage_duration"] = None
        dataset.attrs["time_coverage_resolution"] = None

    @staticmethod
    def add_gridded_global_attributes(dataset):
        # The following dictionary entries have to be supplied by the data generators
        dataset.attrs["geospatial_lat_units"] = None
        dataset.attrs["geospatial_lon_units"] = None
        dataset.attrs["geospatial_lat_resolution"] = None
        dataset.attrs["geospatial_lon_resolution"] = None