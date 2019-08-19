<img alt="FIDUCEO FCDR Tools" align="right" src="http://www.fiduceo.eu/sites/default/files/FIDUCEO-logo.png" />

.

# FCDR Tools

Tools for handling Fiduceo formatted FCDR (Fundamental Climate Data Record) and CDR (Climate Data Record)files.

## Status

[![Build Status](https://travis-ci.org/FIDUCEO/FCDRTools.svg?branch=master)](https://travis-ci.org/FIDUCEO/FCDRTools)
[![codecov.io](https://codecov.io/gh/FIDUCEO/FCDRTools/branch/master/graphs/badge.svg?)](https://codecov.io/gh/FIDUCEO/FCDRTools/branch/master/graphs/badge.svg?)

## Contents

* `reader` - reading utilities
* `writer` - writing utilities
* `test` - test classes

## Dependencies

The FCDR/CDR Tools depend on a number of modules, namely:

* `xarray`
* `netcdf4`
* `numpy`
* `numexpr`
* `dask`
* `gridtools` - (channel "cablab")

## Testing
The FCDR tools are developed using a large set of unit-level and integration tests. To 
check correct functionality please execute in the module root directory:

'python -m unittest discover  -p "*_test.py"'

To execute I/O based (slow) tests, please execute

'python -m unittest discover  -p "*_iotest.py"'
 