# Build and installation instructions for GEDAP

_This document is written with GiHub flavored Markdown (GFM) and can be rendered using the [Atom](https://atom.io/) text editor or the GitHub rendering API (e.g. using the Python [Grip](https://github.com/joeyespo/grip) package)._

## System requirements

gedap3 requires the same software packages to be installed on the target system as the standard GEDAP version, besides for Python which needs to be updated to:

Software package                   | Required version
-----------------------------------|-----------------
Python                             | 3.6

as GEDAP, gedap3 is also shipped with the pybind11 Python binding library. Additionally, gedap3 is also shipped with the FCDRTools library that is used within the FIDUCEO project to generate CDRs.

The Python layer of GEDAP requires the following Python packages be installed on the target system:

* FCDRTools library from Brockmann
* setuptools
* numpy
* scipy
* netCDF4
* matplotlib (testing applications)
* click (testing applications)

FCDRTools requires the following additional packages:
* coverage 
* dask 
* gridtools 
* numexpr 
* xarray

### Installation procedure on Ubuntu

gedap3 is developed with Anaconda. The provided environment.yml file is a merge of the environment containing gedap dependencies and the environment containing FCDRTools dependencies. The NetCDF-CXX4 is now also available in Anaconda.



## Environment variables

The installation and execution of GEDAP makes use of the environment variable `$GEDAP_DIR`, which must point to the root of the GEDAP project tree. This can be automatically set up on shell login by adding the following line to the `.profile` file:

```
$ export GEDAP_DIR=<absolute_path_to_gedap_project>
```

## Building GEDAP

1. Navigate to the GEDAP project root directory:

    ```
    $ cd $GEDAP_DIR
    ```
2. Create a build directory and navigate to it:

    ```
    $ mkdir build
    $ cd build
    ```
3. Execute the CMake script to generate a makefile (**see below for detailed build configuration instructions**):

   ```
   $ cmake ..
   ```

4. Build the C++ layer of GEDAP:

    ```
    $ make
    ```

   If this step of the procedure fails, it is likely that the dependencies were inappropriately installed or registered to CMake. See the next section.
5. Navigate to the root directory of GEDAP:

    ```
    $ cd $GEDAP_DIR
    ```
6. Install GEDAP to the Python environment:

    ```
    $ python setup.py install
    ```

7. Install FCDRTools to the Python environment:

Navigate to the root directory of the FCDRTools 

```
$ cd $GEDAP_DIR/FCDRTools-2_0_0_RELEASE
```
and install the package

```
$ python setup.py install
```


   For an installation in developer mode, see Troubleshooting.

## Fine-tuning the build

> The build of gedap3 has run smoothly on a Ubuntu machine without having to set any additional cmake flag reported below. However, if this is not the case at Eumetsat, don't hesitate emailing us or try debugging using the instructions below.


GEDAP relies on a series of dependencies which have to be carefully pointed to during the build process. The CMake build script automatically looks up compilers and libraries, but it can be mistaken on complex systems with multiple compilers or libraries installed. In such cases, the system can be hinted for searches.

In the following, variables which can point CMake to the right direction for its search will be defined. To configure the CMake build with variable `VAR` having the value `SOME_VALUE`, CMake should be called using the `-D` flag:
```
$ cmake -DVAR=SOME_VALUE ..
```
Variables can also be configured by editing `CMakeCache.txt` or through CMake GUIs `ccmake` and `cmake-gui`.

### Pointing to the right compiler

Compilers can be defined through the following CMake variables:

Language | Variable
---------|-----------------------
Fortran  | CMAKE_Fortran_COMPILER
C        | CMAKE_C_COMPILER
C++      | CMAKE_CXX_COMPILER

### Pointing to the right Python executable

The Python interpreter can be hinted for as well through the `PYTHON_EXECUTABLE` CMake variable. If this variable is not defined, the `GEDAP_PYTHON_EXECUTABLE` environment variable will be looked up. If it is not defined, the `CONDA_PREFIX` environment variable will be looked up.

1                   | 2                              | 3                   | 4
--------------------|--------------------------------|---------------------|---------------------------
`PYTHON_EXECUTABLE` | `ENV{GEDAP_PYTHON_EXECUTABLE}` | `ENV{CONDA_PREFIX}` | Default (system-dependent)

### Pointing to the right libraries

If GEDAP is srunning in a Conda environment, the option `ENABLE_AUTO_CONDA_HINTING` can be activated to have CMake look for libraries and headers in the current Conda environment (i.e., the directory pointed by environment variable `$CONDA_PREFIX`).

Library     | 1                            | 2                             | 3
------------|------------------------------|-------------------------------|---------------------------
gfortran    | `LIBGFORTRAN_LIBRARIES_HINT` | Default (using `gfortran -v`) |
HDF5        | `HDF5_ROOT`                  | `ENV{CONDA_PREFIX}`           | Default (system-dependent)
NetCDF      | `NETCDF_DIR`                 | `ENV{CONDA_PREFIX}`           | Default (system-dependent)
NetCDF-CXX4 | `NetCDF_CXX4_DIR`            | `ENV{CONDA_PREFIX}`           | Default (system-dependent)
Eigen       | `Eigen_DIR`                  | `ENV{CONDA_PREFIX}`           | Default (system-dependent)
CISAR       | `CISAR_DIR`                  | `ENV{CISAR_DIR}`              | None, error

### Build options

The following build options are available:

Option                      | Type | Default | Description
----------------------------|------|---------|----------------------------------------------------------------
`ENABLE_AUTO_CONDA_HINTING` | BOOL | ON      | Use `ENV{CONDA_PREFIX}` automatically for library search hints
`ENABLE_INSTANTANEOUS`      | BOOL | ON      | INSTANTANEOUS build (uses the INSTANTANEOUS CISAR library file)
`ENABLE_OPENMP`             | BOOL | ON      | Link OpenMP library \*
`ENABLE_SIMPLEPREC`         | BOOL | OFF     | CISAR was built with single precision floating point numbers \*
`ENABLE_TESTING`            | BOOL | OFF     | Compile test programs

**\* Must be consistent with CISAR's build settings** (if not, GEDAP will most probably not work)

## Checking the installation

1. Check that the Python package is correctly installed by issuing the following command to the terminal:

    ```
    $ which gedap
    ```

   If this returns a blank string despite the installation procedure was successful, it is likely that the Python binary directory is not included in the `$PATH`. See Troubleshooting for more information.

2. Run the GEDAP script to get the help text by issuing the following command to the terminal:

    ```
    $ gedap -h
    ```

   The program should display the help text then exit without error.

## Troubleshooting

* __Missing NetCDF header file.__ On many distributions, the NetCDF-CXX4 package do not have the development version of the NetCDF-C package as a dependency. Make sure that the development versions of HDF5, NetCDF and NetCDF-CXX4 are installed.

* __Python binary directory missing from `$PATH`.__ On some systems, it can happen that the Python binary directory is not included in the `$PATH` variable. If `gedap -h` returns a blank string, it means that the GEDAP executable cannot be found in the directories contained in `$PATH`. To solve this:

    1. Run the Python setup script again in verbose mode:

        ```
        $ python setup.py --verbose install
        ```
    2. Find out where the GEDAP script was installed: in the output of the install script, look for a line which looks like

        ```
        Installing gedap script to <some_path>
        ```
    3. Append `<some_path>` to `$PATH` (preferably in your `.profile`).

* __Installation without root privileges.__ Installing Python packages does not necessarily require root privileges. Even when not using virtual environments, it is possible to install GEDAP and its dependencies as a regular user by simply adding the `--user` option to the Python setup command:

    ```
    $ python setup.py install --user
    ```

* __Installation in developer mode.__ GEDAP can be installed in developer mode (see the help text of `setup.py`). In that case, the package is installed to the Python environment through references, allowing the user to see their changes immediately instead of having to run again `python setup.py install`. To do so, step #6 of the installation procedure has to be replaced with:

    ```
    $ python setup.py develop
    ```


## Installing GEDAP on Linux (Ubuntu 16.04)

* Errors compiling and installing the netCDF-cxx4 library: it has been impossible to follow the above descripted procedure to install the netCDF-cxx4 library within the conda env on Linux (several ``undefined references" errors). The solution found has been to use the conda environment only for the python related things, while the libraries netcdf, hdf5 and netcdf-cxx4 are the system wide ones. It is important to verify that all the cmake variables point to the right libraries to avoid conflicts.

* Make sure that the compilers used on Linux are the native ones and not the ones shipped with Conda: a known bug is indeed present in the fortran compiler and doesn't allow GEDAP to run correctly (get_unit(): Bad internal unit KIND error).
