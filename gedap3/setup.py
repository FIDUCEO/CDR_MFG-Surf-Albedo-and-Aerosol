"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
import sys
sys.path.append("./gedap/version")
from setuptools import setup, find_packages
from version import __version__ as gedap_version


setup(
    name='gedap',
    version=gedap_version,
    author='Rayference',
    author_email='info@rayference.eu',
    description='',
    long_description='',

    install_requires=[
        'numpy>=1.7',
        'scipy>=0.12',
        'netCDF4>=1.2',
        # Extra requirements for testing applications
        'click>=5.0',
        'matplotlib>=1.5'
    ],
    entry_points={
        'console_scripts': [
            'gedap=gedap.scripts.nodist.main:main_cli'
        ],
    },

    packages=find_packages(),
    package_data={
        'gedap.cpp': ['*.so']
    },
    zip_safe=False,
)
