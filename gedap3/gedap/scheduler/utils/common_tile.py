"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
"""
module common_tile
------------------

:author: Alix Damman
:date: Created on 06 Sep 2016

This module contains useful code to be use be other modules from this package.
"""

rpv_var_names = ['rho_0', 'k', 'theta', 'rho_c'] 
''' names of the 4 RPV parameters '''


class Aerosol_Class(object):
    """
    This class represents an aerosol class (see CISAR)

    :param name: name of the aerosol class
    :type name: string
    :param type: type of the aerosol class (Fine, Coarse)
    :type type: string
    """
    
    def __init__(self, name="", type=""):
        self.name = name
        self.type = type
        
    def is_fine_mode(self):
        return self.type == "Fine"
    
    def is_coarse_mode(self):
        return self.type == "Coarse"
        
    def __repr__(self):
        return "name: " + self.name + " - type: " + self.type
