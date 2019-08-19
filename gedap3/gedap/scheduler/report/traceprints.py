"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
import sys
import traceback
from ..utils.print_tools import colors


class TracePrints(object):
    """This class allows tracking of errant print statements
    (see https://stackoverflow.com/questions/1617494/finding-a-print-statement-in-python).
    Usage: just import this module in the script where you want to activate print stack tracing.
    Warning: the resulting output is extremely verbose."""
    def __init__(self):
        self.stdout = sys.stdout

    def write(self, s):
        self.stdout.write(colors.fg.blue + "Writing %r\n" % s + colors.reset)
        traceback.print_stack(file=self.stdout)


sys.stdout = TracePrints()
