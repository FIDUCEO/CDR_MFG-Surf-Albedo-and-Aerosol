"""
  The copyrights for the GEDAP algorithm and computer codes, remain with 
  Rayference SPRL as an Intellectual Property Right 
 
  Rayference Copyright (c) 
 
"""
import logging.handlers
import os
from ..utils.makedirs import mkdir_p


# cleaning the log file

gedap_env = os.getenv('GEDAP_DIR')

if gedap_env is None:
    raise EnvironmentError("GEDAP_DIR environment variable not set")

#these lines apparently don' do anything since there's no log there
#path_to_log = os.path.join(gedap_env, 'SOURCE/PYTHON_SOURCE/gedap/scripts/nodist/output/gedap.log')
#os.remove(path_to_log) if os.path.exists(path_to_log) else None

LOG_FILENAME = 'output/gedap.log'
mkdir_p(os.path.join(os.getcwd(), 'output/'))


logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
                    datefmt="%Y-%m-%d %H:%M:%S")

logger = logging.getLogger(__name__)


# Add the log message handler to the logger
#handler = logging.handlers.RotatingFileHandler(
#              LOG_FILENAME, maxBytes=20, backupCount=10)

#logger.addHandler(handler)


def cpp_log(level, st_class, st_method, msg):
    """
    Python function that is passed to and used by the C++ code (see the C++ Logger class in folder Common).

    :param level: level of the log message
    :type level: int
    :param st_class: name of the class that sends a message
    :type st_class: string
    :param st_method: name of the method that sends a message
    :type st_method: string
    :param msg: log message
    :type msg: string
    """
    msg = "[" + st_class + "::" + st_method + "] " + msg
    logger.log(int(level), msg)


def print_and_log(msg, log):
    """Prints and logs a message.

    :param msg: String to be printed or logged
    :param log: Logger method to use (e.g. logger.debug)"""
    print(msg)
    log(msg)

