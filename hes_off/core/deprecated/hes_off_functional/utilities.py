## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
##                         ___   ___  ________    ________         ______    ________ ________                        ##
##                        |  |  |  | |   ____|   /       |        /  __  \  |   ____||   ____|                        ##
##                        |  |__|  | |  |__     |   (----` ______|  |  |  | |  |__   |  |__                           ##
##                        |   __   | |   __|     \   \    |______|  |  |  | |   __|  |   __|                          ##
##                        |  |  |  | |  |____.----)   |          |  `--'  | |  |     |  |                             ##
##                        |__|  |__| |_______|_______/            \______/  |__|     |__|                             ##
##                                                                                                                    ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##

# Import packages
import re
import sys
import pdb
import copy
import time
import functools
import numpy as np

## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self, name=None):
        self._start_time = None
        self.name = name

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        if self.name:
            print("Elapsed time for", self.name, "is: {:0.5e} seconds".format(elapsed_time))
        else:
            print("Elapsed time is: {:0.5e} seconds".format(elapsed_time))

    def __enter__(self):
        """Start a new timer as a context manager"""
        self.start()
        return self

    def __exit__(self, *exc_info):
        """Stop the context manager timer"""
        self.stop()

    def __call__(self, func):
        """Support using Timer as a decorator"""
        @functools.wraps(func)
        def wrapper_timer(*args, **kwargs):
            with self:
                return func(*args, **kwargs)
        return wrapper_timer


def read_configuration_file(filename):
    """ Read the configuration file that stores the input data """
    try:
        IN = {}
        with open(filename, 'r') as reader:
            for line in reader:
                line  = re.sub(" ", "", line)                               # Eliminate white spaces
                line  = re.sub("[ \t\n\"()\[\]{}]", "", line)               # Eliminate " [ ] ( ) { } \n characters
                words = re.split("[=,;#%]", line)                           # Split over = , ; # % characters
                if not any([words[0] in symbol for symbol in ['%','#']]):   # Ignore lines starting by % or #
                    words = list(filter(None, words))                       # Variables with no value are empty lists
                    for i in range(len(words)):
                        try:
                            words[i] = float(words[i])                      # Try to convert to float when possible
                        except:
                            words[i] = words[i]                             # Handle as a string otherwise
                    if len(words[1:]) == 1:
                        IN[words[0]] = words[1]                             # Define variable as a float or string
                    else:
                        IN[words[0]] = words[1:]                            # Define variable as a list

            # Add configuration file path to the dictionary
            IN['CONFIGURATION_PATH'] = filename
    except:
        raise Exception('\n\n\n The configuration file could not be read. Exiting the program.')

    # Convert units from MW to W
    for key, value in IN.items():
        if "POWER" in key or "HEAT" in key:
            if isinstance(value, list):
                IN[key] = [item*1e6 for item in value]
            else:
                IN[key] = value*1e6

    return IN


def write_configuration_file(filename, IN):
    """ Write a configuration file from dictionary """
    with open(filename, 'w') as writer:
        for key, value in IN.items():

            if "POWER" in key or "HEAT" in key:
                if isinstance(value, list):
                    value = [item / 1e6 for item in value]
                else:
                    value = value / 1e6

            if not isinstance(value, list):
                if not isfloat(value):
                    value = '"' + str(value) + '"'
            writer.write("%s = %s\n" % (key, value))


def print_dictionary(IN):
    """ Print a configuration dictionary with a readable format """
    for key, value in IN.items():
        print("{:<24}{}".format(key, value))


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def print_progress (iteration, total, prefix ='', suffix ='', decimals = 2, barLength = 100):

    """ Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)

    """

    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    if percents>100.0: percents=100.0
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s [%s] %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")


def print_banner():
    """ Print the HES-OFF banner"""
    print("## ---------------------------------------------------------------------------------- ##")
    print("## ---------------------------------------------------------------------------------- ##")
    print("##         ___   ___  ________    ________         ______    ________ ________        ##")
    print("##        |  |  |  | |   ____|   /       |        /  __  \  |   ____||   ____|        ##")
    print("##        |  |__|  | |  |__     |   (----` ______|  |  |  | |  |__   |  |__           ##")
    print("##        |   __   | |   __|     \   \    |______|  |  |  | |   __|  |   __|          ##")
    print("##        |  |  |  | |  |____.----)   |          |  `--'  | |  |     |  |             ##")
    print("##        |__|  |__| |_______|_______/            \______/  |__|     |__|             ##")
    print("##                                                                                    ##")
    print("## ---------------------------------------------------------------------------------- ##")
    print("## ---------------------------------------------------------------------------------- ##")
    print('\n')

