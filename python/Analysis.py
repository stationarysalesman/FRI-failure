"""
@author: Tyler Camp
Copyright Tyler Camp 2016
Class Analysis: An abstract class used as a template for other kinds of
analyses."""

from abc import ABCMeta, abstractmethod, abstractproperty
import errno
import os
from datetime import datetime
import time
import random

class Analysis:
    __metaclass__ = ABCMeta

    @abstractmethod
    def start_analysis(self):
        """Driver for the analysis work flow"""
        return None

    @abstractmethod
    def cleanup(self):
        """Perform general cleanup after work flow completion"""
        return None

    # Boilerplate methods common to most analyses

    def init_temp_dir(self, u_str):
        """Create a directory in tmp using provided string str and a
        pseudorandom string based on the current date and time"""

        rstr = str(datetime.today()) + str(time.time())
        random.seed(rstr)
        val = str(random.random())
        path = "/tmp/" + u_str + val + "/"
        try:
            temp_dir = open(path)
        except IOError as e:
            if e.errno == errno.EACCES:
                print e
                exit(e.errno)
            elif e.errno == errno.ENOENT:
                os.mkdir(path)
            else:
                raise
        # Already exists/ok
        else:
            temp_dir.close()

        return path

    """input_dir: director containing input"""

    _input_dir = None

    def get_input_dir(self):
        return self._input_dir

    def set_inputdirs(self, val):
        self._input_dir = val

    input_dir = abstractproperty(get_input_dir, set_inputdirs)

    """output_dir: directory containing output"""

    _output_dir = None

    def get_output_dir(self):
        return self._outputDir

    def set_output_dir(self, val):
        self._output_dir = val

    output_dir = abstractproperty(get_output_dir, set_output_dir)

    """log file"""

    @abstractproperty
    def log_file(self):
        return None

    @abstractmethod
    def write_log(self, data):
        return None

    """Analysis metadata"""

    @abstractproperty
    def start_time(self):
        return None

    @abstractproperty
    def end_time(self):
        return None
