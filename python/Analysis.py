"""
@author: Tyler Camp
Copyright Tyler Camp 2016
Class Analysis: An abstract class used as a template for other kinds of
analyses."""

from abc import ABCMeta, abstractmethod, abstractproperty

class Analysis:
    __metaclass__ = ABCMeta

    """start_analysis(): driver for the analysis work flow"""
    @abstractmethod
    def start_analysis(self):
        return None

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






