# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 20:34:52 2015

@author: tyler

This file is part of fri-failure-analysis.

    fri-failure-analysis is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fri-failure-analysis is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fri-failure-analysis.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import datetime

class jobmanager:
    def __init__(self, input_dirs, output_dirs, logfile_name, process_module):
        self.master_output_dir = "./output/"
        self.input_dirs = input_dirs
        self.output_dirs = output_dirs
        self.logfile_name = self.master_output_dir+logfile_name
        self.process_module = __import__(process_module)
        self.__initialize_input__()
        self.__initialize_output__()
        self.__initialize_logfile__()
        return
    @property
    def input_dirs(self):
        """Return a dictionar of input directories."""
        return self.input_dirs
    @property
    def output_subdir_list(self):
        """Return a list of output subdirectories."""
        return self.output_subdir_list
    @property
    def logfile_name(self):
        """Return the name of the logfile."""
        return self.logfile_name
    @property
    def process_module(self):
        """Return the name of the process module."""
        return self.process_module
        
    def __initialize_input__(self):
        """Initialize master input directory and check for access to subdirectories."""
        for key in self.input_dirs.keys():
            if not(os.access(self.input_dirs[key], os.F_OK)):
                print "Error: directory", key, "inaccessible."
                return
        return
    def __initialize_output__(self):
        """Initialize all output directories."""
        if not(os.access(self.master_output_dir, os.F_OK)):
            os.mkdir(self.master_output_dir)
        for key in self.output_dirs.keys():
            if not(os.access(self.master_output_dir+self.output_dirs[key], os.F_OK)):
                os.mkdir(self.master_output_dir+self.output_dirs[key])
            self.output_dirs[key] = self.master_output_dir + self.output_dirs[key]
        return
    def build_input_path(self, path):
        """Return master input directory prepended to given path."""
        for p in self.input_subdir_list:
            if (p == path):
                return self.master_input_dir + p
    
        
    def __initialize_logfile__(self):
        """Create the log file."""
        with open(self.logfile_name, "a") as log:
            log.write("Logfile created " + str(datetime.datetime.today()))
