# $Id$
# Hop --- a framework to analyze solvation dynamics from MD simulations
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Logging class that implements fine grained control over how much
output is produced. It also allows logging to a file.

(As everyone else I am re-inventing the wheel here. Perhaps have a
look at the stdlib class?...)

XXX: use standard logger!
"""

import sys

class Levelcodes:
    __verbosity_levels = {0:'quiet',1:'important',2:'normal',3:'verbose',
                          4:'details',
                          5:'debug_1',6:'debug_2',7:'debug_3',8:'debug_4',
                          9:'debug_5',10:'all'}
    def __init__(self):
        for level,name in self.__verbosity_levels.items():
            self.__dict__[name.upper()] = level            
            
    def get(self,x):
        try:
            if type(x) is int:
                return self.__verbosity_levels[abs(x)]  # abs as <0 means logfile
            else:
                return self.__dict__[x.upper()]
        except KeyError:
            raise ValueError('No level explicitly defined for %r' % x)

    def __getattr__(self,x):
        return self.get(x)

    def __getitem__(self,x):
        return self.get(x)

    def show(self):
        print "Verbosity levels"
        print "%10s %6s" % ('Level','number')
        print "-" * 17
        for level,name in self.__verbosity_levels.items():
            print "%10s %6d" % (name.upper(),level)
        print "-" * 17

levelno = Levelcodes()

class Logger(object):
    def __init__(self,verbosity=levelno.VERBOSE,logfilename='hop.log'):
        super(Logger,self).__init__()
        self.__verbosity = 0                  # need to fill managed variable with number
        self.logfilename = logfilename
        self.logfile = None     # file handle
        # set at the end of init:
        self.verbosity = verbosity       # self.verbosity: managed property 

    def set_verbosity(level=None,logfilename=None):
        """set_verbosity([level],logfilename=<filename>)
        Set the verbosity level
        level < 0 : level <- abs(level) but output is also appended to logfile
        level == 0: minimum
        level == 3: verbose
        level > 3 : debugging"""
        if level is None:
            return self.verbosity
        if logfilename:
            self.logfilename = logfilename
        self.verbosity = level
        return self.verbosity

    def open_log(self,filename):
        """Open the log file <filename> and set this as the log filename."""
        import time
        hrule = "-"*70+'\n'

        self.logfilename = filename
        self.logfile = open(filename,'a')
        self.logfile.write(hrule)
        self.logfile.write("Opened logfile for hop analysis on "+time.strftime('%X %x %Z')+'\n')
        self.logfile.write("verbosity level: %d (%s)\n" % (self.verbosity,levelno[self.verbosity]))

    def close_log(self):
        """Close open logfile"""
        if self.logfile is not None:
            self.logfile.close()

    def verbosity():
        doc = """verbosity level; see hop.logger.levelno.show() for description. level<0 writes to file"""
        def fget(self):
            return self.__verbosity
        def fset(self,level):
            if level >= 0 and self.__verbosity < 0:
                # close the logfile if we don't write to it anymore
                self.close_log()
            self.__verbosity = level   # set the managed attribute!!
            if level < 0:
                logfilename = self.logfilename or "hop.log"  # default log file name
                try:
                    self.open_log(logfilename)
                except:
                    raise IOError("Failed to open logfile.")
        return locals()
    verbosity = property(**verbosity())

    def msg(self,level,m=None):
        """msg(level,[m])
       
        1) Print message string if the level <= verbosity. level describes the
        priority with lower = more important.

        Terminate string with \\n if a newline is desired or \\r to overwrite
        the current line (eg for output progress indication)

        Note that if the global verbosity level is < 0 then the message is also
        written to the logfile.

        2) If called without a string then msg(level) returns True if it would
        print a message, and False otherwise.
        """
        if level <= abs(self.verbosity):
            if m:
                print m, # terminate string with \n if newline desired
                sys.stdout.flush()
                if self.verbosity < 0:
                    # TODO: initial \n mess up formatting but now I also strip initial spaces
                    logmsg = m.strip()   # remove \r (on screen progress) or \n
                    self.logfile.write('[%2d] ' % level + logmsg + '\n')
            return True
        return False

    def __del__(self):
        # make sure that the logfile IS closed properly
        self.close_log()

    def __repr__(self):
        return "Logger(verbosity=%d,logfilename='%s')" % (self.verbosity,self.logfilename)

