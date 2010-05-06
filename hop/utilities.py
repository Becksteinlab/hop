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

"""Random mix of convenience functions that don't fit anywhere
else. For messages I should probably use python's logger module but
this is working so far (even though it's pretty crappy)."""

import sys
import os.path
import cPickle
import warnings
import hop

# unbound methods filename_function(), to be used in other
# classes; the plan is to make all classes that require them
# subclasses of hop.utilities.Saveable and bind them to this super
# class. (load() and save() are already tentatively XXXed out)

# used in many classes for filename handling (not all are Saveable yet)
#   filename = hop.utilities.filename_function
# (adds the _filename attribute to the class self!)
# NOTE: filename_function is not developed anymore and deprecated
#       Try to derive classes from Saveable.
def filename_function(self,filename=None,ext=None,set_default=False,use_my_ext=False):
    """Supply a file name for the object.

    fn = filename()             ---> <default_filename>
    fn = filename('name.ext')   ---> 'name'
    fn = filename(ext='pickle') ---> <default_filename>'.pickle'
    fn = filename('name.inp','pdf') --> 'name.pdf'
    fn = filename('foo.pdf',ext='png',use_my_ext=True) --> 'foo.pdf'

    The returned filename is stripped of the extension (use_my_ext=False) and
    if provided, another extension is appended. Chooses a default if no
    filename is given.  Raises a ValueError exception if no default file name
    is known.

    If set_default=True then the default filename is also set.

    use_my_ext=True lets the suffix of a provided filename take priority over a
    default ext(tension).
    """
    if filename is None:
        if not hasattr(self,'_filename'):
            self._filename = None        # add attribute to class 
        if self._filename:
            filename = self._filename
        else:
            raise ValueError("A file name is required because no default file name was defined.")
        my_ext = None
    else:
        filename, my_ext = os.path.splitext(filename)
        if set_default:                  # replaces existing default file name
            self._filename = filename
    if my_ext and use_my_ext:  
        ext = my_ext
    if ext is not None:
        if ext.startswith('.'):
            ext = ext[1:]  # strip a dot to avoid annoying mistakes
        filename = filename + '.' + ext
    return filename

def fileextension(filename,default=None):
    """Return the file extension without the leading dot or the default."""
    ext = os.path.splitext(filename)[1]
    if len(ext) > 1:
        if ext.startswith('.'):
            ext = ext[1:]
        return ext
    else:
        return default

# load() and save() are unbound plugin methods for other classes
# (Once the code has been cleaned up they will ONLY exist as methods of Saveable)
def XXXsave(self,filename=None):
    """Save as a pickled python object."""
    if not self._saved_attributes:
        warnings.warn("No data save: the object declared empty '_saved_attributes'.")
        return
    # TODO: need to properly design this whole pickling stuff
    data = {}
    if self._saved_attributes == 'all':
        # HACK: manually filter out some attributes such as type('method-wrapper') 
        #       objects that cannot be pickled
        saved_attributes = [x for x in self.__dict__.keys() if x not in self._excluded_attributes]
    else:
        saved_attributes = self._saved_attributes        
    for attr in saved_attributes:
        try:
            data[attr] = self.__dict__[attr]
        except KeyError:
            warnings.warn("Attribute '"+attr+"' has not been computed and will not be saved.",
                          category=hop.MissingDataWarning)
    fh = open(self.filename(filename,'pickle',set_default=True),'wb')  # 2.5: with open(..) as fh:
    try: 
        cPickle.dump(data,fh,cPickle.HIGHEST_PROTOCOL)
    finally:
        fh.close()
    del data

def XXXload(self,filename=None,merge=False):         
    """Reinstantiate class from a pickled file (produced with save())."""
    if not self._saved_attributes:
        warnings.warn("No data loaded: the object declared empty '_saved_attributes'.")
        return
    fh = open(self.filename(filename,'pickle',set_default=True),'rb')  # 2.5: with open(..) as fh:
    try: 
        data = cPickle.load(fh)
    finally:
        fh.close()
    if self._saved_attributes == 'all':
        saved_attributes = data.keys()
    else:
        saved_attributes = self._saved_attributes
    # restore attributes from the temporary instance (works but not elegant...)
    for attr in saved_attributes:
        try:
            if merge and attr in _merge_attributes: # only works for type(attr)==dict
                self.__dict__[attr].update(data[attr])
            else:
                self.__dict__[attr] = data[attr]
        except KeyError:
            warnings.warn("Expected attribute '"+attr+"' was not found in saved file '"+filename+"'.", 
                          category=hop.MissingDataWarning)
    del data


class Saveable(object):
    """Baseclass that supports save()ing and load()ing.
    
    Override the class variables 

      _saved_attributes = [] # list attributes to be pickled
      _merge_attributes = [] # list dicts to be UPDATED from the pickled file with load(merge=True)
      _excluded_attributes = [] # list attributes that should never be pickled
      
    Note: 

      _saved_attributes = 'all' # pickles ALL attributes, equivalent to self.__dict__.keys()
                                # (use _excluded_attributes with 'all'!)

     Use _excluded_attributes to filter out some attributes such as
     type('method-wrapper') objects that cannot be pickled (e.g. when using properties).
    """

    _saved_attributes = []
    _merge_attributes = []
    _excluded_attributes = []

    # TODO: could I use __new__ to load the pickled file?
    def __init__(self,*args,**kwargs):
        # TODO: should initialize _XXX_attributes[] via __init__() and use super(cls,Saveable).__init__()
        #       in subclasses
        kwargs.setdefault('filename',None)
        # Multiple Inheritance is probably NOT going to work...                  
        super(Saveable,self).__init__()        # XXX: ... shouldn't this take *args,**kwargs ?? OB-2009-06-10
        if kwargs['filename'] is not None:
            self.load(kwargs['filename'],'pickle')   # sets _saved_attributes in __dict__
        else:
            pass # do the initialization from the data

    # basic pickle protocol
    def __getstate__(self):
        if not self._saved_attributes:
            warnings.warn("No data saved: the class declared empty '_saved_attributes'.")
            return False
        data = {}
        if self._saved_attributes == 'all':
            # HACK: filter out some attributes such as type('method-wrapper') objects that
            #       cannot be pickled
            saved_attributes = [x for x in self.__dict__.keys() if x not in self._excluded_attributes]
        else:
            saved_attributes = self._saved_attributes        
        for attr in saved_attributes:
            try:
                data[attr] = self.__dict__[attr]
            except KeyError:
                warnings.warn("Attribute '"+attr+"' has not been computed and will not be saved.",
                              category=hop.MissingDataWarning)
        return data

    def __setstate__(self,data):
        # Simple:
        self.__dict__.update(data)  # could check for _saved_attributes

    def filename(self,filename=None,ext=None,set_default=False,use_my_ext=False):
        """Supply a file name for the object.

        fn = filename()             ---> <default_filename>
        fn = filename('name.ext')   ---> 'name'
        fn = filename(ext='pickle') ---> <default_filename>'.pickle'
        fn = filename('name.inp','pdf') --> 'name.pdf'
        fn = filename('foo.pdf',ext='png',use_my_ext=True) --> 'foo.pdf'

        The returned filename is stripped of the extension (use_my_ext=False) and
        if provided, another extension is appended. Chooses a default if no
        filename is given.  Raises a ValueError exception if no default file name
        is known.

        If set_default=True then the default filename is also set.

        use_my_ext=True lets the suffix of a provided filename take priority over a
        default ext(tension).
        """
        if filename is None:
            if not hasattr(self,'_filename'):
                self._filename = None        # add attribute to class 
            if self._filename:
                filename = self._filename
            else:
                raise ValueError("A file name is required because no default file name was defined.")
            my_ext = None
        else:
            filename, my_ext = os.path.splitext(filename)
            if set_default:                  # replaces existing default file name
                self._filename = filename
        if my_ext and use_my_ext:  
            ext = my_ext
        if ext is not None:
            if ext.startswith('.'):
                ext = ext[1:]  # strip a dot to avoid annoying mistakes
            filename = filename + '.' + ext
        return filename

    def save(self,filename=None):
        """Save class to a pickled file."""
        fh = open(self.filename(filename,'pickle',set_default=True),'wb') # 2.5: with open(..) as fh:
        try:
            cPickle.dump(self,fh,cPickle.HIGHEST_PROTOCOL)
        finally:
            fh.close()

    def load(self,filename=None,merge=False):         
        """Reinstantiate class from a pickled file (produced with save())."""
        fh = open(self.filename(filename,'pickle',set_default=True),'rb')
        try:
            tmp = cPickle.load(fh)
        finally:
            fh.close()
        try:
            data = tmp.__dict__
        except AttributeError:
            warnings.warn("Loading an old-style save file.",category=DeprecationWarning)
            data = tmp  # just hope we got all attributes...
            # backwards-compatibility hacks:
            try:    
                data['P'] = data['parameters']
            except KeyError:
                pass
        del tmp
        # restore attributes from the temporary instance
        if not self._saved_attributes:
            warnings.warn("No data loaded: the object declared empty '_saved_attributes'.",
                          category=hop.MissingDataWarning)
            return False
        if self._saved_attributes == 'all':
            saved_attributes = [x for x in data.keys() if x not in self._excluded_attributes]
        else:
            saved_attributes = self._saved_attributes
        for attr in saved_attributes:
            try:
                if merge and attr in self._merge_attributes: # only works for type(attr)==dict
                    self.__dict__[attr].update(data[attr])
                else:
                    self.__dict__[attr] = data[attr]
            except KeyError:
                warnings.warn("Expected attribute '"+attr+"' was not found in saved file '"+filename+"'.", 
                              category=hop.MissingDataWarning)
        del data
        return True

def easy_load(names,baseclass,keymethod):
    """Instantiate a class either from an existing instance or a pickled file.
    
    instance_list = easy_load(names,baseclass,keymethod)
    
    >>> x = easy_load(<filename>,Xclass,'my_method_name')
    >>> [x1,x2,...] = easy_load([<filename1>, <fn2>,...], Xclass,'my_method_name')
    >>> [x1,x2,...] = easy_load([x1, x2, ..], Xclass,'my_method_name')
    
    If the argument does not implement the keymethod then try loading
    from a file. 

    API:

    For this to work, the baseclass (eg Saveable) must be able to instantiate 
    itself using

    x = baseclass(filename=name)

    If a single name is given, a singlet is returned, otherwise a list of instances.

    (The docs are longer than the code...)
    """
    def load(name):
        if hasattr(name,keymethod):
            x = name
        else:
            x = baseclass(filename=name)
        return x
    if not iterable(names):
        return load(names)
    return [load(x) for x in names]

#------------------------------------------------------------
# status message functions
#------------------------------------------------------------
# Set the global verbosity = 0 to define the default verbosity level.
# I haven't found a way to make this a class -- unless I let every
# other class inherit from the base class (eg 'interactive') and make
# verbosity a class-level variable (or whatever it is called) .
#
# (NOTE: this will be removed once we use logger)

verbosity = 3
logfile = None
LOGFILENAME = None

def set_verbosity(level=None,logfilename=None):
    """set_verbosity([level],logfilename=<filename>)
    Set the verbosity level
    level < 0 : level <- abs(level) but output is also appended to logfile
    level == 0: minimum
    level == 3: verbose
    level > 3 : debugging"""
    global verbosity,logfile,LOGFILENAME
    if level is None:
        return verbosity
    if level < 0:
        # this really should be a class, with the destructor closing the file
        if logfilename is None:
            logfilename = LOGFILENAME
        else:
            LOGFILENAME = logfilename
        try:
            logfile = open(logfilename,'a')
        except:
            raise IOError("Failed to open logfile; provide a filename when level<0.")
    if level > 0 and verbosity < 0:
        # close the logfile if we don't write to it anymore
        try:
            close_log()
        except:
            pass
    verbosity = level
    return verbosity

def close_log():
    """Close open logfile; must be done manually."""
    try:
        logfile.close()
    except AttributeError:
        raise ValueError("no open logfile; use negative verbosity.")

def get_verbosity():
    return verbosity

def msg(level,m=None):
    """msg(level,[m])
       
    1) Print message string if the level <= verbose. level describes the
    priority with lower = more important.

    Terminate string with \\n if a newline is desired or \\r to overwrite
    the current line (eg for output progress indication)

    Note that if the global verbosity level is < 0 then the message is also
    written to the logfile.

    2) If called without a string then msg(level) returns True if it would
    print a message, and False otherwise.
    """
    if level <= abs(verbosity):
        if m:
            print(m), # terminate string with \n if newline desired
            sys.stdout.flush()
            if verbosity < 0:
                logfile.write(m)
        return True
    return False

def fixedwidth_bins(delta,xmin,xmax):
    """Return bins of width delta that cover xmin,xmax (or a larger range).

    dict = fixedwidth_bins(delta,xmin,xmax)

    The dict contains 'Nbins', 'delta', 'min', and 'max'.
    """
    import numpy
    if not numpy.all(xmin < xmax):
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = numpy.asarray(delta,dtype=numpy.float_)
    _xmin = numpy.asarray(xmin,dtype=numpy.float_)
    _xmax = numpy.asarray(xmax,dtype=numpy.float_)
    _length = _xmax - _xmin
    N = numpy.ceil(_length/_delta).astype(numpy.int_)      # number of bins
    dx = 0.5 * (N*_delta - _length)   # add half of the excess to each end
    return {'Nbins':N, 'delta':_delta,'min':_xmin-dx, 'max':_xmax+dx}
    
    

def flatiter(seq):
    """Returns an iterator that flattens a sequence of sequences of sequences...
    (c) 2005 Peter Otten, at  http://www.thescripts.com/forum/thread23631.html
    """
    # avoid infinite recursion with strings
    if isinstance(seq, basestring):
        yield seq
    else:
        try:
            for i in seq:
                for k in flatiter(i):
                    yield k
        except TypeError:
            yield seq

def flatten(x):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]

    From http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
    """
    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def matplotlib_interactive(interactive=False):
    import matplotlib
    if not interactive:
        matplotlib.use('Agg')  # allows running without X11 on compute nodes
    matplotlib.interactive(interactive)            
    return interactive

def iterable(obj):
    """Returns True if obj can be iterated over and is NOT a string."""
    try: len(obj)
    except: return False
    if type(obj) is str:
        return False    # avoid iterating over characters of a string
    return True

def asiterable(obj):
    """Return an object that is an iterable: object itself or wrapepd in a list.

    iterable <-- asiterable(something)
    
    Treats strings as NOT-iterable.
    """
    if not iterable(obj):
        return [obj]
    return obj

def Pearson_r(x,y):
    """Pearson's r (correlation coefficient)

    r = Pearson(x,y)

    x and y are arrays of same length
    
    Historical note:
    Naive implementation of Pearson's r:

    Ex = scipy.stats.mean(x)
    Ey = scipy.stats.mean(y)

    covxy = numpy.sum((x-Ex)*(y-Ey))
    r = covxy/math.sqrt(numpy.sum((x-Ex)**2)*numpy.sum((y-Ey)**2))
    return r
    """
    import numpy
    return numpy.corrcoef(x,y)[1,0]

def linfit(x,y,dy=[]):
    """Fit a straight line y = a + bx to the data in x and y; errors
    on y should be provided in dy in order to assess the goodness of
    the fit and derive errors on the parameters.

    result_dict = linfit(x,y[,dy])

    Fit y = a + bx to the data in x and y by analytically minimizing
    chi^2. dy holds the standard deviations of the individual y_i. If
    dy is not given, they are assumed to be constant (note that in
    this case Q is set to 1 and it is meaningless and chi2 is
    normalised to unit standard deviation on all points!).

    Returns the parameters a and b, their uncertainties sigma_a and
    sigma_b, and their correlation coefficient r_ab; it also returns
    the chi-squared statistic and the goodness-of-fit probability Q
    (that the fit would have chi^2 this large or larger; Q < 10^-2
    indicates that the model is bad --- Q is the probability that a
    value of chi-square as _poor_ as the calculated statistic chi2
    should occur by chance.)

    result_dict =
       intercept, sigma_intercept    a +/- sigma_a
       slope, sigma_slope            b +/- sigma_b
       parameter_correlation         correlation coefficient r_ab
                                     between a and b
       chi_square                    chi^2 test statistic
       Q_fit                         goodness-of-fit probability

    Based on 'Numerical Recipes in C', Ch 15.2.
    """
    import math
    import numpy
    import scipy.stats

    n = len(x)
    m = len(y)
    if n != m:
        raise ValueError("lengths of x and y must match: %s != %s" % (n, m))
    
    try:
        have_dy = (len(dy) > 0)
    except TypeError:
        have_dy = False

    if not have_dy:
        dy = numpy.ones((n),numpy.float)

    x  = numpy.asarray(x)
    y  = numpy.asarray(y)
    dy = numpy.asarray(dy)

    s2  = dy*dy
    S   = numpy.add.reduce(1/s2)
    Sx  = numpy.add.reduce(x/s2)
    Sy  = numpy.add.reduce(y/s2)
    Sxx = numpy.add.reduce(x*x/s2)
    Sxy = numpy.add.reduce(x*y/s2)

    t   = (x - Sx/S)/dy
    Stt = numpy.add.reduce(t*t)

    b = numpy.add.reduce(t*y/dy)/Stt
    a = (Sy - Sx*b)/S

    sa = math.sqrt((1 + (Sx*Sx)/(S*Stt))/S)
    sb = math.sqrt(1/Stt)

    covab = -Sx/(S*Stt)
    r = covab/(sa*sb)

    chi2 = numpy.add.reduce(((y-a-b*x)/dy)**2)
    if not have_dy:
        # estimate error if none were provided
        sigmadata = math.sqrt(chi2/(n-2))
        sa *= sigmadata
        sb *= sigmadata
        Q = 1.0
    else:
        Q = scipy.stats.chisqprob(chi2,n-2)

    return {"intercept":a,"slope":b,
            "sigma_intercept":sa,"sigma_slope":sb,
            "parameter_correlation":r, "chi_square":chi2, "Q":Q}

def autocorrelation_fft(series,include_mean=False,periodic=False,
                        start=None,stop=None,**kwargs):
    """Calculate the auto correlation function.
    
    acf = autocorrelation_fft(series,include_mean=False,**kwargs)

    The time series is correlated with itself across its whole length. It is 0-padded
    and the ACF is corrected for the 0-padding (the values for larger lags are
    increased) unless mode='valid' (see below). 
    Only the [0,len(series)[ interval is returned. The series is normalized to ots 0-th 
    element.

    Note that the series for mode='same'|'full' is inaccurate for long times and
    should probably be truncated at 1/2*len(series). Alternatively, only sample a
    subseries with the stop keyword.

    :Arguments:
    series        (time) series, a 1D numpy array
    include_mean  False: subtract mean(series) from series
    periodic      False: corrected for 0-padding
                  True: return as is
    start,stop    If set, calculate the ACF of series[start:stop] with series;
                  in this case mode='valid' is enforced
    kwargs        keyword arguments for scipy.signal.fftconvolve
                  mode = 'full' | 'same' | 'valid'  (see there)
    """
    import numpy
    import scipy.signal
    kwargs.setdefault('mode','full')
    series = numpy.squeeze(series.astype(float)).copy()   # must copy because de-meaning modifies it
    if not include_mean:
        mean = series.mean()
        series -= mean
    if start or stop:
        kwargs['mode'] = 'valid'   # makes sense to only return the true unpolluted ACF
        start = start or 0
        stop = stop or len(series)
        if start >= stop:
            raise ValueError('Must be start < stop but start = %(start)d >= stop = %(stop)d.' 
                             % locals())

    ac = scipy.signal.fftconvolve(series,series[stop:start:-1,...],**kwargs)

    if kwargs['mode'] == 'valid':
        # origin at start+1
        norm = ac[start+1] or 1.0   # to guard against ACFs of zero arrays
        # Note that this is periodic (and symmetric) over result[0,stop-start+1] and so we 
        # only return one half:
        ##return numpy.concatenate( (ac[start+1:], ac[:start+1]) )[:len(ac)/2] / norm
        # ac is symmetric around start+1 so we average the two halves (just in case):
        ac[:] = numpy.concatenate( (ac[start+1:], ac[:start+1]) ) / norm
        ac = numpy.resize(ac,len(ac)+1)   # make space for replicated 0-th element
        ac[-1] = ac[0]
        if len(ac) % 2 == 1:   
            # orig ac was even
            return 0.5*(ac[:len(ac)/2] + ac[:len(ac)/2:-1])
        else:
            # orig ac was odd: replicate the least important datapoint for second half
            return 0.5*(ac[:len(ac)/2] + ac[:len(ac)/2-1:-1])
    else:
        origin = ac.shape[0]/2        # should work for both odd and even len(series)
        ac = ac[origin:]
        assert len(ac) <= len(series), "Oops: len(ac)=%d  len(series)=%d" % (len(ac),len(series))
        if not periodic:
            ac *= len(series)/(len(series) - 1.0*scipy.arange(len(ac)))   # correct for 0 padding
    norm = ac[0] or 1.0  # to guard against ACFs of zero arrays
    return ac/norm

def averaged_autocorrelation(series,length=None,sliding_window=None,**kwargs):
    """Calculates the averaged ACF of a series.

    mean(acf), std(acf) = averaged_autocorrelation(series,length=None,sliding_window=None):

    Calculate the ACF of a series for only a fraction of the total length, <length> but
    repeat the calculation by setting the origin progressively every <sliding_window>
    steps and average over all the ACFs.

    :Arguments:
    series          time series (by default, mean will be removed)
    length          length (in frames) of the ACF (default: 1/2*len(series))
    sliding_window  repeat ACF calculation every N frames (default: len(series)/100)
    kwargs          additional arguments to autocorrelation_fft()    
    """
    import numpy
    kwargs.pop('start',None) # must filter those kwargs as we set them ourselves
    kwargs.pop('stop',None)
    nframes = len(series)
    length = length or nframes/2
    _length = nframes - length   # _length is the length of the comparison series
    sliding_window = sliding_window or nframes/100
    # note: do NOT be tempted to change nframes-_length to nframes-_length+1
    #       (this will make the last acf 1 step longer, see series[stop:start:-1] !)
    acfs = numpy.array([autocorrelation_fft(series,start=start,stop=start+_length,**kwargs) 
                 for start in xrange(0,nframes-_length,sliding_window)])
    return acfs.mean(axis=0), acfs.std(axis=0)


# Compatibility layer for running in python 2.3


# sorted()
# --- (can also use numpy.sort function instead) ---
try:
    sorted([])
    sorted = sorted
except NameError:
    def sorted(iterable, cmp=None, key=None, reverse=False):
        """sorted(iterable, cmp=None, key=None, reverse=False) --> new sorted list
        Naive pre python 2.4 compatibility fudge.

        With key, cmp must make use of triplets (key,int,value).
        It's a fudge after all. 
        """
        L = list(iterable)
        args = ()
        if cmp is not None:
            args = (cmp,)
        if key is None:
            L.sort(*args)
        else:
            # decorate-sort-undecorate
            deco = [(key(x),i,x) for i,x in enumerate(L)]
            deco.sort(*args)
            L[:] = [y[2] for y in deco] 
        if reverse:
            L.reverse()
        return L

try:
    import collections
    class Fifo(collections.deque):
        pop = collections.deque.popleft

    class Ringbuffer(Fifo):
        """Ring buffer of size capacity; 'pushes' data from left and discards
        on the right.
        """
        #  See http://mail.python.org/pipermail/tutor/2005-March/037149.html.
        def __init__(self,capacity,iterable=None):
            if iterable is None: iterable = []
            super(Ringbuffer,self).__init__(iterable)
            assert capacity > 0
            self.capacity = capacity
            while len(self) > self.capacity:
                super(Ringbuffer,self).pop()   # prune initial loading
        def append(self,x):
            while len(self) >= self.capacity:
                super(Ringbuffer,self).pop()
            super(Ringbuffer,self).append(x)
        def __repr__(self):
            return "Ringbuffer(capacity="+str(self.capacity)+", "+str(list(self))+")"
except ImportError:
    class Ringbuffer(list):
        """Ringbuffer that can be treated as a list. Note that the real queuing
        order is only obtained with the tolist() method.

        Based on
        http://www.onlamp.com/pub/a/python/excerpt/pythonckbk_chap1/index1.html
        """
        def __init__(self, capacity, iterable=None):
            assert capacity > 0
            self.capacity = capacity
            if iterable is None:
                self = []
            else:
                self[:] = list(iterable)[-self.capacity:]

        class __Full(list):
            def append(self, x):
                self[self.cur] = x
                self.cur = (self.cur+1) % self.capacity
            def tolist(self):
                """Return a list of elements from the oldest to the newest."""
                return self[self.cur:] + self[:self.cur]

        def append(self, x):
            super(Ringbuffer,self).append(x)
            if len(self) >= self.capacity:
                self[:] = self[-self.capacity:]
                self.cur = 0
                # Permanently change self's class from non-full to full
                self.__class__ = self.__Full
        def extend(self,iterable):
            for x in list(iterable)[-self.capacity:]:
                self.append(x)
        def tolist(self):
            """Return a list of elements from the oldest to the newest."""
            return self
        def __repr__(self):
            return "Ringbuffer(capacity="+str(self.capacity)+", "+str(list(self))+")"

class DefaultDict(dict):
    """Dictionary based on defaults and updated with keys/values from user."""
    def __init__(self,defaultdict,userdict=None,**kwargs):
        super(DefaultDict,self).__init__(**defaultdict)
        if userdict is not None:
            self.update(userdict)
        self.update(kwargs)

class IntrospectiveDict(dict):
    """A dictionary that contains its keys as attributes for easier introspection.

    Keys that collide with dict methods or attributes are _not_ added as attributes.

    The implementation is simple and certainly not optimized for larger dictionaries
    or ones which are often accessed. Only use it for 'final results' collections
    that you are likely to investigate interactively.

    ARGH: This cannot be pickled safely.
    """
    def __init__(self,*args,**kwargs):
        super(IntrospectiveDict,self).__init__(*args,**kwargs)
        self.__reserved = dict.__dict__.keys()  # don't know how to use super here?
        self._update()

    def __getattr__(self,x):
        self._update()   # maybe something changed? Can't be bothered to subclass update,etc
        try:
            return self.__dict__[x]
        except KeyError:
            return super(IntrospectiveDict,self).__getattr__(x)

    def _update(self):
        for k,v in self.items():  # initialisation of the attributes from the keys
            self._set(k,v)
        
    def _set(self,k,v):
        if k not in self.__reserved:
            self.__setattr__(k,v)

