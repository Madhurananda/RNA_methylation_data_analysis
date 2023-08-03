# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _lrst
else:
    import _lrst

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _lrst.delete_SwigPyIterator

    def value(self):
        return _lrst.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _lrst.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _lrst.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _lrst.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _lrst.SwigPyIterator_equal(self, x)

    def copy(self):
        return _lrst.SwigPyIterator_copy(self)

    def next(self):
        return _lrst.SwigPyIterator_next(self)

    def __next__(self):
        return _lrst.SwigPyIterator___next__(self)

    def previous(self):
        return _lrst.SwigPyIterator_previous(self)

    def advance(self, n):
        return _lrst.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _lrst.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _lrst.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _lrst.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _lrst.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _lrst.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _lrst.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _lrst:
_lrst.SwigPyIterator_swigregister(SwigPyIterator)

class IntVector(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _lrst.IntVector_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _lrst.IntVector___nonzero__(self)

    def __bool__(self):
        return _lrst.IntVector___bool__(self)

    def __len__(self):
        return _lrst.IntVector___len__(self)

    def __getslice__(self, i, j):
        return _lrst.IntVector___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _lrst.IntVector___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _lrst.IntVector___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _lrst.IntVector___delitem__(self, *args)

    def __getitem__(self, *args):
        return _lrst.IntVector___getitem__(self, *args)

    def __setitem__(self, *args):
        return _lrst.IntVector___setitem__(self, *args)

    def pop(self):
        return _lrst.IntVector_pop(self)

    def append(self, x):
        return _lrst.IntVector_append(self, x)

    def empty(self):
        return _lrst.IntVector_empty(self)

    def size(self):
        return _lrst.IntVector_size(self)

    def swap(self, v):
        return _lrst.IntVector_swap(self, v)

    def begin(self):
        return _lrst.IntVector_begin(self)

    def end(self):
        return _lrst.IntVector_end(self)

    def rbegin(self):
        return _lrst.IntVector_rbegin(self)

    def rend(self):
        return _lrst.IntVector_rend(self)

    def clear(self):
        return _lrst.IntVector_clear(self)

    def get_allocator(self):
        return _lrst.IntVector_get_allocator(self)

    def pop_back(self):
        return _lrst.IntVector_pop_back(self)

    def erase(self, *args):
        return _lrst.IntVector_erase(self, *args)

    def __init__(self, *args):
        _lrst.IntVector_swiginit(self, _lrst.new_IntVector(*args))

    def push_back(self, x):
        return _lrst.IntVector_push_back(self, x)

    def front(self):
        return _lrst.IntVector_front(self)

    def back(self):
        return _lrst.IntVector_back(self)

    def assign(self, n, x):
        return _lrst.IntVector_assign(self, n, x)

    def resize(self, *args):
        return _lrst.IntVector_resize(self, *args)

    def insert(self, *args):
        return _lrst.IntVector_insert(self, *args)

    def reserve(self, n):
        return _lrst.IntVector_reserve(self, n)

    def capacity(self):
        return _lrst.IntVector_capacity(self)
    __swig_destroy__ = _lrst.delete_IntVector

# Register IntVector in _lrst:
_lrst.IntVector_swigregister(IntVector)

class DoubleVector(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _lrst.DoubleVector_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _lrst.DoubleVector___nonzero__(self)

    def __bool__(self):
        return _lrst.DoubleVector___bool__(self)

    def __len__(self):
        return _lrst.DoubleVector___len__(self)

    def __getslice__(self, i, j):
        return _lrst.DoubleVector___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _lrst.DoubleVector___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _lrst.DoubleVector___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _lrst.DoubleVector___delitem__(self, *args)

    def __getitem__(self, *args):
        return _lrst.DoubleVector___getitem__(self, *args)

    def __setitem__(self, *args):
        return _lrst.DoubleVector___setitem__(self, *args)

    def pop(self):
        return _lrst.DoubleVector_pop(self)

    def append(self, x):
        return _lrst.DoubleVector_append(self, x)

    def empty(self):
        return _lrst.DoubleVector_empty(self)

    def size(self):
        return _lrst.DoubleVector_size(self)

    def swap(self, v):
        return _lrst.DoubleVector_swap(self, v)

    def begin(self):
        return _lrst.DoubleVector_begin(self)

    def end(self):
        return _lrst.DoubleVector_end(self)

    def rbegin(self):
        return _lrst.DoubleVector_rbegin(self)

    def rend(self):
        return _lrst.DoubleVector_rend(self)

    def clear(self):
        return _lrst.DoubleVector_clear(self)

    def get_allocator(self):
        return _lrst.DoubleVector_get_allocator(self)

    def pop_back(self):
        return _lrst.DoubleVector_pop_back(self)

    def erase(self, *args):
        return _lrst.DoubleVector_erase(self, *args)

    def __init__(self, *args):
        _lrst.DoubleVector_swiginit(self, _lrst.new_DoubleVector(*args))

    def push_back(self, x):
        return _lrst.DoubleVector_push_back(self, x)

    def front(self):
        return _lrst.DoubleVector_front(self)

    def back(self):
        return _lrst.DoubleVector_back(self)

    def assign(self, n, x):
        return _lrst.DoubleVector_assign(self, n, x)

    def resize(self, *args):
        return _lrst.DoubleVector_resize(self, *args)

    def insert(self, *args):
        return _lrst.DoubleVector_insert(self, *args)

    def reserve(self, n):
        return _lrst.DoubleVector_reserve(self, n)

    def capacity(self):
        return _lrst.DoubleVector_capacity(self)
    __swig_destroy__ = _lrst.delete_DoubleVector

# Register DoubleVector in _lrst:
_lrst.DoubleVector_swigregister(DoubleVector)

class Int64Vector(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _lrst.Int64Vector_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _lrst.Int64Vector___nonzero__(self)

    def __bool__(self):
        return _lrst.Int64Vector___bool__(self)

    def __len__(self):
        return _lrst.Int64Vector___len__(self)

    def __getslice__(self, i, j):
        return _lrst.Int64Vector___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _lrst.Int64Vector___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _lrst.Int64Vector___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _lrst.Int64Vector___delitem__(self, *args)

    def __getitem__(self, *args):
        return _lrst.Int64Vector___getitem__(self, *args)

    def __setitem__(self, *args):
        return _lrst.Int64Vector___setitem__(self, *args)

    def pop(self):
        return _lrst.Int64Vector_pop(self)

    def append(self, x):
        return _lrst.Int64Vector_append(self, x)

    def empty(self):
        return _lrst.Int64Vector_empty(self)

    def size(self):
        return _lrst.Int64Vector_size(self)

    def swap(self, v):
        return _lrst.Int64Vector_swap(self, v)

    def begin(self):
        return _lrst.Int64Vector_begin(self)

    def end(self):
        return _lrst.Int64Vector_end(self)

    def rbegin(self):
        return _lrst.Int64Vector_rbegin(self)

    def rend(self):
        return _lrst.Int64Vector_rend(self)

    def clear(self):
        return _lrst.Int64Vector_clear(self)

    def get_allocator(self):
        return _lrst.Int64Vector_get_allocator(self)

    def pop_back(self):
        return _lrst.Int64Vector_pop_back(self)

    def erase(self, *args):
        return _lrst.Int64Vector_erase(self, *args)

    def __init__(self, *args):
        _lrst.Int64Vector_swiginit(self, _lrst.new_Int64Vector(*args))

    def push_back(self, x):
        return _lrst.Int64Vector_push_back(self, x)

    def front(self):
        return _lrst.Int64Vector_front(self)

    def back(self):
        return _lrst.Int64Vector_back(self)

    def assign(self, n, x):
        return _lrst.Int64Vector_assign(self, n, x)

    def resize(self, *args):
        return _lrst.Int64Vector_resize(self, *args)

    def insert(self, *args):
        return _lrst.Int64Vector_insert(self, *args)

    def reserve(self, n):
        return _lrst.Int64Vector_reserve(self, n)

    def capacity(self):
        return _lrst.Int64Vector_capacity(self)
    __swig_destroy__ = _lrst.delete_Int64Vector

# Register Int64Vector in _lrst:
_lrst.Int64Vector_swigregister(Int64Vector)

class Int2DVector(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _lrst.Int2DVector_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _lrst.Int2DVector___nonzero__(self)

    def __bool__(self):
        return _lrst.Int2DVector___bool__(self)

    def __len__(self):
        return _lrst.Int2DVector___len__(self)

    def __getslice__(self, i, j):
        return _lrst.Int2DVector___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _lrst.Int2DVector___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _lrst.Int2DVector___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _lrst.Int2DVector___delitem__(self, *args)

    def __getitem__(self, *args):
        return _lrst.Int2DVector___getitem__(self, *args)

    def __setitem__(self, *args):
        return _lrst.Int2DVector___setitem__(self, *args)

    def pop(self):
        return _lrst.Int2DVector_pop(self)

    def append(self, x):
        return _lrst.Int2DVector_append(self, x)

    def empty(self):
        return _lrst.Int2DVector_empty(self)

    def size(self):
        return _lrst.Int2DVector_size(self)

    def swap(self, v):
        return _lrst.Int2DVector_swap(self, v)

    def begin(self):
        return _lrst.Int2DVector_begin(self)

    def end(self):
        return _lrst.Int2DVector_end(self)

    def rbegin(self):
        return _lrst.Int2DVector_rbegin(self)

    def rend(self):
        return _lrst.Int2DVector_rend(self)

    def clear(self):
        return _lrst.Int2DVector_clear(self)

    def get_allocator(self):
        return _lrst.Int2DVector_get_allocator(self)

    def pop_back(self):
        return _lrst.Int2DVector_pop_back(self)

    def erase(self, *args):
        return _lrst.Int2DVector_erase(self, *args)

    def __init__(self, *args):
        _lrst.Int2DVector_swiginit(self, _lrst.new_Int2DVector(*args))

    def push_back(self, x):
        return _lrst.Int2DVector_push_back(self, x)

    def front(self):
        return _lrst.Int2DVector_front(self)

    def back(self):
        return _lrst.Int2DVector_back(self)

    def assign(self, n, x):
        return _lrst.Int2DVector_assign(self, n, x)

    def resize(self, *args):
        return _lrst.Int2DVector_resize(self, *args)

    def insert(self, *args):
        return _lrst.Int2DVector_insert(self, *args)

    def reserve(self, n):
        return _lrst.Int2DVector_reserve(self, n)

    def capacity(self):
        return _lrst.Int2DVector_capacity(self)
    __swig_destroy__ = _lrst.delete_Int2DVector

# Register Int2DVector in _lrst:
_lrst.Int2DVector_swigregister(Int2DVector)

MAX_INPUT_FILES = _lrst.MAX_INPUT_FILES
class Input_Para(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    threads = property(_lrst.Input_Para_threads_get, _lrst.Input_Para_threads_set)
    num_input_files = property(_lrst.Input_Para_num_input_files_get, _lrst.Input_Para_num_input_files_set)
    out_prefix = property(_lrst.Input_Para_out_prefix_get, _lrst.Input_Para_out_prefix_set)
    other_flags = property(_lrst.Input_Para_other_flags_get, _lrst.Input_Para_other_flags_set)
    rdm_seed = property(_lrst.Input_Para_rdm_seed_get, _lrst.Input_Para_rdm_seed_set)
    downsample_percentage = property(_lrst.Input_Para_downsample_percentage_get, _lrst.Input_Para_downsample_percentage_set)
    user_defined_fastq_base_qual_offset = property(_lrst.Input_Para_user_defined_fastq_base_qual_offset_get, _lrst.Input_Para_user_defined_fastq_base_qual_offset_set)
    output_folder = property(_lrst.Input_Para_output_folder_get, _lrst.Input_Para_output_folder_set)
    input_files = property(_lrst.Input_Para_input_files_get, _lrst.Input_Para_input_files_set)

    def add_input_file(self, _ip_file):
        return _lrst.Input_Para_add_input_file(self, _ip_file)

    def __init__(self):
        _lrst.Input_Para_swiginit(self, _lrst.new_Input_Para())
    __swig_destroy__ = _lrst.delete_Input_Para

# Register Input_Para in _lrst:
_lrst.Input_Para_swigregister(Input_Para)

MAX_READ_LENGTH = _lrst.MAX_READ_LENGTH
MAX_MAP_QUALITY = _lrst.MAX_MAP_QUALITY
MAX_BASE_QUALITY = _lrst.MAX_BASE_QUALITY
MAX_READ_QUALITY = _lrst.MAX_READ_QUALITY
MAX_SIGNAL_VALUE = _lrst.MAX_SIGNAL_VALUE
PERCENTAGE_ARRAY_SIZE = _lrst.PERCENTAGE_ARRAY_SIZE
ZeroDefault = _lrst.ZeroDefault
MoneDefault = _lrst.MoneDefault
class Output_Info(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    error_flag = property(_lrst.Output_Info_error_flag_get, _lrst.Output_Info_error_flag_set)
    error_str = property(_lrst.Output_Info_error_str_get, _lrst.Output_Info_error_str_set)

    def __init__(self):
        _lrst.Output_Info_swiginit(self, _lrst.new_Output_Info())
    __swig_destroy__ = _lrst.delete_Output_Info

# Register Output_Info in _lrst:
_lrst.Output_Info_swigregister(Output_Info)

class Basic_Seq_Statistics(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    total_num_reads = property(_lrst.Basic_Seq_Statistics_total_num_reads_get, _lrst.Basic_Seq_Statistics_total_num_reads_set)
    total_num_bases = property(_lrst.Basic_Seq_Statistics_total_num_bases_get, _lrst.Basic_Seq_Statistics_total_num_bases_set)
    longest_read_length = property(_lrst.Basic_Seq_Statistics_longest_read_length_get, _lrst.Basic_Seq_Statistics_longest_read_length_set)
    n50_read_length = property(_lrst.Basic_Seq_Statistics_n50_read_length_get, _lrst.Basic_Seq_Statistics_n50_read_length_set)
    n95_read_length = property(_lrst.Basic_Seq_Statistics_n95_read_length_get, _lrst.Basic_Seq_Statistics_n95_read_length_set)
    n05_read_length = property(_lrst.Basic_Seq_Statistics_n05_read_length_get, _lrst.Basic_Seq_Statistics_n05_read_length_set)
    mean_read_length = property(_lrst.Basic_Seq_Statistics_mean_read_length_get, _lrst.Basic_Seq_Statistics_mean_read_length_set)
    nx_read_length = property(_lrst.Basic_Seq_Statistics_nx_read_length_get, _lrst.Basic_Seq_Statistics_nx_read_length_set)
    NXX_read_length = property(_lrst.Basic_Seq_Statistics_NXX_read_length_get, _lrst.Basic_Seq_Statistics_NXX_read_length_set)
    median_read_length = property(_lrst.Basic_Seq_Statistics_median_read_length_get, _lrst.Basic_Seq_Statistics_median_read_length_set)
    total_a_cnt = property(_lrst.Basic_Seq_Statistics_total_a_cnt_get, _lrst.Basic_Seq_Statistics_total_a_cnt_set)
    total_c_cnt = property(_lrst.Basic_Seq_Statistics_total_c_cnt_get, _lrst.Basic_Seq_Statistics_total_c_cnt_set)
    total_g_cnt = property(_lrst.Basic_Seq_Statistics_total_g_cnt_get, _lrst.Basic_Seq_Statistics_total_g_cnt_set)
    total_tu_cnt = property(_lrst.Basic_Seq_Statistics_total_tu_cnt_get, _lrst.Basic_Seq_Statistics_total_tu_cnt_set)
    total_n_cnt = property(_lrst.Basic_Seq_Statistics_total_n_cnt_get, _lrst.Basic_Seq_Statistics_total_n_cnt_set)
    gc_cnt = property(_lrst.Basic_Seq_Statistics_gc_cnt_get, _lrst.Basic_Seq_Statistics_gc_cnt_set)
    read_length_count = property(_lrst.Basic_Seq_Statistics_read_length_count_get, _lrst.Basic_Seq_Statistics_read_length_count_set)
    read_gc_content_count = property(_lrst.Basic_Seq_Statistics_read_gc_content_count_get, _lrst.Basic_Seq_Statistics_read_gc_content_count_set)
    read_lengths = property(_lrst.Basic_Seq_Statistics_read_lengths_get, _lrst.Basic_Seq_Statistics_read_lengths_set)

    def reset(self):
        return _lrst.Basic_Seq_Statistics_reset(self)

    def add(self, t_seq_st):
        return _lrst.Basic_Seq_Statistics_add(self, t_seq_st)

    def global_sum(self):
        return _lrst.Basic_Seq_Statistics_global_sum(self)

    def global_sum_no_gc(self):
        return _lrst.Basic_Seq_Statistics_global_sum_no_gc(self)

    def resize(self):
        return _lrst.Basic_Seq_Statistics_resize(self)

    def __init__(self, *args):
        _lrst.Basic_Seq_Statistics_swiginit(self, _lrst.new_Basic_Seq_Statistics(*args))
    __swig_destroy__ = _lrst.delete_Basic_Seq_Statistics

# Register Basic_Seq_Statistics in _lrst:
_lrst.Basic_Seq_Statistics_swigregister(Basic_Seq_Statistics)

class Basic_Seq_Quality_Statistics(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    base_quality_distribution = property(_lrst.Basic_Seq_Quality_Statistics_base_quality_distribution_get, _lrst.Basic_Seq_Quality_Statistics_base_quality_distribution_set)
    read_average_base_quality_distribution = property(_lrst.Basic_Seq_Quality_Statistics_read_average_base_quality_distribution_get, _lrst.Basic_Seq_Quality_Statistics_read_average_base_quality_distribution_set)
    min_base_quality = property(_lrst.Basic_Seq_Quality_Statistics_min_base_quality_get, _lrst.Basic_Seq_Quality_Statistics_min_base_quality_set)
    max_base_quality = property(_lrst.Basic_Seq_Quality_Statistics_max_base_quality_get, _lrst.Basic_Seq_Quality_Statistics_max_base_quality_set)
    pos_quality_distribution = property(_lrst.Basic_Seq_Quality_Statistics_pos_quality_distribution_get, _lrst.Basic_Seq_Quality_Statistics_pos_quality_distribution_set)
    pos_quality_distribution_dev = property(_lrst.Basic_Seq_Quality_Statistics_pos_quality_distribution_dev_get, _lrst.Basic_Seq_Quality_Statistics_pos_quality_distribution_dev_set)
    pos_quality_distribution_count = property(_lrst.Basic_Seq_Quality_Statistics_pos_quality_distribution_count_get, _lrst.Basic_Seq_Quality_Statistics_pos_quality_distribution_count_set)
    max_length = property(_lrst.Basic_Seq_Quality_Statistics_max_length_get, _lrst.Basic_Seq_Quality_Statistics_max_length_set)
    read_quality_distribution = property(_lrst.Basic_Seq_Quality_Statistics_read_quality_distribution_get, _lrst.Basic_Seq_Quality_Statistics_read_quality_distribution_set)
    min_read_quality = property(_lrst.Basic_Seq_Quality_Statistics_min_read_quality_get, _lrst.Basic_Seq_Quality_Statistics_min_read_quality_set)
    max_read_quality = property(_lrst.Basic_Seq_Quality_Statistics_max_read_quality_get, _lrst.Basic_Seq_Quality_Statistics_max_read_quality_set)

    def reset(self):
        return _lrst.Basic_Seq_Quality_Statistics_reset(self)

    def add(self, t_qual_st):
        return _lrst.Basic_Seq_Quality_Statistics_add(self, t_qual_st)

    def global_sum(self):
        return _lrst.Basic_Seq_Quality_Statistics_global_sum(self)

    def __init__(self, *args):
        _lrst.Basic_Seq_Quality_Statistics_swiginit(self, _lrst.new_Basic_Seq_Quality_Statistics(*args))
    __swig_destroy__ = _lrst.delete_Basic_Seq_Quality_Statistics

# Register Basic_Seq_Quality_Statistics in _lrst:
_lrst.Basic_Seq_Quality_Statistics_swigregister(Basic_Seq_Quality_Statistics)

class Output_FA(Output_Info):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    long_read_info = property(_lrst.Output_FA_long_read_info_get, _lrst.Output_FA_long_read_info_set)

    def __init__(self):
        _lrst.Output_FA_swiginit(self, _lrst.new_Output_FA())
    __swig_destroy__ = _lrst.delete_Output_FA

# Register Output_FA in _lrst:
_lrst.Output_FA_swigregister(Output_FA)

class Output_FQ(Output_FA):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    seq_quality_info = property(_lrst.Output_FQ_seq_quality_info_get, _lrst.Output_FQ_seq_quality_info_set)

    def __init__(self):
        _lrst.Output_FQ_swiginit(self, _lrst.new_Output_FQ())
    __swig_destroy__ = _lrst.delete_Output_FQ

# Register Output_FQ in _lrst:
_lrst.Output_FQ_swigregister(Output_FQ)

class Output_BAM(Output_FQ):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    num_primary_alignment = property(_lrst.Output_BAM_num_primary_alignment_get, _lrst.Output_BAM_num_primary_alignment_set)
    num_secondary_alignment = property(_lrst.Output_BAM_num_secondary_alignment_get, _lrst.Output_BAM_num_secondary_alignment_set)
    num_reads_with_secondary_alignment = property(_lrst.Output_BAM_num_reads_with_secondary_alignment_get, _lrst.Output_BAM_num_reads_with_secondary_alignment_set)
    num_supplementary_alignment = property(_lrst.Output_BAM_num_supplementary_alignment_get, _lrst.Output_BAM_num_supplementary_alignment_set)
    num_reads_with_supplementary_alignment = property(_lrst.Output_BAM_num_reads_with_supplementary_alignment_get, _lrst.Output_BAM_num_reads_with_supplementary_alignment_set)
    num_reads_with_both_secondary_supplementary_alignment = property(_lrst.Output_BAM_num_reads_with_both_secondary_supplementary_alignment_get, _lrst.Output_BAM_num_reads_with_both_secondary_supplementary_alignment_set)
    forward_alignment = property(_lrst.Output_BAM_forward_alignment_get, _lrst.Output_BAM_forward_alignment_set)
    reverse_alignment = property(_lrst.Output_BAM_reverse_alignment_get, _lrst.Output_BAM_reverse_alignment_set)
    map_quality_distribution = property(_lrst.Output_BAM_map_quality_distribution_get, _lrst.Output_BAM_map_quality_distribution_set)
    min_map_quality = property(_lrst.Output_BAM_min_map_quality_get, _lrst.Output_BAM_min_map_quality_set)
    max_map_quality = property(_lrst.Output_BAM_max_map_quality_get, _lrst.Output_BAM_max_map_quality_set)
    num_matched_bases = property(_lrst.Output_BAM_num_matched_bases_get, _lrst.Output_BAM_num_matched_bases_set)
    num_mismatched_bases = property(_lrst.Output_BAM_num_mismatched_bases_get, _lrst.Output_BAM_num_mismatched_bases_set)
    num_ins_bases = property(_lrst.Output_BAM_num_ins_bases_get, _lrst.Output_BAM_num_ins_bases_set)
    num_del_bases = property(_lrst.Output_BAM_num_del_bases_get, _lrst.Output_BAM_num_del_bases_set)
    num_clip_bases = property(_lrst.Output_BAM_num_clip_bases_get, _lrst.Output_BAM_num_clip_bases_set)
    accuracy_per_read = property(_lrst.Output_BAM_accuracy_per_read_get, _lrst.Output_BAM_accuracy_per_read_set)
    mapped_long_read_info = property(_lrst.Output_BAM_mapped_long_read_info_get, _lrst.Output_BAM_mapped_long_read_info_set)
    unmapped_long_read_info = property(_lrst.Output_BAM_unmapped_long_read_info_get, _lrst.Output_BAM_unmapped_long_read_info_set)
    mapped_seq_quality_info = property(_lrst.Output_BAM_mapped_seq_quality_info_get, _lrst.Output_BAM_mapped_seq_quality_info_set)
    unmapped_seq_quality_info = property(_lrst.Output_BAM_unmapped_seq_quality_info_get, _lrst.Output_BAM_unmapped_seq_quality_info_set)

    def reset(self):
        return _lrst.Output_BAM_reset(self)

    def add(self, t_output_bam):
        return _lrst.Output_BAM_add(self, t_output_bam)

    def global_sum(self):
        return _lrst.Output_BAM_global_sum(self)

    def __init__(self):
        _lrst.Output_BAM_swiginit(self, _lrst.new_Output_BAM())
    __swig_destroy__ = _lrst.delete_Output_BAM

# Register Output_BAM in _lrst:
_lrst.Output_BAM_swigregister(Output_BAM)

class Basic_SeqTxt_Statistics(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    long_read_info = property(_lrst.Basic_SeqTxt_Statistics_long_read_info_get, _lrst.Basic_SeqTxt_Statistics_long_read_info_set)
    seq_quality_info = property(_lrst.Basic_SeqTxt_Statistics_seq_quality_info_get, _lrst.Basic_SeqTxt_Statistics_seq_quality_info_set)
    signal_range = property(_lrst.Basic_SeqTxt_Statistics_signal_range_get, _lrst.Basic_SeqTxt_Statistics_signal_range_set)
    min_signal = property(_lrst.Basic_SeqTxt_Statistics_min_signal_get, _lrst.Basic_SeqTxt_Statistics_min_signal_set)
    max_signal = property(_lrst.Basic_SeqTxt_Statistics_max_signal_get, _lrst.Basic_SeqTxt_Statistics_max_signal_set)

    def reset(self):
        return _lrst.Basic_SeqTxt_Statistics_reset(self)

    def add(self, output_data):
        return _lrst.Basic_SeqTxt_Statistics_add(self, output_data)

    def global_sum(self):
        return _lrst.Basic_SeqTxt_Statistics_global_sum(self)

    def __init__(self):
        _lrst.Basic_SeqTxt_Statistics_swiginit(self, _lrst.new_Basic_SeqTxt_Statistics())
    __swig_destroy__ = _lrst.delete_Basic_SeqTxt_Statistics

# Register Basic_SeqTxt_Statistics in _lrst:
_lrst.Basic_SeqTxt_Statistics_swigregister(Basic_SeqTxt_Statistics)

class Output_SeqTxt(Output_Info):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    all_long_read_info = property(_lrst.Output_SeqTxt_all_long_read_info_get, _lrst.Output_SeqTxt_all_long_read_info_set)
    passed_long_read_info = property(_lrst.Output_SeqTxt_passed_long_read_info_get, _lrst.Output_SeqTxt_passed_long_read_info_set)
    failed_long_read_info = property(_lrst.Output_SeqTxt_failed_long_read_info_get, _lrst.Output_SeqTxt_failed_long_read_info_set)

    def reset(self):
        return _lrst.Output_SeqTxt_reset(self)

    def add(self, output_data):
        return _lrst.Output_SeqTxt_add(self, output_data)

    def global_sum(self):
        return _lrst.Output_SeqTxt_global_sum(self)

    def __init__(self):
        _lrst.Output_SeqTxt_swiginit(self, _lrst.new_Output_SeqTxt())
    __swig_destroy__ = _lrst.delete_Output_SeqTxt

# Register Output_SeqTxt in _lrst:
_lrst.Output_SeqTxt_swigregister(Output_SeqTxt)

class Base_Signals(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    read_name = property(_lrst.Base_Signals_read_name_get, _lrst.Base_Signals_read_name_set)
    base_count = property(_lrst.Base_Signals_base_count_get, _lrst.Base_Signals_base_count_set)
    sequence_data_str = property(_lrst.Base_Signals_sequence_data_str_get, _lrst.Base_Signals_sequence_data_str_set)
    basecall_signals = property(_lrst.Base_Signals_basecall_signals_get, _lrst.Base_Signals_basecall_signals_set)

    def getBaseCount(self):
        return _lrst.Base_Signals_getBaseCount(self)

    def getReadName(self):
        return _lrst.Base_Signals_getReadName(self)

    def getSequenceString(self):
        return _lrst.Base_Signals_getSequenceString(self)

    def getDataVector(self):
        return _lrst.Base_Signals_getDataVector(self)

    def __init__(self, read_name, sequence_data_str, basecall_signals):
        _lrst.Base_Signals_swiginit(self, _lrst.new_Base_Signals(read_name, sequence_data_str, basecall_signals))
    __swig_destroy__ = _lrst.delete_Base_Signals

# Register Base_Signals in _lrst:
_lrst.Base_Signals_swigregister(Base_Signals)

class Output_FAST5(Output_FA):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    seq_quality_info = property(_lrst.Output_FAST5_seq_quality_info_get, _lrst.Output_FAST5_seq_quality_info_set)
    read_count = property(_lrst.Output_FAST5_read_count_get, _lrst.Output_FAST5_read_count_set)
    base_count = property(_lrst.Output_FAST5_base_count_get, _lrst.Output_FAST5_base_count_set)
    read_base_signals = property(_lrst.Output_FAST5_read_base_signals_get, _lrst.Output_FAST5_read_base_signals_set)

    def getReadCount(self):
        return _lrst.Output_FAST5_getReadCount(self)

    def getTotalBaseCount(self):
        return _lrst.Output_FAST5_getTotalBaseCount(self)

    def getNthReadName(self, read_index):
        return _lrst.Output_FAST5_getNthReadName(self, read_index)

    def getNthReadSequence(self, read_index):
        return _lrst.Output_FAST5_getNthReadSequence(self, read_index)

    def addReadBaseSignals(self, values):
        return _lrst.Output_FAST5_addReadBaseSignals(self, values)

    def getNthReadBaseSignals(self, read_index):
        return _lrst.Output_FAST5_getNthReadBaseSignals(self, read_index)

    def getNthReadBaseMeans(self, read_index):
        return _lrst.Output_FAST5_getNthReadBaseMeans(self, read_index)

    def getNthReadBaseStds(self, read_index):
        return _lrst.Output_FAST5_getNthReadBaseStds(self, read_index)

    def getNthReadBaseMedians(self, read_index):
        return _lrst.Output_FAST5_getNthReadBaseMedians(self, read_index)

    def getNthReadPearsonSkewnessCoeff(self, read_index):
        return _lrst.Output_FAST5_getNthReadPearsonSkewnessCoeff(self, read_index)

    def getNthReadKurtosis(self, read_index):
        return _lrst.Output_FAST5_getNthReadKurtosis(self, read_index)

    def __init__(self):
        _lrst.Output_FAST5_swiginit(self, _lrst.new_Output_FAST5())
    __swig_destroy__ = _lrst.delete_Output_FAST5

# Register Output_FAST5 in _lrst:
_lrst.Output_FAST5_swigregister(Output_FAST5)


def callBAMModule(_input_data, py_output_bam):
    return _lrst.callBAMModule(_input_data, py_output_bam)

def callFASTQModule(_input_data, py_output_fq):
    return _lrst.callFASTQModule(_input_data, py_output_fq)

def callFASTAModule(_input_data, py_output_fa):
    return _lrst.callFASTAModule(_input_data, py_output_fa)

def callSeqTxtModule(_input_data, py_output_seqtxt):
    return _lrst.callSeqTxtModule(_input_data, py_output_seqtxt)

def callFAST5Module(_input_data, py_output_FAST5):
    return _lrst.callFAST5Module(_input_data, py_output_FAST5)


