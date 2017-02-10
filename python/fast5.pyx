from cython.operator cimport dereference as deref

from libc.stdint cimport int16_t
from libcpp cimport bool as cbool
from libcpp.map cimport map as cmap
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "fast5.hpp" namespace "fast5":

    cdef string cpp_version "fast5::version"

    ctypedef cmap[string, string] Attr_Map

    struct Channel_Id_Params:
        string channel_number
        double digitisation
        double offset
        double range
        double sampling_rate

    ctypedef Attr_Map Tracking_Id_Params

    ctypedef Attr_Map Sequences_Params

    struct Raw_Samples_Params:
        string read_id
        long long read_number
        long long start_mux
        long long start_time
        long long duration

    ctypedef float Raw_Sample

    ctypedef int16_t Raw_Int_Sample

    struct EventDetection_Events_Params:
        string read_id
        long long read_number
        long long scaling_used
        long long start_mux
        long long start_time
        long long duration
        double median_before
        unsigned abasic_found

    struct EventDetection_Event:
        double mean
        double stdv
        long long start
        long long length

    struct Basecall_Model_Params:
        double scale
        double shift
        double drift
        double var
        double scale_sd
        double var_sd

    struct Basecall_Model_State:
        long long variant
        double level_mean
        double level_stdv
        double sd_mean
        double sd_stdv
        double weight
        #char kmer[8]

    struct Basecall_Events_Params:
        double start_time
        double duration

    struct Basecall_Event:
        double mean
        double stdv
        double start
        double length
        double p_model_state
        double p_mp_state
        double p_A
        double p_C
        double p_G
        double p_T
        long long move
        #char model_state[8]
        #char mp_state[8]

    struct Basecall_Alignment_Entry:
        long long template_index
        long long complement_index
        #char kmer[8]

    cppclass Cpp_File "fast5::File":
        Cpp_File()
        Cpp_File(string)
        Cpp_File(string, cbool)

        cbool is_open()
        cbool is_rw()
        string file_name()
        void open(string)
        void open(string, cbool)
        void create(string)
        void create(string, cbool)
        void close()
        @staticmethod
        cbool is_valid_file(string)

        string file_version()

        cbool have_channel_id_params()
        Channel_Id_Params get_channel_id_params()
        cbool have_sampling_rate()
        double get_sampling_rate()

        cbool have_tracking_id_params()
        Tracking_Id_Params get_tracking_id_params()

        cbool have_sequences_params()
        Sequences_Params get_sequences_params()

        vector[string] get_raw_samples_read_name_list()
        cbool have_raw_samples()
        cbool have_raw_samples(string)
        Raw_Samples_Params get_raw_samples_params()
        Raw_Samples_Params get_raw_samples_params(string)
        vector[Raw_Int_Sample] get_raw_int_samples()
        vector[Raw_Int_Sample] get_raw_int_samples(string)
        vector[Raw_Sample] get_raw_samples()
        vector[Raw_Sample] get_raw_samples(string)

        vector[string] get_eventdetection_group_list()
        cbool have_eventdetection_group()
        cbool have_eventdetection_group(string)
        vector[string] get_eventdetection_read_name_list()
        vector[string] get_eventdetection_read_name_list(string)
        cbool have_eventdetection_events()
        cbool have_eventdetection_events(string)
        cbool have_eventdetection_events(string, string)
        Attr_Map get_eventdetection_params()
        Attr_Map get_eventdetection_params(string)
        EventDetection_Events_Params get_eventdetection_events_params()
        EventDetection_Events_Params get_eventdetection_events_params(string)
        EventDetection_Events_Params get_eventdetection_events_params(string, string)
        vector[EventDetection_Event] get_eventdetection_events()
        vector[EventDetection_Event] get_eventdetection_events(string)
        vector[EventDetection_Event] get_eventdetection_events(string, string)

        vector[string] get_basecall_group_list()
        cbool have_basecall_group()
        cbool have_basecall_group(string)
        vector[string] get_basecall_strand_group_list(unsigned)
        cbool have_basecall_strand_group(unsigned)
        cbool have_basecall_strand_group(unsigned, string)
        string get_basecall_1d_group(string)
        string get_basecall_eventdetection_group(string)
        Attr_Map get_basecall_params(string)
        cbool have_basecall_log(string)
        string get_basecall_log(string)

        cbool have_basecall_fastq(unsigned)
        cbool have_basecall_fastq(unsigned, string)
        string get_basecall_fastq(unsigned)
        string get_basecall_fastq(unsigned, string)
        cbool have_basecall_seq(unsigned)
        cbool have_basecall_seq(unsigned, string)
        string get_basecall_seq(unsigned)
        string get_basecall_seq(unsigned, string)

        cbool have_basecall_model(unsigned)
        cbool have_basecall_model(unsigned, string)
        string get_basecall_model_file(unsigned)
        string get_basecall_model_file(unsigned, string)
        Basecall_Model_Params get_basecall_model_params(unsigned)
        Basecall_Model_Params get_basecall_model_params(unsigned, string)
        vector[Basecall_Model_State] get_basecall_model(unsigned)
        vector[Basecall_Model_State] get_basecall_model(unsigned, string)

        cbool have_basecall_events(unsigned)
        cbool have_basecall_events(unsigned, string)
        Basecall_Events_Params get_basecall_events_params(unsigned)
        Basecall_Events_Params get_basecall_events_params(unsigned, string)
        vector[Basecall_Event] get_basecall_events(unsigned)
        vector[Basecall_Event] get_basecall_events(unsigned, string)

        cbool have_basecall_alignment()
        cbool have_basecall_alignment(string)
        vector[Basecall_Alignment_Entry] get_basecall_alignment()
        vector[Basecall_Alignment_Entry] get_basecall_alignment(string)

__version__ = cpp_version

cdef class File:
    cdef unique_ptr[Cpp_File] thisptr

    def __init__(self, name=None, rw=None):
        if name is None:
            self.thisptr.reset(new Cpp_File())
        elif rw is None:
            self.thisptr.reset(new Cpp_File(name))
        else:
            self.thisptr.reset(new Cpp_File(name, rw))

    def is_open(self):
        return deref(self.thisptr).is_open()
    def is_rw(self):
        return deref(self.thisptr).is_rw()
    def file_name(self):
        return deref(self.thisptr).file_name()
    def open(self, file_name, rw=None):
        if rw is None:
            return deref(self.thisptr).open(file_name)
        else:
            return deref(self.thisptr).open(file_name, rw)
    def create(self, file_name, trunc=None):
        if trunc is None:
            return deref(self.thisptr).open(file_name)
        else:
            return deref(self.thisptr).open(file_name, trunc)
    def close(self):
        return deref(self.thisptr).close()
    @staticmethod
    def is_valid_file(s):
        return File.is_valid_file(s)

    def file_version(self):
        return deref(self.thisptr).file_version()

    def have_channel_id_params(self):
        return deref(self.thisptr).have_channel_id_params()
    def get_channel_id_params(self):
        return deref(self.thisptr).get_channel_id_params()

    def have_tracking_id_params(self):
        return deref(self.thisptr).have_tracking_id_params()
    def get_tracking_id_params(self):
        return deref(self.thisptr).get_tracking_id_params()

    def have_sequences_params(self):
        return deref(self.thisptr).have_sequences_params()
    def get_sequences_params(self):
        return deref(self.thisptr).get_sequences_params()

    def get_raw_samples_read_name_list(self):
        return deref(self.thisptr).get_raw_samples_read_name_list()
    def have_raw_samples(self, rn=None):
        if rn is None:
            return deref(self.thisptr).have_raw_samples()
        else:
            return deref(self.thisptr).have_raw_samples(rn)
    def get_raw_samples_params(self, rn=None):
        if rn is None:
            return deref(self.thisptr).get_raw_samples_params()
        else:
            return deref(self.thisptr).get_raw_samples_params(rn)
    def get_raw_int_samples(self, rn=None):
        if rn is None:
            return deref(self.thisptr).get_raw_int_samples()
        else:
            return deref(self.thisptr).get_raw_int_samples(rn)
    def get_raw_samples(self, rn=None):
        if rn is None:
            return deref(self.thisptr).get_raw_samples()
        else:
            return deref(self.thisptr).get_raw_samples(rn)

    def get_eventdetection_group_list(self):
        return deref(self.thisptr).get_eventdetection_group_list()
    def have_eventdetection_group(self, gr=None):
        if gr is None:
            return deref(self.thisptr).have_eventdetection_group()
        else:
            return deref(self.thisptr).have_eventdetection_group(gr)
    def get_eventdetection_params(self, gr=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_params()
        else:
            return deref(self.thisptr).get_eventdetection_params(gr)
    def get_eventdetection_read_name_list(self, gr=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_read_name_list()
        else:
            return deref(self.thisptr).get_eventdetection_read_name_list(gr)
    def have_eventdetection_events(self, gr=None, rn=None):
        if gr is None:
            return deref(self.thisptr).have_eventdetection_events()
        elif rn is None:
            return deref(self.thisptr).have_eventdetection_events(gr)
        else:
            return deref(self.thisptr).have_eventdetection_events(gr, rn)
    def get_eventdetection_events_params(self, gr=None, rn=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_events_params()
        elif rn is None:
            return deref(self.thisptr).get_eventdetection_events_params(gr)
        else:
            return deref(self.thisptr).get_eventdetection_events_params(gr, rn)
    def get_eventdetection_events(self, gr=None, rn=None):
        if gr is None:
            return deref(self.thisptr).get_eventdetection_events()
        elif rn is None:
            return deref(self.thisptr).get_eventdetection_events(gr)
        else:
            return deref(self.thisptr).get_eventdetection_events(gr, rn)

    def get_basecall_group_list(self):
        return deref(self.thisptr).get_basecall_group_list()
    def get_basecall_strand_group_list(self, st):
        return deref(self.thisptr).get_basecall_strand_group_list(st)
    def have_basecall_group(self, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_group()
        else:
            return deref(self.thisptr).have_basecall_group(gr)
    def have_basecall_strand_group(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_strand_group(st)
        else:
            return deref(self.thisptr).have_basecall_strand_group(st, gr)
    def get_basecall_1d_group(self, gr):
        return deref(self.thisptr).get_basecall_1d_group(gr)
    def get_basecall_eventdetection_group(self, gr):
        return deref(self.thisptr).get_basecall_eventdetection_group(gr)
    def get_basecall_params(self, gr):
        return deref(self.thisptr).get_basecall_params(gr)
    def get_basecall_log(self, gr):
        return deref(self.thisptr).get_basecall_log(gr)

    def have_basecall_fastq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_fastq(st)
        else:
            return deref(self.thisptr).have_basecall_fastq(st, gr)
    def get_basecall_fastq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_fastq(st)
        else:
            return deref(self.thisptr).get_basecall_fastq(st, gr)
    def have_basecall_seq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_seq(st)
        else:
            return deref(self.thisptr).have_basecall_seq(st, gr)
    def get_basecall_seq(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_seq(st)
        else:
            return deref(self.thisptr).get_basecall_seq(st, gr)

    def have_basecall_model(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_model(st)
        else:
            return deref(self.thisptr).have_basecall_model(st, gr)
    def get_basecall_model_file(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_model_file(st)
        else:
            return deref(self.thisptr).get_basecall_model_file(st, gr)
    def get_basecall_model_params(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_model_params(st)
        else:
            return deref(self.thisptr).get_basecall_model_params(st, gr)
    def get_basecall_model(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_model(st)
        else:
            return deref(self.thisptr).get_basecall_model(st, gr)

    def have_basecall_events(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_events(st)
        else:
            return deref(self.thisptr).have_basecall_events(st, gr)
    def get_basecall_events_params(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_events_params(st)
        else:
            return deref(self.thisptr).get_basecall_events_params(st, gr)
    def get_basecall_events(self, st, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_events(st)
        else:
            return deref(self.thisptr).get_basecall_events(st, gr)

    def have_basecall_alignment(self, gr=None):
        if gr is None:
            return deref(self.thisptr).have_basecall_alignment()
        else:
            return deref(self.thisptr).have_basecall_alignment(gr)
    def get_basecall_alignment(self, gr=None):
        if gr is None:
            return deref(self.thisptr).get_basecall_alignment()
        else:
            return deref(self.thisptr).get_basecall_alignment(gr)
