#ifndef __FAST5_HPP
#define __FAST5_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>

#include "hdf5_tools.hpp"
#include "fast5_pack.hpp"
#define MAX_K_LEN 8

namespace
{
    inline static std::string array_to_string(std::array< char, MAX_K_LEN > const & a)
    {
        return std::string(a.begin(), std::find(a.begin(), a.end(), '\0'));
    }
}

namespace fast5
{

typedef std::map< std::string, std::string > Attr_Map;

struct Channel_Id_Parameters
{
    std::string channel_number;
    double digitisation;
    double offset;
    double range;
    double sampling_rate;
    Channel_Id_Parameters()
        : channel_number(""),
          digitisation(0.0),
          offset(0.0),
          range(0.0),
          sampling_rate(0.0) {}
}; // struct Channel_Id_Parameters

typedef Attr_Map Tracking_Id_Parameters;

typedef Attr_Map Sequences_Parameters;

typedef float Raw_Samples_Entry;
typedef int16_t Raw_Samples_Int_Entry;

struct Raw_Samples_Parameters
{
    std::string read_id;
    long long read_number;
    long long start_mux;
    long long start_time;
    long long duration;
}; // struct Raw_Samples_Parameters

struct Raw_Samples_Pack
{
    fast5_pack::Huffman_Coder::Code_Type signal;
    Attr_Map signal_param;
}; // struct Raw_Samples_Pack

struct EventDetection_Event_Entry
{
    double mean;
    double stdv;
    long long start;
    long long length;
    friend bool operator == (EventDetection_Event_Entry const & lhs, EventDetection_Event_Entry const & rhs)
    {
        return lhs.mean == rhs.mean
            and lhs.stdv == rhs.stdv
            and lhs.start == rhs.start
            and lhs.length == rhs.length;
    }
}; // struct EventDetection_Event

struct EventDetection_Event_Parameters
{
    std::string read_id;
    long long read_number;
    long long scaling_used;
    long long start_mux;
    long long start_time;
    long long duration;
    double median_before;
    unsigned abasic_found;
}; // struct EventDetection_Event_Parameters

struct EventDetection_Events_Pack
{
    fast5_pack::Huffman_Coder::Code_Type skip;
    Attr_Map skip_param;
    fast5_pack::Huffman_Coder::Code_Type len;
    Attr_Map len_param;
}; // struct EventDetection_Events_Pack

//
// This struct represents the expected signal measured
// given the kmer sequence that is in the pore when the
// the observations are made. A pore model consists
// of 1024 of these entries (one per 5-mer) and global
// shift/scaling parameters.
//
struct Model_Entry
{
    long long variant;
    double level_mean;
    double level_stdv;
    double sd_mean;
    double sd_stdv;
    double weight;
    std::array< char, MAX_K_LEN > kmer;
    std::string get_kmer() const { return array_to_string(kmer); }
    friend bool operator == (Model_Entry const & lhs, Model_Entry const & rhs)
    {
        return lhs.variant == rhs.variant
            and lhs.level_mean == rhs.level_mean
            and lhs.level_stdv == rhs.level_stdv
            and lhs.sd_mean == rhs.sd_mean
            and lhs.sd_stdv == rhs.sd_stdv
            and lhs.weight == rhs.weight
            and lhs.kmer == rhs.kmer;
    }
}; // struct Model_Entry

//
// This struct represents the global transformations
// that must be applied to each Model_Entry
//
struct Model_Parameters
{
    double scale;
    double shift;
    double drift;
    double var;
    double scale_sd;
    double var_sd;
}; // struct Model_Parameters


struct Basecall_Fastq_Pack
{
    fast5_pack::Huffman_Coder::Code_Type bp;
    Attr_Map bp_param;
    fast5_pack::Huffman_Coder::Code_Type qv;
    Attr_Map qv_param;
    std::string read_name;
    std::uint8_t qv_bits;
}; // struct Basecall_Fastq_Pack

//
// This struct represents an observed event.
// The members of the struct are the same as
// the fields encoded in the FAST5 file.
//
struct Event_Entry
{
    double mean;
    double stdv;
    double start;
    double length;
    double p_model_state;
    double p_mp_state;
    double p_A;
    double p_C;
    double p_G;
    double p_T;
    long long move;
    std::array< char, MAX_K_LEN > model_state;
    std::array< char, MAX_K_LEN > mp_state;
    std::string get_model_state() const { return array_to_string(model_state); }
    std::string get_mp_state() const { return array_to_string(mp_state); }
    friend bool operator == (Event_Entry const & lhs, Event_Entry const & rhs)
    {
        return lhs.mean == rhs.mean
            and lhs.stdv == rhs.stdv
            and lhs.start == rhs.start
            and lhs.length == rhs.length
            and lhs.p_model_state == rhs.p_model_state
            and lhs.p_mp_state == rhs.p_mp_state
            and lhs.p_A == rhs.p_A
            and lhs.p_C == rhs.p_C
            and lhs.p_G == rhs.p_G
            and lhs.p_T == rhs.p_T
            and lhs.move == rhs.move
            and lhs.model_state == rhs.model_state
            and lhs.mp_state == rhs.mp_state;
    }
}; // struct Event_Entry

struct Basecall_Event_Parameters
{
    double start_time;
    double duration;
};

struct Basecall_Events_Pack
{
    fast5_pack::Huffman_Coder::Code_Type skip;
    Attr_Map skip_param;
    fast5_pack::Huffman_Coder::Code_Type move;
    Attr_Map move_param;
    fast5_pack::Huffman_Coder::Code_Type bases;
    Attr_Map bases_param;
    //
    std::string ed_gr;
    long long start_time;
    long long duration;
    unsigned state_size;
}; // struct Basecall_Events_Pack

//
// This struct represents a template-to-complement
// match that is emitted by ONT's 2D basecaller
//
struct Event_Alignment_Entry
{
    long long template_index;
    long long complement_index;
    std::array< char, MAX_K_LEN > kmer;
    std::string get_kmer() const { return array_to_string(kmer); }
    friend bool operator == (Event_Alignment_Entry const & lhs, Event_Alignment_Entry const & rhs)
    {
        return lhs.template_index == rhs.template_index
            and lhs.complement_index == rhs.complement_index
            and lhs.kmer == rhs.kmer;
    }
}; // struct Event_Alignment_Entry


class File
    : private hdf5_tools::File
{
private:
    typedef hdf5_tools::File Base;
public:
    //using Base::is_open;
    //using Base::is_rw;
    //using Base::file_name;
    //using Base::create;
    //using Base::close;
    using Base::get_object_count;
    using Base::is_valid_file;
    //using Base::write;

    File() = default;
    File(std::string const & file_name, bool rw = false) { open(file_name, rw); }

    bool is_open() const { return static_cast< const Base* >(this)->is_open(); }
    bool is_rw() const { return static_cast< const Base* >(this)->is_rw(); }
    std::string const & file_name() const { return static_cast< const Base* >(this)->file_name(); }
    void create(std::string const & file_name, bool truncate = false) { static_cast< Base* >(this)->create(file_name, truncate); }
    void close() { static_cast< Base* >(this)->close(); }

    void open(std::string const & file_name, bool rw = false)
    {
        Base::open(file_name, rw);
        if (is_open())
        {
            // get channel_id_params
            load_channel_id_params();
            // detect raw samples read name
            detect_raw_samples_read_name_list();
            // detect eventdetection groups
            detect_eventdetection_group_list();
            // detect basecall groups
            detect_basecall_group_list();
        }
    }

    /**
     * Extract "/file_version" attribute. This must exist.
     */
    std::string file_version() const
    {
        std::string res;
        assert(Base::exists(file_version_path()));
        Base::read(file_version_path(), res);
        return res;
    }

    /**
     * Check if "/UniqueGlobalKey/channel_id" attributes exist.
     */
    bool have_channel_id_params() const
    {
        return _channel_id_params.sampling_rate > 0.0;
    }
    /**
     * Extract "/UniqueGlobalKey/channel_id" attributes.
     */
    Channel_Id_Parameters get_channel_id_params() const
    {
        return _channel_id_params;
    }
    /**
     * Check if sampling rate exists.
     */
    bool have_sampling_rate() const
    {
        return have_channel_id_params();
    }
    /**
     * Get sampling rate.
     */
    double get_sampling_rate() const
    {
        return _channel_id_params.sampling_rate;
    }

    /**
     * Check if "/UniqueGlobalKey/tracking_id" attributes exist.
     */
    bool have_tracking_id_params() const
    {
        return Base::group_exists(tracking_id_path());
    }
    /**
     * Extract "/UniqueGlobalKey/tracking_id" attributes.
     */
    Tracking_Id_Parameters get_tracking_id_params() const
    {
        return get_attr_map(tracking_id_path());
    }

    /**
     * Check if sequences attributes exists.
     */
    bool have_sequences_params() const
    {
        return Base::group_exists(sequences_path());
    }
    /**
     * Get sequences attributes.
     */
    Sequences_Parameters get_sequences_params() const
    {
        return get_attr_map(sequences_path());
    }

    /**
     * Get list of raw samples read names.
     */
    std::vector< std::string > const & get_raw_samples_read_name_list() const
    {
        return _raw_samples_read_name_list;
    }
    /**
     * Check if raw samples exist.
     * If _rn non-empty, check if raw samples exist for given read.
     */
    bool have_raw_samples(std::string const & _rn = std::string()) const
    {
        if (not have_channel_id_params())
        {
            return false;
        }
        auto rn_l = get_raw_samples_read_name_list();
        if (_rn.empty())
        {
            return not rn_l.empty();
        }
        else
        {
            std::set< std::string > rn_d(rn_l.begin(), rn_l.end());
            return rn_d.count(_rn) > 0;
        }
    }
    bool have_raw_samples_unpack(std::string const & rn) const
    {
        return Base::dataset_exists(raw_samples_path(rn));
    }
    bool have_raw_samples_pack(std::string const & rn) const
    {
        return Base::group_exists(raw_samples_pack_path(rn));
    }
    /**
     * Get raw samples attributes for given read name (default: first read name).
     */
    Raw_Samples_Parameters get_raw_samples_params(std::string const & _rn = std::string()) const
    {
        Raw_Samples_Parameters res;
        std::string const rn = not _rn.empty()? _rn : get_raw_samples_read_name_list().front();
        std::string p = raw_samples_params_path(rn);
        Base::read(p + "/read_id", res.read_id);
        Base::read(p + "/read_number", res.read_number);
        Base::read(p + "/start_mux", res.start_mux);
        Base::read(p + "/start_time", res.start_time);
        Base::read(p + "/duration", res.duration);
        return res;
    }
    void add_raw_samples_params(std::string const & rn, Raw_Samples_Parameters const & params) const
    {
        std::string p = raw_samples_params_path(rn);
        Base::write_attribute(p + "/read_id", params.read_id);
        Base::write_attribute(p + "/read_number", params.read_number);
        Base::write_attribute(p + "/start_mux", params.start_mux);
        Base::write_attribute(p + "/start_time", params.start_time);
        Base::write_attribute(p + "/duration", params.duration);
    }
    /**
     * Get raw samples for given read name as ints (default: first read name).
     */
    std::vector< Raw_Samples_Int_Entry > get_raw_samples_int(std::string const & _rn = std::string()) const
    {
        // get raw samples
        std::vector< Raw_Samples_Int_Entry > res;
        std::string const rn = not _rn.empty()? _rn : get_raw_samples_read_name_list().front();
        if (Base::dataset_exists(raw_samples_path(rn)))
        {
            Base::read(raw_samples_path(rn), res);
        }
        else
        {
            auto rsp = get_raw_samples_pack(rn);
            res = unpack_rw(rsp);
        }
        return res;
    }
    /**
     * Get raw samples for given read name (default: first read name).
     */
    std::vector< Raw_Samples_Entry > get_raw_samples(std::string const & _rn = std::string()) const
    {
        // get raw samples
        auto raw_samples_int = get_raw_samples_int(_rn);
        // decode levels
        std::vector< Raw_Samples_Entry > res;
        res.reserve(raw_samples_int.size());
        for (auto int_level : raw_samples_int)
        {
            res.push_back(raw_sample_to_float(int_level));
        }
        return res;
    }
    Raw_Samples_Pack get_raw_samples_pack(std::string const & rn) const
    {
        Raw_Samples_Pack rsp;
        Base::read(raw_samples_pack_path(rn) + "/Signal", rsp.signal);
        rsp.signal_param = get_attr_map(raw_samples_pack_path(rn) + "/Signal");
        return rsp;
    }
    /**
     * Add raw samples.
     */
    void add_raw_samples_int(std::string const & rn, std::vector< Raw_Samples_Int_Entry > const & rs) const
    {
        Base::write_dataset(raw_samples_path(rn), rs);
    }
    /**
     * Add packed raw smaples.
     */
    void add_raw_samples_pack(std::string const & rn, Raw_Samples_Pack const & rsp) const
    {
        Base::write_dataset(raw_samples_pack_path(rn) + "/Signal", rsp.signal);
        add_attr_map(raw_samples_pack_path(rn) + "/Signal", rsp.signal_param);
    }

    /**
     * Get list of EventDetection groups.
     */
    std::vector< std::string > const & get_eventdetection_group_list() const
    {
        return _eventdetection_group_list;
    }
    /**
     * Check if any EventDetection groups exist.
     */
    bool have_eventdetection_groups() const
    {
        return not get_eventdetection_group_list().empty();
    }
    /**
     * Get list of reads for given EventDetection group (default: first EventDetection group).
     */
    std::vector< std::string > get_eventdetection_read_name_list(std::string const & _gr = std::string()) const
    {
        std::string const & gr = not _gr.empty()? _gr : get_eventdetection_group_list().front();
        return detect_eventdetection_read_name_list(gr);
    }
    /**
     * Check if EventDetection events exist.
     * If _gr given: check if events exist for given group; else: check first EventDetection group.
     * If _rn given: check if events exist for given group and read name.
     */
    bool have_eventdetection_events(
        std::string const & _gr = std::string(),
        std::string const & _rn = std::string()) const
    {
        std::string gr;
        if (_gr.empty())
        {
            auto gr_l = get_eventdetection_group_list();
            if (gr_l.empty()) return false;
            gr = gr_l.front();
        }
        else
        {
            gr = _gr;
        }
        auto rn_l = get_eventdetection_read_name_list(gr);
        if (_rn.empty())
        {
            return not rn_l.empty();
        }
        else
        {
            std::set< std::string > rn_d(rn_l.begin(), rn_l.end());
            return rn_d.count(_rn) > 0;
        }
    }
    bool have_eventdetection_events_unpack(std::string const & gr, std::string const & rn) const
    {
        return Base::dataset_exists(eventdetection_events_path(gr, rn));
    }
    bool have_eventdetection_events_pack(std::string const & gr, std::string const & rn) const
    {
        return Base::group_exists(eventdetection_events_pack_path(gr, rn));
    }
    /**
     * Get EventDetection params for given EventDetection group (default: first EventDetection group).
     */
    Attr_Map get_eventdetection_params(std::string const & _gr = std::string()) const
    {
        std::string const & gr = not _gr.empty()? _gr : get_eventdetection_group_list().front();
        return get_attr_map(eventdetection_params_path(gr));
    }
    void add_eventdetection_params(std::string const & gr, Attr_Map const & am) const
    {
        add_attr_map(eventdetection_params_path(gr), am);
    }
    /**
     * Get EventDetection event params for given EventDetection group, and given read name
     * (default: first EventDetection group, and first read name in it).
     */
    EventDetection_Event_Parameters get_eventdetection_event_params(
        std::string const & _gr = std::string(), std::string const & _rn = std::string()) const
    {
        EventDetection_Event_Parameters res;
        std::string const & gr = not _gr.empty()? _gr : get_eventdetection_group_list().front();
        std::string const rn = not _rn.empty()? _rn : get_eventdetection_read_name_list(gr).front();
        auto p = eventdetection_event_params_path(gr, rn);
        auto a_v = Base::get_attr_list(p);
        std::set< std::string > a_s(a_v.begin(), a_v.end());
        Base::read(p + "/read_number", res.read_number);
        Base::read(p + "/scaling_used", res.scaling_used);
        Base::read(p + "/start_mux", res.start_mux);
        Base::read(p + "/start_time", res.start_time);
        Base::read(p + "/duration", res.duration);
        // optional fields
        if (a_s.count("read_id"))
        {
            Base::read(p + "/read_id", res.read_id);
        }
        if (a_s.count("median_before"))
        {
            Base::read(p + "/median_before", res.median_before);
        }
        else
        {
            res.median_before = -1;
        }
        if (a_s.count("abasic_found"))
        {
            Base::read(p + "/abasic_found", res.abasic_found);
        }
        else
        {
            res.abasic_found = 2;
        }
        return res;
    }
    void add_eventdetection_event_params(std::string const & gr, std::string const & rn,
                                         EventDetection_Event_Parameters const & ede_params) const
    {
        auto p = eventdetection_event_params_path(gr, rn);
        if (not ede_params.read_id.empty()) Base::write_attribute(p + "/read_id", ede_params.read_id);
        Base::write_attribute(p + "/read_number", ede_params.read_number);
        Base::write_attribute(p + "/scaling_used", ede_params.scaling_used);
        Base::write_attribute(p + "/start_mux", ede_params.start_mux);
        Base::write_attribute(p + "/start_time", ede_params.start_time);
        Base::write_attribute(p + "/duration", ede_params.duration);
        if (ede_params.median_before > 0.0) Base::write_attribute(p + "/median_before", ede_params.median_before);
        if (ede_params.abasic_found < 2) Base::write_attribute(p + "/abasic_found", ede_params.abasic_found);
    }
    /**
     * Get EventDetection events for given EventDetection group, and given read name.
     */
    std::vector< EventDetection_Event_Entry > get_eventdetection_events(
        std::string const & _gr = std::string(), std::string const & _rn = std::string()) const
    {
        std::vector< EventDetection_Event_Entry > res;
        std::string const & gr = not _gr.empty()? _gr : get_eventdetection_group_list().front();
        std::string const rn = not _rn.empty()? _rn : get_eventdetection_read_name_list(gr).front();
        if (have_eventdetection_events_unpack(gr, rn))
        {
            auto p = eventdetection_events_path(gr, rn);
            auto struct_member_names = Base::get_struct_members(p);
            assert(struct_member_names.size() >= 4);
            bool have_stdv = false;
            bool have_variance = false;
            for (auto const & s : struct_member_names)
            {
                if (s == "stdv") have_stdv = true;
                else if (s == "variance") have_variance = true;
            }
            hdf5_tools::Compound_Map m;
            m.add_member("mean", &EventDetection_Event_Entry::mean);
            m.add_member("start", &EventDetection_Event_Entry::start);
            m.add_member("length", &EventDetection_Event_Entry::length);
            if (have_stdv)
            {
                m.add_member("stdv", &EventDetection_Event_Entry::stdv);
            }
            else if (have_variance)
            {
                m.add_member("variance", &EventDetection_Event_Entry::stdv);
            }
            else
            {
                // must have stdv or variance
                abort();
            }
            Base::read(p, res, m);
            if (not have_stdv)
            {
                // have read variances
                for (auto& e : res)
                {
                    e.stdv = std::sqrt(e.stdv);
                }
            }
        } // have ed unpack
        else // have pack
        {
            auto ed_pack = get_eventdetection_events_pack(gr, rn);
            auto ed_param = get_eventdetection_event_params(gr, rn);
            auto rs = get_raw_samples(rn);
            auto rs_param = get_raw_samples_params(rn);
            res = unpack_ed(ed_pack, ed_param, rs, rs_param);
        }
        return res;
    } // get_eventdetection_events()
    void add_eventdetection_events(
        std::string const & gr, std::string const & rn,
        std::vector< EventDetection_Event_Entry > const & ed) const
    {
        hdf5_tools::Compound_Map m;
        m.add_member("mean", &EventDetection_Event_Entry::mean);
        m.add_member("start", &EventDetection_Event_Entry::start);
        m.add_member("length", &EventDetection_Event_Entry::length);
        m.add_member("stdv", &EventDetection_Event_Entry::stdv);
        Base::write_dataset(eventdetection_events_path(gr, rn), ed, m);
    }
    EventDetection_Events_Pack get_eventdetection_events_pack(
        std::string const & gr, std::string const & rn) const
    {
        EventDetection_Events_Pack edp;
        Base::read(eventdetection_events_pack_path(gr, rn) + "/Skip", edp.skip);
        edp.skip_param = get_attr_map(eventdetection_events_pack_path(gr, rn) + "/Skip");
        Base::read(eventdetection_events_pack_path(gr, rn) + "/Len", edp.len);
        edp.len_param = get_attr_map(eventdetection_events_pack_path(gr, rn) + "/Len");
        return edp;
    }
    void add_eventdetection_events_pack(
        std::string const & gr, std::string const & rn,
        EventDetection_Events_Pack const & edp) const
    {
        Base::write_dataset(eventdetection_events_pack_path(gr, rn) + "/Skip", edp.skip);
        add_attr_map(eventdetection_events_pack_path(gr, rn) + "/Skip", edp.skip_param);
        Base::write_dataset(eventdetection_events_pack_path(gr, rn) + "/Len", edp.len);
        add_attr_map(eventdetection_events_pack_path(gr, rn) + "/Len", edp.len_param);
    }

    /**
     * Get list of all Basecall groups.
     */
    std::vector< std::string > const & get_basecall_group_list() const
    {
        return _basecall_group_list;
    }
    /**
     * Check if any Basecall groups exist.
     */
    bool have_basecall_groups() const
    {
        return not get_basecall_group_list().empty();
    }
    /**
     * Get list of Basecall groups for given strand.
     */
    std::vector< std::string > const & get_basecall_strand_group_list(unsigned st) const
    {
        return _basecall_strand_group_list[st];
    }
    /**
     * Check if any Basecall groups exist for given strand.
     */
    bool have_basecall_strand_groups(unsigned st) const
    {
        return not get_basecall_strand_group_list(st).empty();
    }
    /**
     * Get Basecall group params for given Basecall group.
     */
    Attr_Map get_basecall_params(std::string const & gr) const
    {
        return get_attr_map(basecall_root_path() + "/" + basecall_group_prefix() + gr);
    }
    void add_basecall_params(std::string const & gr, Attr_Map const & am) const
    {
        add_attr_map(basecall_root_path() + "/" + basecall_group_prefix() + gr, am);
    }
    /**
     * Check if Basecall log exists for given Basecall group.
     */
    bool have_basecall_log(std::string const & gr) const
    {
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + gr + "/Log";
        return Base::exists(path);
    }
    /**
     * Get Basecall log for given Basecall group.
     */
    std::string get_basecall_log(std::string const & gr) const
    {
        std::string res;
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + gr + "/Log";
        Base::read(path, res);
        return res;
    }
    /**
     * Check if Basecall fastq exists for given Basecall group and given strand.
     */
    bool have_basecall_fastq(unsigned st, std::string const & _gr = std::string()) const
    {
        if (_gr.empty() and get_basecall_strand_group_list(st).empty()) return false;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        return have_basecall_fastq_unpack(st, gr) or have_basecall_fastq_pack(st, gr);
    }
    bool have_basecall_fastq_unpack(unsigned st, std::string const & gr) const
    {
        return Base::dataset_exists(basecall_fastq_path(gr, st));
    }
    bool have_basecall_fastq_pack(unsigned st, std::string const & gr) const
    {
        return Base::group_exists(basecall_fastq_pack_path(gr, st));
    }
    /**
     * Get Basecall fastq for given Basecall group and given strand.
     */
    std::string get_basecall_fastq(unsigned st, std::string const & _gr = std::string()) const
    {
        std::string res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        if (have_basecall_fastq_unpack(st, gr))
        {
            Base::read(basecall_fastq_path(gr, st), res);
        }
        else
        {
            auto fq_pack = get_basecall_fastq_pack(st, gr);
            res = unpack_fq(fq_pack);
        }
        return res;
    }
    /**
     * Add Basecall fastq
     */
    void add_basecall_fastq(unsigned st, std::string const & gr, std::string const & fq) const
    {
        Base::write(basecall_fastq_path(gr, st), true, fq);
    }
    Basecall_Fastq_Pack get_basecall_fastq_pack(unsigned st, std::string const & gr) const
    {
        Basecall_Fastq_Pack fq_pack;
        auto p = basecall_fastq_pack_path(gr, st);
        Base::read(p + "/BP", fq_pack.bp);
        fq_pack.bp_param = get_attr_map(p + "/BP");
        Base::read(p + "/QV", fq_pack.qv);
        fq_pack.qv_param = get_attr_map(p + "/QV");
        Base::read(p + "/read_name", fq_pack.read_name);
        Base::read(p + "/qv_bits", fq_pack.qv_bits);
        return fq_pack;        
    }
    void add_basecall_fastq_pack(unsigned st, std::string const & gr, Basecall_Fastq_Pack const & fq_pack) const
    {
        auto p = basecall_fastq_pack_path(gr, st);
        Base::write_dataset(p + "/BP", fq_pack.bp);
        add_attr_map(p + "/BP", fq_pack.bp_param);
        Base::write_dataset(p + "/QV", fq_pack.qv);
        add_attr_map(p + "/QV", fq_pack.qv_param);
        Base::write_attribute(p + "/read_name", fq_pack.read_name);
        Base::write_attribute(p + "/qv_bits", fq_pack.qv_bits);
    }
    /**
     * Check if Basecall seq exists for given Basecall group and given strand.
     */
    bool have_basecall_seq(unsigned st, std::string const & _gr = std::string()) const
    {
        return have_basecall_fastq(st, _gr);
    }
    /**
     * Get Basecall sequence for given Basecall group and given strand.
     */
    std::string get_basecall_seq(unsigned st, std::string const & _gr = std::string()) const
    {
        return fq2seq(get_basecall_fastq(st, _gr));
    }
    /**
     * Add Basecall seq
     */
    void add_basecall_seq(unsigned st, std::string const & gr,
                          std::string const & name, std::string const & seq, int default_qual = 33) const
    {
        std::ostringstream oss;
        oss << '@' << name << std::endl
            << seq << std::endl
            << '+' << std::endl
            << std::string(seq.size(), static_cast< char >(default_qual));
        add_basecall_fastq(st, gr, oss.str());
    }
    /**
     * Check if Basecall model exist for given Basecall group and given strand.
     */
    bool have_basecall_model(unsigned st, std::string const & _gr = std::string()) const
    {
        if (_gr.empty() and get_basecall_strand_group_list(st).empty()) return false;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        auto gr_1d = get_basecall_group_1d(gr);
        return Base::dataset_exists(basecall_model_path(gr_1d, st));
    }
    /**
     * Get Basecall model file name for given Basecall group and given strand.
     */
    std::string get_basecall_model_file(unsigned st, std::string const & _gr = std::string()) const
    {
        std::string res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        auto gr_1d = get_basecall_group_1d(gr);
        assert(Base::exists(basecall_model_file_path(gr_1d, st)));
        Base::read(basecall_model_file_path(gr_1d, st), res);
        return res;
    }
    void add_basecall_model_file(unsigned st, std::string const & gr, std::string const & file_name) const
    {
        auto gr_1d = get_basecall_group_1d(gr);
        std::string path = basecall_model_file_path(gr_1d, st);
        Base::write(path, false, file_name);
    }
    /**
     * Get Basecall model parameters for given Basecall group and given strand.
     */
    Model_Parameters get_basecall_model_params(unsigned st, std::string const & _gr = std::string()) const
    {
        Model_Parameters res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        auto gr_1d = get_basecall_group_1d(gr);
        std::string path = basecall_model_path(gr_1d, st);
        Base::read(path + "/scale", res.scale);
        Base::read(path + "/shift", res.shift);
        Base::read(path + "/drift", res.drift);
        Base::read(path + "/var", res.var);
        Base::read(path + "/scale_sd", res.scale_sd);
        Base::read(path + "/var_sd", res.var_sd);
        return res;
    }
    template < typename T >
    void add_basecall_model_params(unsigned st, std::string const & gr, T const & params) const
    {
        auto gr_1d = get_basecall_group_1d(gr);
        std::string path = basecall_model_path(gr_1d, st);
        Base::write(path + "/scale", false, params.scale);
        Base::write(path + "/shift", false, params.shift);
        Base::write(path + "/drift", false, params.drift);
        Base::write(path + "/var", false, params.var);
        Base::write(path + "/scale_sd", false, params.scale_sd);
        Base::write(path + "/var_sd", false, params.var_sd);
    }
    /**
     * Get Basecall model for given Basecall group and given strand.
     */
    std::vector< Model_Entry > get_basecall_model(unsigned st, std::string const & _gr = std::string()) const
    {
        std::vector< Model_Entry > res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        hdf5_tools::Compound_Map m;
        m.add_member("kmer", &Model_Entry::kmer);
        m.add_member("level_mean", &Model_Entry::level_mean);
        m.add_member("level_stdv", &Model_Entry::level_stdv);
        m.add_member("sd_mean", &Model_Entry::sd_mean);
        m.add_member("sd_stdv", &Model_Entry::sd_stdv);
        auto gr_1d = get_basecall_group_1d(gr);
        Base::read(basecall_model_path(gr_1d, st), res, m);
        return res;
    }
    /**
     * Add Basecall model
     */
    template < typename T >
    void add_basecall_model(unsigned st, std::string const & gr, std::vector< T > const & m) const
    {
        hdf5_tools::Compound_Map cm;
        cm.add_member("kmer", &T::kmer);
        cm.add_member("level_mean", &T::level_mean);
        cm.add_member("level_stdv", &T::level_stdv);
        cm.add_member("sd_mean", &T::sd_mean);
        cm.add_member("sd_stdv", &T::sd_stdv);
        auto gr_1d = get_basecall_group_1d(gr);
        Base::write(basecall_model_path(gr_1d, st), true, m, cm);
    }
    /**
     * Check if Basecall events exist for given Basecall group and given strand.
     */
    bool have_basecall_events(unsigned st, std::string const & _gr = std::string()) const
    {
        if (_gr.empty() and get_basecall_strand_group_list(st).empty()) return false;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        auto gr_1d = get_basecall_group_1d(gr);
        return (Base::dataset_exists(basecall_events_path(gr_1d, st))
                or Base::group_exists(basecall_events_pack_path(gr_1d, st)));
    }
    bool have_basecall_events_unpack(unsigned st, std::string const & gr) const
    {
        return Base::dataset_exists(basecall_events_path(get_basecall_group_1d(gr), st));
    }
    bool have_basecall_events_pack(unsigned st, std::string const & gr) const
    {
        return Base::group_exists(basecall_events_pack_path(get_basecall_group_1d(gr), st));
    }
    /**
     * Get Basecall events for given Basecall group and given strand.
     */
    std::vector< Event_Entry > get_basecall_events(unsigned st, std::string const & _gr = std::string()) const
    {
        std::vector< Event_Entry > res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        auto gr_1d = get_basecall_group_1d(gr);
        if (have_basecall_events_unpack(st, gr))
        {
            hdf5_tools::Compound_Map m;
            m.add_member("mean", &Event_Entry::mean);
            m.add_member("start", &Event_Entry::start);
            m.add_member("stdv", &Event_Entry::stdv);
            m.add_member("length", &Event_Entry::length);
            m.add_member("p_model_state", &Event_Entry::p_model_state);
            m.add_member("model_state", &Event_Entry::model_state);
            m.add_member("move", &Event_Entry::move);
            Base::read(basecall_events_path(gr_1d, st), res, m);
        }
        else // have pack
        {
            auto ev_pack = get_basecall_events_pack(st, gr);
            if (not have_basecall_fastq(st, gr))
            {
                throw std::invalid_argument("missing fastq for basecall events unpacking");
            }
            auto fq = get_basecall_fastq(st, gr);
            if (not have_eventdetection_events(ev_pack.ed_gr))
            {
                throw std::invalid_argument("missing evendetection events for basecall events unpacking");
            }
            auto ed = get_eventdetection_events(ev_pack.ed_gr);
            res = unpack_ev(ev_pack, fq, ed, _channel_id_params.sampling_rate);
        }
        return res;
    }
    /**
     * Add Basecall events
     */
    template < typename T >
    void add_basecall_events(unsigned st, std::string const & gr, std::vector< T > const & ev) const
    {
        hdf5_tools::Compound_Map cm;
        cm.add_member("mean", &T::mean);
        cm.add_member("start", &T::start);
        cm.add_member("stdv", &T::stdv);
        cm.add_member("length", &T::length);
        cm.add_member("p_model_state", &T::p_model_state);
        cm.add_member("model_state", &T::model_state);
        cm.add_member("move", &T::move);
        auto gr_1d = get_basecall_group_1d(gr);
        Base::write(basecall_events_path(gr_1d, st), true, ev, cm);
    }
    Basecall_Events_Pack get_basecall_events_pack(unsigned st, std::string const & gr) const
    {
        auto p = basecall_events_pack_path(gr, st);
        Basecall_Events_Pack ev_pack;
        Base::read(p + "/Skip", ev_pack.skip);
        ev_pack.skip_param = get_attr_map(p + "/Skip");
        Base::read(p + "/Move", ev_pack.move);
        ev_pack.move_param = get_attr_map(p + "/Move");
        //Base::read(p + "/Bases", ev_pack.bases);
        //ev_pack.bases_param = get_attr_map(p + "/Bases");
        Base::read(p + "/ed_gr", ev_pack.ed_gr);
        Base::read(p + "/start_time", ev_pack.start_time);
        Base::read(p + "/duration", ev_pack.duration);
        Base::read(p + "/state_size", ev_pack.state_size);
        return ev_pack;
    }
    void add_basecall_events_pack(unsigned st, std::string const & gr, Basecall_Events_Pack const & ev_pack) const
    {
        auto p = basecall_events_pack_path(gr, st);
        Base::write_dataset(p + "/Skip", ev_pack.skip);
        add_attr_map(p + "/Skip", ev_pack.skip_param);
        Base::write_dataset(p + "/Move", ev_pack.move);
        add_attr_map(p + "/Move", ev_pack.move_param);
        //Base::write_dataset(p + "/Bases", ev_pack.bases);
        //add_attr_map(p + "/Bases", ev_pack.bases_param);
        Base::write_attribute(p + "/ed_gr", ev_pack.ed_gr);
        Base::write_attribute(p + "/start_time", ev_pack.start_time);
        Base::write_attribute(p + "/duration", ev_pack.duration);
        Base::write_attribute(p + "/state_size", ev_pack.state_size);
    }
    Basecall_Event_Parameters get_basecall_event_params(unsigned st, std::string const & _gr = std::string()) const
    {
        Basecall_Event_Parameters res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(st).front();
        auto gr_1d = get_basecall_group_1d(gr);
        if (have_basecall_events_unpack(st, gr_1d))
        {
            auto p = basecall_events_path(gr_1d, st);
            auto params = get_attr_map(p);
            if (params.count("start_time")) std::istringstream(params["start_time"]) >> res.start_time;
            else res.start_time = 0.0;
            if (params.count("duration")) std::istringstream(params["duration"]) >> res.duration;
            else res.duration = 0.0;
        }
        else if (have_basecall_events_pack(st, gr_1d))
        {
            auto p = basecall_events_pack_path(gr_1d, st);
            Basecall_Events_Pack ev_pack;
            Base::read(p + "/start_time", ev_pack.start_time);
            Base::read(p + "/duration", ev_pack.duration);
            res.start_time = time_to_float(ev_pack.start_time, _channel_id_params.sampling_rate);
            res.duration = time_to_float(ev_pack.duration, _channel_id_params.sampling_rate);
        }
        return res;
    }
    void add_basecall_event_params(unsigned st, std::string const & gr,
                                   Basecall_Event_Parameters const & bce_param) const
    {
        auto gr_1d = get_basecall_group_1d(gr);
        auto p = basecall_events_path(gr_1d, st);
        if (not Base::dataset_exists(p))
        {
            throw std::invalid_argument("basecall events must be added before their params");
        }
        Base::write_attribute(p + "/start_time", bce_param.start_time);
        Base::write_attribute(p + "/duration", bce_param.duration);
    }
    /**
     * Check if Basecall event alignment exist for given Basecall group.
     */
    bool have_basecall_event_alignment(std::string const & _gr = std::string()) const
    {
        if (_gr.empty() and get_basecall_strand_group_list(2).empty()) return false;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(2).front();
        return Base::dataset_exists(basecall_event_alignment_path(gr));
    }
    /**
     * Get Basecall events for given Basecall group.
     */
    std::vector< Event_Alignment_Entry > get_basecall_event_alignment(std::string const & _gr = std::string()) const
    {
        std::vector< Event_Alignment_Entry > res;
        std::string const & gr = not _gr.empty()? _gr : get_basecall_strand_group_list(2).front();
        hdf5_tools::Compound_Map m;
        m.add_member("template", &Event_Alignment_Entry::template_index);
        m.add_member("complement", &Event_Alignment_Entry::complement_index);
        m.add_member("kmer", &Event_Alignment_Entry::kmer);
        Base::read(basecall_event_alignment_path(gr), res, m);
        return res;
    }

    /**
     * Get basecall group holding 1d calls.
     */
    std::string get_basecall_group_1d(std::string const & gr) const
    {
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + gr + "/basecall_1d";
        if (Base::attribute_exists(path))
        {
            std::string tmp;
            Base::read(path, tmp);
            auto tmp1 = tmp.substr(0, 18);
            auto tmp2 = tmp.substr(18);
            if (tmp1 == "Analyses/Basecall_"
                and Base::group_exists(basecall_root_path() + "/" + basecall_group_prefix() + tmp2))
            {
                return tmp2;
            }
        }
        return gr;
    }
    /**
     * Get EventDetection group for given Basecall group, if available.
     */
    std::string get_basecall_eventdetection_group(std::string const & gr) const
    {
        if (not have_eventdetection_events()) return "";
        auto bc_param = get_basecall_params(gr);
        if (bc_param.count("event_detection"))
        {
            std::string tmp = bc_param["event_detection"];
            auto pos = tmp.find(eventdetection_group_prefix());
            if (pos != std::string::npos)
            {
                pos += eventdetection_group_prefix().size();
                auto end_pos = tmp.find("/", pos);
                if (end_pos == std::string::npos)
                {
                    end_pos = tmp.size();
                }
                auto ed_gr = tmp.substr(pos, end_pos - pos);
                return (have_eventdetection_events(ed_gr)? ed_gr : "");
            }
        }
        // the hard way
        if (not have_basecall_events(0, gr)) return "";
        auto ev = get_basecall_events(0, gr);
        EventDetection_Event_Entry ede;
        ede.mean = 0.0;
        ede.stdv = 1.0;
        ede.start = time_to_int(ev[0].start, get_sampling_rate());
        ede.length = time_to_int(ev[0].length, get_sampling_rate());
        auto ed_gr_l = get_eventdetection_group_list();
        for (auto & ed_gr : ed_gr_l)
        {
            auto ed = get_eventdetection_events(ed_gr);
            auto p = std::equal_range(
                ed.begin(), ed.end(), ede,
                [] (EventDetection_Event_Entry const & lhs, EventDetection_Event_Entry const & rhs) {
                    return lhs.start < rhs.start;
                });
            if (p.first != p.second and p.first->length == ede.length)
            {
                return ed_gr;
            }
        }
        return "";
    }

    static std::string fq2seq(std::string const & fq)
    {
        return split_fq(fq)[1];
    }
    static std::array< std::string, 4 > split_fq(std::string const & fq)
    {
        std::array< std::string, 4 > res = {{"", "", "", ""}};
        size_t i = 0;
        for (int k = 0; k < 4; ++k)
        {
            if (k % 2 == 0) ++i;
            size_t j = fq.find_first_of('\n', i);
            if (j == std::string::npos)
            {
                if (k == 3)
                {
                    j = fq.size();
                }
                else
                {
                    return {{"", "", "", ""}};
                }
            }
            res[k] = fq.substr(i, j - i);
            i = j + 1;
        }
        return res;
    }

    /**
     * Copy all attributes from one file to another.
     */
    static void copy_attributes(File const & src_f, File const & dst_f, std::string const & p = std::string())
    {
        assert(src_f.is_open());
        assert(dst_f.is_open());
        assert(dst_f.is_rw());
        rec_copy_attributes(src_f, dst_f, p);
    }

    /**
     * Pack raw samples
     */
    static Raw_Samples_Pack pack_rw(std::vector< Raw_Samples_Int_Entry > const & rsi)
    {
        Raw_Samples_Pack rsp;
        std::tie(rsp.signal, rsp.signal_param) = rw_coder().encode(rsi, true);
        return rsp;
    }
    /**
     * Unpack raw samples
     */
    static std::vector< Raw_Samples_Int_Entry > unpack_rw(Raw_Samples_Pack const & rsp)
    {
        return rw_coder().decode< Raw_Samples_Int_Entry >(rsp.signal, rsp.signal_param);
    }
    /**
     * Pack eventdetection events
     */
    static EventDetection_Events_Pack pack_ed(
        std::vector< EventDetection_Event_Entry > const & ed,
        EventDetection_Event_Parameters const & ed_param)
    {
        EventDetection_Events_Pack ed_pack;
        std::vector< long long > skip;
        std::vector< long long > len;
        long long last_end = ed_param.start_time;
        for (unsigned i = 0; i < ed.size(); ++i)
        {
            skip.push_back(ed[i].start - last_end);
            len.push_back(ed[i].length);
            last_end = ed[i].start + ed[i].length;
        }
        std::tie(ed_pack.skip, ed_pack.skip_param) = ed_skip_coder().encode(skip, false);
        std::tie(ed_pack.len, ed_pack.len_param) = ed_len_coder().encode(len, false);
        return ed_pack;
    }
    /**
     * Unpack eventdetection events
     */
    static std::vector< EventDetection_Event_Entry > unpack_ed(
        EventDetection_Events_Pack const & ed_pack,
        EventDetection_Event_Parameters const & ed_param,
        std::vector< Raw_Samples_Entry > const & rs,
        Raw_Samples_Parameters const & rs_param)
    {
        auto skip = ed_skip_coder().decode< long long >(ed_pack.skip, ed_pack.skip_param);
        auto len = ed_len_coder().decode< long long >(ed_pack.len, ed_pack.len_param);
        if (skip.size() != len.size())
        {
            throw std::invalid_argument("unpack_ed failure: skip and length of different size");
        }
        std::vector< EventDetection_Event_Entry > ed(skip.size());
        long long last_end = ed_param.start_time;
        long long off_by_one = ed_param.start_time == rs_param.start_time; // hack
        for (unsigned i = 0; i < skip.size(); ++i)
        {
            ed[i].start = last_end + skip[i];
            ed[i].length = len[i];
            last_end = ed[i].start + ed[i].length;
            // use rs to reconstruct mean and stdv
            long long rs_start_idx = ed[i].start - rs_param.start_time + off_by_one;
            if (rs_start_idx < 0 or rs_start_idx + ed[i].length > (long long)rs.size() + off_by_one)
            {
                throw std::invalid_argument("unpack_ed failure: bad rs_start_idx");
            }
            double s = 0.0;
            double s2 = 0.0;
            unsigned n = ed[i].length;
            if (off_by_one and rs_start_idx + n == (long long)rs.size()) --n;
            for (unsigned j = 0; j < n; ++j)
            {
                auto x = rs[rs_start_idx + j];
                s += x;
                s2 += x * x;
            }
            ed[i].mean = s / n;
            ed[i].stdv = n > 1
                ? std::sqrt((s2 - s*s/n)/(n))
                : 0;
        }
        return ed;
    }

    /**
     * Pack basecall fastq
     */
    static Basecall_Fastq_Pack pack_fq(std::string const & fq, unsigned qv_bits = 5)
    {
        static unsigned const max_qv_bits = 5;
        static std::uint8_t const max_qv = ((std::uint8_t)1 << max_qv_bits) - 1;
        Basecall_Fastq_Pack fq_pack;
        auto fqa = split_fq(fq);
        fq_pack.read_name = fqa[0];
        std::vector< std::int8_t > bp(fqa[1].begin(), fqa[1].end());
        qv_bits = std::min(qv_bits, max_qv_bits);
        auto qv_mask = max_qv & (max_qv << (max_qv_bits - qv_bits));
        fq_pack.qv_bits = qv_bits;
        std::vector< std::uint8_t > qv;
        for (auto c : fqa[3])
        {
            std::uint8_t val = (std::uint8_t)(c - 33);
            val = std::min(val, max_qv);
            val &= qv_mask;
            qv.push_back(val);
        }
        std::tie(fq_pack.bp, fq_pack.bp_param) = fq_bp_coder().encode(bp, false);
        std::tie(fq_pack.qv, fq_pack.qv_param) = fq_qv_coder().encode(qv, false);
        return fq_pack;
    }
    static std::string unpack_fq(Basecall_Fastq_Pack const & fq_pack)
    {
        std::string res;
        res += "@";
        res += fq_pack.read_name;
        res += "\n";
        auto bp = fq_bp_coder().decode< std::int8_t >(fq_pack.bp, fq_pack.bp_param);
        for (auto c : bp) res += c;
        res += "\n+\n";
        auto qv = fq_qv_coder().decode< std::uint8_t >(fq_pack.qv, fq_pack.qv_param);
        for (auto c : qv) res += (char)33 + c;
        res += "\n";
        return res;
    }

    /**
     * Pack basecall events
     */
    static Basecall_Events_Pack pack_ev(
        std::vector< Event_Entry > const & ev,
        std::string const & fq,
        Basecall_Event_Parameters const & ev_param,
        std::vector< EventDetection_Event_Entry > const & ed,
        std::string const & ed_gr,
        double sampling_rate)
    {
        Basecall_Events_Pack ev_pack;
        ev_pack.ed_gr = ed_gr;
        ev_pack.start_time = time_to_int(ev_param.start_time, sampling_rate);
        ev_pack.duration = time_to_int(ev_param.duration, sampling_rate);
        ev_pack.state_size = ev[0].get_model_state().size();

        auto fqa = split_fq(fq);
        std::vector< long long > skip;
        std::vector< std::uint8_t > mv;
        std::string bases;
        long long j = 0;
        for (unsigned i = 0; i < ev.size(); ++i)
        {
            // skip
            auto ti = time_to_int(ev[i].start, sampling_rate);
            auto last_j = j;
            while (j < (long long)ed.size() and ed[j].start < ti) ++j;
            if (j == (long long)ed.size())
            {
                throw std::invalid_argument("pack_ev failed: no matching ed event");
            }
            skip.push_back(j - last_j - 1);
            // move
            if (ev[i].move < 0 or ev[i].move > std::numeric_limits< uint8_t >::max())
            {
                throw std::invalid_argument("pack_ev failed: invalid move");
            }
            mv.push_back(ev[i].move);
            // state
            auto s = ev[i].get_model_state();
            if (s.size() != ev_pack.state_size)
            {
                throw std::invalid_argument("pack_ev failed: different state sizes");
            }
            int bases_to_skip = 0;
            if (i > 0 and ev[i].move < (int)ev_pack.state_size)
            {
                bases_to_skip = (int)ev_pack.state_size - ev[i].move;
            }
            for (auto c : s.substr(bases_to_skip))
            {
                bases.push_back(c);
            }
        }
        if (bases != fqa[1])
        {
            throw std::invalid_argument("pack_ev failed: sequence of states does not match fastq seq");
        }

        std::tie(ev_pack.skip, ev_pack.skip_param) = ev_skip_coder().encode(skip, false);
        std::tie(ev_pack.move, ev_pack.move_param) = ev_move_coder().encode(mv, false);
        return ev_pack;
    }

    static std::vector< Event_Entry >
    unpack_ev(Basecall_Events_Pack const & ev_pack,
              std::string const & fq,
              std::vector< EventDetection_Event_Entry > const & ed,
              double sampling_rate)
    {
        std::vector< Event_Entry > res;
        auto skip = ev_skip_coder().decode< long long >(ev_pack.skip, ev_pack.skip_param);
        auto mv = ev_move_coder().decode< std::uint8_t >(ev_pack.move, ev_pack.move_param);
        //auto bases = ev_bases_coder().decode< std::int8_t >(ev_pack.bases, ev_pack.bases_param);
        auto fqa = split_fq(fq);
        auto const & bases = fqa[1];
        if (skip.size() != mv.size())
        {
            throw std::invalid_argument("unpack_ev failed: skip and move have different sizes");
        }
        res.resize(skip.size());
        long long j = 0;
        std::string s;
        unsigned bases_pos = 0;
        for (unsigned i = 0; i < res.size(); ++i)
        {
            j += skip[i] + 1;
            res[i].start = time_to_float(ed[j].start, sampling_rate);
            res[i].length = time_to_float(ed[j].length, sampling_rate);
            res[i].mean = ed[j].mean;
            res[i].stdv = ed[j].stdv;
            res[i].move = mv[i];
            if (i > 0) s = s.substr(mv[i]); // apply move
            while (s.size() < ev_pack.state_size) s += bases[bases_pos++];
            std::copy(s.begin(), s.end(), res[i].model_state.begin());
            if (ev_pack.state_size < MAX_K_LEN) res[i].model_state[ev_pack.state_size] = 0;
        }
        return res;
    }

    static long long time_to_int(double tf, double sampling_rate)
    {
        return tf * sampling_rate + .5;
    }

    static double time_to_float(long long ti, double sampling_rate)
    {
        return (long double)ti / sampling_rate;
    }

private:
    // channel id params, including sampling rate
    Channel_Id_Parameters _channel_id_params;

    // list of read names for which we have raw samples
    std::vector< std::string > _raw_samples_read_name_list;

    // list of EventDetection groups
    std::vector< std::string > _eventdetection_group_list;

    // list of Basecall groups
    std::vector< std::string > _basecall_group_list;

    // list of per-strand Basecall groups; 0/1/2 = template/complement/2d
    std::array< std::vector< std::string >, 3 > _basecall_strand_group_list;

    void load_channel_id_params()
    {
        if (not Base::group_exists(channel_id_path())) return;
        Base::read(channel_id_path() + "/channel_number", _channel_id_params.channel_number);
        Base::read(channel_id_path() + "/digitisation", _channel_id_params.digitisation);
        Base::read(channel_id_path() + "/offset", _channel_id_params.offset);
        Base::read(channel_id_path() + "/range", _channel_id_params.range);
        Base::read(channel_id_path() + "/sampling_rate", _channel_id_params.sampling_rate);
    }

    void detect_raw_samples_read_name_list()
    {
        if (not Base::group_exists(raw_samples_root_path())) return;
        auto rn_list = Base::list_group(raw_samples_root_path());
        for (auto const & rn : rn_list)
        {
            if (Base::dataset_exists(raw_samples_path(rn))
                or Base::group_exists(raw_samples_pack_path(rn)))
            {
                _raw_samples_read_name_list.push_back(rn);
            }
        }
    }

    void detect_eventdetection_group_list()
    {
        if (not Base::group_exists(eventdetection_root_path())) return;
        auto g_list = Base::list_group(eventdetection_root_path());
        for (auto const & g : g_list)
        {
            if (g.size() <= eventdetection_group_prefix().size()) continue;
            auto p = std::mismatch(eventdetection_group_prefix().begin(),
                                   eventdetection_group_prefix().end(),
                                   g.begin());
            if (p.first != eventdetection_group_prefix().end()) continue;
            _eventdetection_group_list.emplace_back(p.second, g.end());
        }
    }

    std::vector< std::string > detect_eventdetection_read_name_list(std::string const & gr) const
    {
        std::vector< std::string > res;
        std::string p = eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr + "/Reads";
        if (not Base::group_exists(p)) return res;
        auto rn_list = Base::list_group(p);
        for (auto const & rn : rn_list)
        {
            if (Base::dataset_exists(eventdetection_events_path(gr, rn))
                or Base::group_exists(eventdetection_events_pack_path(gr, rn)))
            {
                res.push_back(rn);
            }
        }
        return res;
    }

    void detect_basecall_group_list()
    {
        if (not Base::group_exists(basecall_root_path())) return;
        auto g_list = Base::list_group(basecall_root_path());
        for (auto const & g : g_list)
        {
            if (g.size() <= basecall_group_prefix().size()) continue;
            auto p = std::mismatch(basecall_group_prefix().begin(),
                                   basecall_group_prefix().end(),
                                   g.begin());
            if (p.first != basecall_group_prefix().end()) continue;
            _basecall_group_list.emplace_back(p.second, g.end());
            for (unsigned st = 0; st < 3; ++st)
            {
                if (Base::group_exists(basecall_root_path() + "/" + g + "/" + basecall_strand_subgroup(st)))
                {
                    _basecall_strand_group_list[st].emplace_back(p.second, g.end());
                }
            }
        }
    }

    Attr_Map get_attr_map(std::string const & path) const
    {
        Attr_Map res;
        auto a_list = Base::get_attr_list(path);
        for (auto const & a : a_list)
        {
            std::string tmp;
            Base::read(path + "/" + a, tmp);
            res[a] = tmp;
        }
        return res;
    }

    void add_attr_map(std::string const & path, Attr_Map const & attr_m) const
    {
        for (auto const & p : attr_m)
        {
            Base::write_attribute(path + "/" + p.first, p.second);
        }
    }

    float raw_sample_to_float(int si) const
    {
        assert(have_channel_id_params());
        return ((float)si + _channel_id_params.offset)
            * _channel_id_params.range / _channel_id_params.digitisation;
    }

    // static paths
    static std::string const & file_version_path()
    {
        static const std::string _file_version_path = "/file_version";
        return _file_version_path;
    }

    static std::string const & channel_id_path()
    {
        static const std::string _channel_id_path = "/UniqueGlobalKey/channel_id";
        return _channel_id_path;
    }
    static std::string const & tracking_id_path()
    {
        static const std::string _tracking_id_path = "/UniqueGlobalKey/tracking_id";
        return _tracking_id_path;
    }
    static std::string const & raw_samples_root_path()
    {
        static const std::string _raw_samples_root_path = "/Raw/Reads";
        return _raw_samples_root_path;
    }
    static std::string raw_samples_params_path(std::string const & rn)
    {
        return raw_samples_root_path() + "/" + rn;
    }
    static std::string raw_samples_path(std::string const & rn)
    {
        return raw_samples_root_path() + "/" + rn + "/Signal";
    }
    static std::string raw_samples_pack_path(std::string const & rn)
    {
        return raw_samples_root_path() + "/" + rn + "/Signal_Pack";
    }
    static std::string const & sequences_path()
    {
        static const std::string _sequences_path = "/Sequences/Meta";
        return _sequences_path;
    }
    static std::string const & eventdetection_root_path()
    {
        static const std::string _eventdetection_root_path = "/Analyses";
        return _eventdetection_root_path;
    }
    static std::string const & eventdetection_group_prefix()
    {
        static const std::string _eventdetection_group_prefix = "EventDetection_";
        return _eventdetection_group_prefix;
    }
    static std::string eventdetection_params_path(std::string const & gr)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr;
    }
    static std::string eventdetection_event_params_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr + "/Reads/" + rn;
    }
    static std::string eventdetection_events_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr + "/Reads/" + rn + "/Events";
    }
    static std::string eventdetection_events_pack_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr + "/Reads/" + rn + "/Events_Pack";
    }

    static std::string const & basecall_root_path()
    {
        static const std::string _basecall_root_path = "/Analyses";
        return _basecall_root_path;
    }
    static std::string const & basecall_group_prefix()
    {
        static const std::string _basecall_group_prefix = "Basecall_";
        return _basecall_group_prefix;
    }
    static std::string const & basecall_strand_subgroup(unsigned st)
    {
        static const std::array< std::string, 3 > _basecall_strand_subgroup =
            {{ "BaseCalled_template", "BaseCalled_complement", "BaseCalled_2D" }};
        return _basecall_strand_subgroup[st];
    }
    static std::string basecall_fastq_path(std::string const & gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr + "/"
            + basecall_strand_subgroup(st) + "/Fastq";
    }
    static std::string basecall_fastq_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr + "/"
            + basecall_strand_subgroup(st) + "/Fastq_Pack";
    }
    static std::string basecall_model_path(std::string const & gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr + "/"
            + basecall_strand_subgroup(st) + "/Model";
    }
    static std::string basecall_model_file_path(std::string const & gr, unsigned st)
    {
        assert(st < 2);
        return basecall_root_path() + "/" + basecall_group_prefix() + gr
            + "/Summary/basecall_1d_" + (st == 0? "template" : "complement") + "/model_file";
    }
    static std::string basecall_events_path(std::string const & gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr + "/"
            + basecall_strand_subgroup(st) + "/Events";
    }
    static std::string basecall_events_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr + "/"
            + basecall_strand_subgroup(st) + "/Events_Pack";
    }
    static std::string basecall_event_alignment_path(std::string const & gr)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr + "/"
            + basecall_strand_subgroup(2) + "/Alignment";
    }

    static void rec_copy_attributes(File const & src_f, File const & dst_f, std::string const & path)
    {
        Base const & src_fb = src_f;
        Base const & dst_fb = dst_f;
        auto a_l = src_fb.get_attr_list(not path.empty()? path : std::string("/"));
        for (auto const & a : a_l)
        {
            Base::copy_attribute(src_fb, dst_fb, path + "/" + a);
        }
        auto sg_l = src_fb.list_group(not path.empty()? path : std::string("/"));
        for (auto const & sg : sg_l)
        {
            if (src_fb.group_exists(path + "/" + sg))
            {
                rec_copy_attributes(src_f, dst_f, path + "/" + sg);
            }
            /*
            else
            {
                std::clog << "skipping dataset: " << path + "/" + sg << std::endl;
            }
            */
        }
    } // rec_copy_attributes()

    static fast5_pack::Huffman_Coder & rw_coder()
    {
        static fast5_pack::Huffman_Coder _rw_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_rw_1.inl"
                ;
            _rw_coder.load_codeword_map(tmp, "fast5_cwmap_rw_1.inl");
            initialized = true;
        }
        return _rw_coder;
    } // rw_coder()
    
    static fast5_pack::Huffman_Coder & ed_skip_coder()
    {
        static fast5_pack::Huffman_Coder _ed_skip_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_ed_skip_1.inl"
                ;
            _ed_skip_coder.load_codeword_map(tmp, "fast5_cwmap_ed_skip_1.inl");
            initialized = true;
        }
        return _ed_skip_coder;
    } // ed_skip_coder()

    static fast5_pack::Huffman_Coder & ed_len_coder()
    {
        static fast5_pack::Huffman_Coder _ed_len_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_ed_len_1.inl"
                ;
            _ed_len_coder.load_codeword_map(tmp, "fast5_cwmap_ed_len_1.inl");
            initialized = true;
        }
        return _ed_len_coder;
    } // ed_len_coder()

    static fast5_pack::Huffman_Coder & ev_skip_coder()
    {
        static fast5_pack::Huffman_Coder _ev_skip_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_ev_skip_1.inl"
                ;
            _ev_skip_coder.load_codeword_map(tmp, "fast5_cwmap_ev_skip_1.inl");
            initialized = true;
        }
        return _ev_skip_coder;
    } // ev_skip_coder()

    static fast5_pack::Huffman_Coder & ev_move_coder()
    {
        static fast5_pack::Huffman_Coder _ev_move_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_ev_move_1.inl"
                ;
            _ev_move_coder.load_codeword_map(tmp, "fast5_cwmap_ev_move_1.inl");
            initialized = true;
        }
        return _ev_move_coder;
    } // ev_move_coder()

    static fast5_pack::Huffman_Coder & fq_bp_coder()
    {
        static fast5_pack::Huffman_Coder _fq_bp_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_fq_bp_1.inl"
                ;
            _fq_bp_coder.load_codeword_map(tmp, "fast5_cwmap_fq_bp_1.inl");
            initialized = true;
        }
        return _fq_bp_coder;
    } // fq_bp_coder()

    static fast5_pack::Huffman_Coder & fq_qv_coder()
    {
        static fast5_pack::Huffman_Coder _fq_qv_coder;
        static bool initialized = false;
        if (not initialized)
        {
            std::vector< std::string > tmp =
#include "fast5_cwmap_fq_qv_1.inl"
                ;
            _fq_qv_coder.load_codeword_map(tmp, "fast5_cwmap_fq_qv_1.inl");
            initialized = true;
        }
        return _fq_qv_coder;
    } // fq_qv_coder()

}; // class File

} // namespace fast5

#endif
