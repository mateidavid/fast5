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
#include "Huffman_Packer.hpp"
#include "Bit_Packer.hpp"

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

typedef hdf5_tools::File::Attr_Map Attr_Map;

struct Channel_Id_Params
{
    std::string channel_number;
    double digitisation;
    double offset;
    double range;
    double sampling_rate;
    Channel_Id_Params()
        : channel_number(""),
          digitisation(0.0),
          offset(0.0),
          range(0.0),
          sampling_rate(0.0) {}
}; // struct Channel_Id_Params

typedef Attr_Map Tracking_Id_Params;

typedef Attr_Map Sequences_Params;

typedef float Raw_Sample;
typedef int16_t Raw_Int_Sample;

struct Raw_Samples_Params
{
    std::string read_id;
    long long read_number;
    long long start_mux;
    long long start_time;
    long long duration;
}; // struct Raw_Samples_Params

struct Raw_Samples_Pack
{
    Huffman_Packer::Code_Type signal;
    Attr_Map signal_params;
}; // struct Raw_Samples_Pack

struct EventDetection_Event
{
    double mean;
    double stdv;
    long long start;
    long long length;
    friend bool operator == (EventDetection_Event const & lhs, EventDetection_Event const & rhs)
    {
        return lhs.mean == rhs.mean
            and lhs.stdv == rhs.stdv
            and lhs.start == rhs.start
            and lhs.length == rhs.length;
    }
}; // struct EventDetection_Event

struct EventDetection_Events_Params
{
    std::string read_id;
    long long read_number;
    long long scaling_used;
    long long start_mux;
    long long start_time;
    long long duration;
    double median_before;
    unsigned abasic_found;
}; // struct EventDetection_Events_Params

struct EventDetection_Events_Pack
{
    Huffman_Packer::Code_Type skip;
    Attr_Map skip_params;
    Huffman_Packer::Code_Type len;
    Attr_Map len_params;
}; // struct EventDetection_Events_Pack

//
// This struct represents the expected signal measured
// given the kmer sequence that is in the pore when the
// the observations are made. A pore model consists
// of 1024 of these entries (one per 5-mer) and global
// shift/scaling params.
//
struct Basecall_Model_State
{
    long long variant;
    double level_mean;
    double level_stdv;
    double sd_mean;
    double sd_stdv;
    double weight;
    std::array< char, MAX_K_LEN > kmer;
    std::string get_kmer() const { return array_to_string(kmer); }
    friend bool operator == (Basecall_Model_State const & lhs, Basecall_Model_State const & rhs)
    {
        return lhs.variant == rhs.variant
            and lhs.level_mean == rhs.level_mean
            and lhs.level_stdv == rhs.level_stdv
            and lhs.sd_mean == rhs.sd_mean
            and lhs.sd_stdv == rhs.sd_stdv
            and lhs.weight == rhs.weight
            and lhs.kmer == rhs.kmer;
    }
}; // struct Basecall_Model_State

//
// This struct represents the global transformations
// that must be applied to each Basecall_Model_State
//
struct Basecall_Model_Params
{
    double scale;
    double shift;
    double drift;
    double var;
    double scale_sd;
    double var_sd;
}; // struct Basecall_Model_Params

struct Basecall_Fastq_Pack
{
    Huffman_Packer::Code_Type bp;
    Attr_Map bp_params;
    Huffman_Packer::Code_Type qv;
    Attr_Map qv_params;
    std::string read_name;
    std::uint8_t qv_bits;
}; // struct Basecall_Fastq_Pack

//
// This struct represents an observed event.
// The members of the struct are the same as
// the fields encoded in the FAST5 file.
//
struct Basecall_Event
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
    friend bool operator == (Basecall_Event const & lhs, Basecall_Event const & rhs)
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
}; // struct Basecall_Event

struct Basecall_Events_Params
{
    double start_time;
    double duration;
};

struct Basecall_Events_Pack
{
    Huffman_Packer::Code_Type skip;
    Attr_Map skip_params;
    Huffman_Packer::Code_Type move;
    Attr_Map move_params;
    Bit_Packer::Code_Type p_model_state;
    Attr_Map p_model_state_params;
    //
    std::string ed_gr;
    long long start_time;
    long long duration;
    unsigned state_size;
    unsigned p_model_state_bits;
}; // struct Basecall_Events_Pack

//
// This struct represents a template-to-complement
// match that is emitted by ONT's 2D basecaller
//
struct Basecall_Alignment_Entry
{
    long long template_index;
    long long complement_index;
    std::array< char, MAX_K_LEN > kmer;
    std::string get_kmer() const { return array_to_string(kmer); }
    friend bool operator == (Basecall_Alignment_Entry const & lhs, Basecall_Alignment_Entry const & rhs)
    {
        return lhs.template_index == rhs.template_index
            and lhs.complement_index == rhs.complement_index
            and lhs.kmer == rhs.kmer;
    }
}; // struct Basecall_Alignment_Entry

struct Basecall_Alignment_Pack
{
    Bit_Packer::Code_Type template_step;
    Bit_Packer::Code_Params_Type template_step_params;
    Bit_Packer::Code_Type complement_step;
    Bit_Packer::Code_Params_Type complement_step_params;
    Huffman_Packer::Code_Type move;
    Huffman_Packer::Code_Params_Type move_params;
    unsigned template_index_start;
    unsigned complement_index_start;
    unsigned kmer_size;
};

class File
    : private hdf5_tools::File
{
private:
    typedef hdf5_tools::File Base;
public:
    //
    // Constructors
    //
    File() = default;
    File(std::string const & file_name, bool rw = false) { open(file_name, rw); }

    //
    // Base methods
    //
    using Base::is_open;
    using Base::is_rw;
    using Base::file_name;
    using Base::create;
    using Base::close;
    using Base::get_object_count;
    using Base::is_valid_file;

    //
    // Base method wrappers
    //
    void
    open(std::string const & file_name, bool rw = false)
    {
        Base::open(file_name, rw);
        reload();
    }
    static void
    copy_attributes(File const & src_f, File const & dst_f, std::string const & p, bool recurse = false)
    {
        Base::copy_attributes(src_f, dst_f, p, recurse);
    }

    //
    // Access /file_version
    //
    std::string
    file_version() const
    {
        std::string res;
        Base::read(file_version_path(), res);
        return res;
    }

    //
    // Access /UniqueGlobalKey/channel_id
    //
    bool
    have_channel_id_params() const
    {
        return _channel_id_params.sampling_rate > 0.0;
    }
    Channel_Id_Params
    get_channel_id_params() const
    {
        return _channel_id_params;
    }
    void
    add_channel_id_params(Channel_Id_Params const & channel_id_params)
    {
        _channel_id_params = channel_id_params;
        Base::write_attribute(channel_id_path() + "/channel_number", _channel_id_params.channel_number);
        Base::write_attribute(channel_id_path() + "/digitisation", _channel_id_params.digitisation);
        Base::write_attribute(channel_id_path() + "/offset", _channel_id_params.offset);
        Base::write_attribute(channel_id_path() + "/range", _channel_id_params.range);
        Base::write_attribute(channel_id_path() + "/sampling_rate", _channel_id_params.sampling_rate);
    }
    bool
    have_sampling_rate() const { return have_channel_id_params(); }
    double
    get_sampling_rate() const { return _channel_id_params.sampling_rate; }

    //
    // Access /UniqueGlobalKey/tracking_id
    //
    bool
    have_tracking_id_params() const
    {
        return Base::group_exists(tracking_id_path());
    }
    Tracking_Id_Params
    get_tracking_id_params() const
    {
        return get_attr_map(tracking_id_path());
    }
    void
    add_tracking_id_params(Tracking_Id_Params const & tracking_id_params) const
    {
        add_attr_map(tracking_id_path(), tracking_id_params);
    }

    //
    // Access /Sequences
    //
    bool
    have_sequences_params() const
    {
        return Base::group_exists(sequences_path());
    }
    Sequences_Params
    get_sequences_params() const
    {
        return get_attr_map(sequences_path());
    }
    void
    add_sequences_params(Sequences_Params const & sequences_params) const
    {
        add_attr_map(sequences_path(), sequences_params);
    }

    //
    // Access Raw Samples
    //
    std::vector< std::string > const &
    get_raw_samples_read_name_list() const
    {
        return _raw_samples_read_names;
    }
    bool
    have_raw_samples(std::string const & rn = std::string()) const
    {
        auto && rn_l = get_raw_samples_read_name_list();
        return (rn.empty()
                ? not rn_l.empty()
                : std::find(rn_l.begin(), rn_l.end(), rn) != rn_l.end());
    }
    bool
    have_raw_samples_unpack(std::string const & rn) const
    {
        return Base::dataset_exists(raw_samples_path(rn));
    }
    bool
    have_raw_samples_pack(std::string const & rn) const
    {
        return Base::group_exists(raw_samples_pack_path(rn));
    }
    Raw_Samples_Params
    get_raw_samples_params(std::string const & rn = std::string()) const
    {
        Raw_Samples_Params res;
        auto && _rn = fill_raw_samples_read_name(rn);
        std::string p = raw_samples_params_path(_rn);
        Base::read(p + "/read_id", res.read_id);
        Base::read(p + "/read_number", res.read_number);
        Base::read(p + "/start_mux", res.start_mux);
        Base::read(p + "/start_time", res.start_time);
        Base::read(p + "/duration", res.duration);
        return res;
    }
    void
    add_raw_samples_params(std::string const & rn, Raw_Samples_Params const & params) const
    {
        std::string p = raw_samples_params_path(rn);
        Base::write_attribute(p + "/read_id", params.read_id);
        Base::write_attribute(p + "/read_number", params.read_number);
        Base::write_attribute(p + "/start_mux", params.start_mux);
        Base::write_attribute(p + "/start_time", params.start_time);
        Base::write_attribute(p + "/duration", params.duration);
    }
    std::vector< Raw_Int_Sample >
    get_raw_int_samples(std::string const & rn = std::string()) const
    {
        std::vector< Raw_Int_Sample > res;
        auto && _rn = fill_raw_samples_read_name(rn);
        if (have_raw_samples_unpack(_rn))
        {
            Base::read(raw_samples_path(_rn), res);
        }
        else if (have_raw_samples_pack(_rn))
        {
            auto rs_pack = get_raw_samples_pack(_rn);
            res = unpack_rw(rs_pack);
        }
        return res;
    }
    void
    add_raw_samples(std::string const & rn, std::vector< Raw_Int_Sample > const & rsi)
    {
        Base::write_dataset(raw_samples_path(rn), rsi);
        reload();
    }
    std::vector< Raw_Sample >
    get_raw_samples(std::string const & rn = std::string()) const
    {
        // get raw samples
        auto rsi = get_raw_int_samples(rn);
        // decode levels
        std::vector< Raw_Sample > res;
        res.reserve(rsi.size());
        for (auto int_level : rsi)
        {
            res.push_back(raw_int_sample_to_float(int_level, _channel_id_params));
        }
        return res;
    }
    Raw_Samples_Pack
    get_raw_samples_pack(std::string const & rn) const
    {
        Raw_Samples_Pack rs_pack;
        Base::read(raw_samples_pack_path(rn) + "/Signal", rs_pack.signal);
        rs_pack.signal_params = get_attr_map(raw_samples_pack_path(rn) + "/Signal");
        return rs_pack;
    }
    void
    add_raw_samples(std::string const & rn, Raw_Samples_Pack const & rs_pack)
    {
        Base::write_dataset(raw_samples_pack_path(rn) + "/Signal", rs_pack.signal);
        add_attr_map(raw_samples_pack_path(rn) + "/Signal", rs_pack.signal_params);
        reload();
    }

    //
    // Access EventDetection groups
    //
    std::vector< std::string > const &
    get_eventdetection_group_list() const
    {
        return _eventdetection_groups;
    }
    bool
    have_eventdetection_group(std::string const & gr = std::string()) const
    {
        return (gr.empty()
                ? not _eventdetection_groups.empty()
                : _eventdetection_read_names.count(gr));
    }
    std::vector< std::string > const &
    get_eventdetection_read_name_list(std::string const & gr = std::string()) const
    {
        static const std::vector< std::string > _empty;
        auto && _gr = fill_eventdetection_group(gr);
        return (_eventdetection_read_names.count(_gr)
                ? _eventdetection_read_names.at(_gr)
                : _empty);
    }
    bool
    have_eventdetection_events(
        std::string const & gr = std::string(), std::string const & rn = std::string()) const
    {
        auto && _gr = fill_eventdetection_group(gr);
        auto && _rn = fill_eventdetection_read_name(_gr, rn);
        return (_eventdetection_read_names.count(_gr)
                and std::find(
                    _eventdetection_read_names.at(_gr).begin(),
                    _eventdetection_read_names.at(_gr).end(),
                    _rn)
                != _eventdetection_read_names.at(_gr).end());
    }

    //
    // Access EventDetection group params
    //
    Attr_Map
    get_eventdetection_params(std::string const & gr = std::string()) const
    {
        auto && _gr = fill_eventdetection_group(gr);
        return get_attr_map(eventdetection_group_path(_gr));
    }
    void
    add_eventdetection_params(std::string const & gr, Attr_Map const & am) const
    {
        add_attr_map(eventdetection_group_path(gr), am);
    }

    //
    // Access EventDetection events params
    //
    EventDetection_Events_Params
    get_eventdetection_events_params(
        std::string const & gr = std::string(), std::string const & rn = std::string()) const
    {
        EventDetection_Events_Params res;
        auto && _gr = fill_eventdetection_group(gr);
        auto && _rn = fill_eventdetection_read_name(_gr, rn);
        auto p = eventdetection_events_params_path(_gr, _rn);
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
    void
    add_eventdetection_events_params(
        std::string const & gr, std::string const & rn,
        EventDetection_Events_Params const & ede_params) const
    {
        auto p = eventdetection_events_params_path(gr, rn);
        if (not ede_params.read_id.empty()) Base::write_attribute(p + "/read_id", ede_params.read_id);
        Base::write_attribute(p + "/read_number", ede_params.read_number);
        Base::write_attribute(p + "/scaling_used", ede_params.scaling_used);
        Base::write_attribute(p + "/start_mux", ede_params.start_mux);
        Base::write_attribute(p + "/start_time", ede_params.start_time);
        Base::write_attribute(p + "/duration", ede_params.duration);
        if (ede_params.median_before > 0.0) Base::write_attribute(p + "/median_before", ede_params.median_before);
        if (ede_params.abasic_found < 2) Base::write_attribute(p + "/abasic_found", ede_params.abasic_found);
    }

    //
    // Access EventDetection events
    //
    bool
    have_eventdetection_events_unpack(std::string const & gr, std::string const & rn) const
    {
        return Base::dataset_exists(eventdetection_events_path(gr, rn));
    }
    bool
    have_eventdetection_events_pack(std::string const & gr, std::string const & rn) const
    {
        return Base::group_exists(eventdetection_events_pack_path(gr, rn));
    }
    std::vector< EventDetection_Event >
    get_eventdetection_events(
        std::string const & gr = std::string(), std::string const & rn = std::string()) const
    {
        std::vector< EventDetection_Event > res;
        auto && _gr = fill_eventdetection_group(gr);
        auto && _rn = fill_eventdetection_read_name(_gr, rn);
        if (have_eventdetection_events_unpack(_gr, _rn))
        {
            auto p = eventdetection_events_path(_gr, _rn);
            auto struct_member_names = Base::get_struct_members(p);
            bool have_stdv = false;
            bool have_variance = false;
            for (auto const & s : struct_member_names)
            {
                if (s == "stdv") have_stdv = true;
                else if (s == "variance") have_variance = true;
            }
            hdf5_tools::Compound_Map m;
            m.add_member("mean", &EventDetection_Event::mean);
            m.add_member("start", &EventDetection_Event::start);
            m.add_member("length", &EventDetection_Event::length);
            if (have_stdv)
            {
                m.add_member("stdv", &EventDetection_Event::stdv);
            }
            else if (have_variance)
            {
                m.add_member("variance", &EventDetection_Event::stdv);
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
            auto ed_pack = get_eventdetection_events_pack(_gr, _rn);
            auto ed_params = get_eventdetection_events_params(_gr, _rn);
            auto rs = get_raw_samples(_rn);
            auto rs_params = get_raw_samples_params(_rn);
            res = unpack_ed(ed_pack, ed_params, rs, rs_params);
        }
        return res;
    } // get_eventdetection_events()
    void
    add_eventdetection_events(
        std::string const & gr, std::string const & rn,
        std::vector< EventDetection_Event > const & ed)
    {
        hdf5_tools::Compound_Map m;
        m.add_member("mean", &EventDetection_Event::mean);
        m.add_member("start", &EventDetection_Event::start);
        m.add_member("length", &EventDetection_Event::length);
        m.add_member("stdv", &EventDetection_Event::stdv);
        Base::write_dataset(eventdetection_events_path(gr, rn), ed, m);
        reload();
    }
    EventDetection_Events_Pack
    get_eventdetection_events_pack(
        std::string const & gr, std::string const & rn) const
    {
        EventDetection_Events_Pack ed_pack;
        Base::read(eventdetection_events_pack_path(gr, rn) + "/Skip", ed_pack.skip);
        ed_pack.skip_params = get_attr_map(eventdetection_events_pack_path(gr, rn) + "/Skip");
        Base::read(eventdetection_events_pack_path(gr, rn) + "/Len", ed_pack.len);
        ed_pack.len_params = get_attr_map(eventdetection_events_pack_path(gr, rn) + "/Len");
        return ed_pack;
    }
    void
    add_eventdetection_events(
        std::string const & gr, std::string const & rn,
        EventDetection_Events_Pack const & ed_pack)
    {
        Base::write_dataset(eventdetection_events_pack_path(gr, rn) + "/Skip", ed_pack.skip);
        add_attr_map(eventdetection_events_pack_path(gr, rn) + "/Skip", ed_pack.skip_params);
        Base::write_dataset(eventdetection_events_pack_path(gr, rn) + "/Len", ed_pack.len);
        add_attr_map(eventdetection_events_pack_path(gr, rn) + "/Len", ed_pack.len_params);
        reload();
    }

    //
    // Access Basecall groups
    //
    std::vector< std::string > const &
    get_basecall_group_list() const
    {
        return _basecall_groups;
    }
    bool
    have_basecall_group(std::string const & gr = std::string()) const
    {
        auto && gr_l = get_basecall_group_list();
        return (gr.empty()
                ? not gr_l.empty()
                : std::find(gr_l.begin(), gr_l.end(), gr) != gr_l.end());
    }
    std::vector< std::string > const &
    get_basecall_strand_group_list(unsigned st) const
    {
        return _basecall_strand_groups.at(st);
    }
    bool
    have_basecall_strand_group(unsigned st, std::string const & gr = std::string()) const
    {
        auto && gr_l = get_basecall_strand_group_list(st);
        return (gr.empty()
                ? not gr_l.empty()
                : std::find(gr_l.begin(), gr_l.end(), gr) != gr_l.end());
    }
    std::string const &
    get_basecall_1d_group(std::string const & gr) const
    {
        return (_basecall_1d_group.count(gr)
                ? _basecall_1d_group.at(gr)
                : gr);
    }
    std::string const &
    get_basecall_eventdetection_group(std::string const & gr) const
    {
        static std::string const empty;
        return (_basecall_eventdetection_group.count(gr)
                ? _basecall_eventdetection_group.at(gr)
                : empty);
    }

    //
    // Access Basecall group params
    //
    Attr_Map
    get_basecall_params(std::string const & gr) const
    {
        return get_attr_map(basecall_group_path(gr));
    }
    void
    add_basecall_params(std::string const & gr, Attr_Map const & am) const
    {
        add_attr_map(basecall_group_path(gr) + gr, am);
    }
    //
    // Access Basecall group log
    //
    bool
    have_basecall_log(std::string const & gr) const
    {
        return Base::exists(basecall_log_path(gr));
    }
    std::string
    get_basecall_log(std::string const & gr) const
    {
        std::string res;
        Base::read(basecall_log_path(gr), res);
        return res;
    }

    //
    // Access Basecall fastq
    //
    bool
    have_basecall_fastq(unsigned st, std::string const & gr = std::string()) const
    {
        auto && _gr = fill_basecall_group(st, gr);
        return have_basecall_fastq_unpack(st, _gr) or have_basecall_fastq_pack(st, _gr);
    }
    bool
    have_basecall_fastq_unpack(unsigned st, std::string const & gr) const
    {
        return Base::dataset_exists(basecall_fastq_path(gr, st));
    }
    bool
    have_basecall_fastq_pack(unsigned st, std::string const & gr) const
    {
        return Base::group_exists(basecall_fastq_pack_path(gr, st));
    }
    std::string
    get_basecall_fastq(unsigned st, std::string const & gr = std::string()) const
    {
        std::string res;
        auto && _gr = fill_basecall_group(st, gr);
        if (have_basecall_fastq_unpack(st, _gr))
        {
            Base::read(basecall_fastq_path(_gr, st), res);
        }
        else
        {
            auto fq_pack = get_basecall_fastq_pack(st, _gr);
            res = unpack_fq(fq_pack);
        }
        return res;
    }
    void
    add_basecall_fastq(unsigned st, std::string const & gr, std::string const & fq)
    {
        Base::write(basecall_fastq_path(gr, st), true, fq);
        reload();
    }
    Basecall_Fastq_Pack
    get_basecall_fastq_pack(unsigned st, std::string const & gr) const
    {
        Basecall_Fastq_Pack fq_pack;
        auto p = basecall_fastq_pack_path(gr, st);
        Base::read(p + "/BP", fq_pack.bp);
        fq_pack.bp_params = get_attr_map(p + "/BP");
        Base::read(p + "/QV", fq_pack.qv);
        fq_pack.qv_params = get_attr_map(p + "/QV");
        Base::read(p + "/read_name", fq_pack.read_name);
        Base::read(p + "/qv_bits", fq_pack.qv_bits);
        return fq_pack;        
    }
    void
    add_basecall_fastq(unsigned st, std::string const & gr, Basecall_Fastq_Pack const & fq_pack)
    {
        auto p = basecall_fastq_pack_path(gr, st);
        Base::write_dataset(p + "/BP", fq_pack.bp);
        add_attr_map(p + "/BP", fq_pack.bp_params);
        Base::write_dataset(p + "/QV", fq_pack.qv);
        add_attr_map(p + "/QV", fq_pack.qv_params);
        Base::write_attribute(p + "/read_name", fq_pack.read_name);
        Base::write_attribute(p + "/qv_bits", fq_pack.qv_bits);
        reload();
    }
    bool
    have_basecall_seq(unsigned st, std::string const & _gr = std::string()) const
    {
        return have_basecall_fastq(st, _gr);
    }
    std::string
    get_basecall_seq(unsigned st, std::string const & _gr = std::string()) const
    {
        return fq2seq(get_basecall_fastq(st, _gr));
    }
    void
    add_basecall_seq(unsigned st, std::string const & gr,
                     std::string const & name, std::string const & seq, int default_qual = 33)
    {
        std::ostringstream oss;
        oss << "@" << name << "\n"
            << seq << "\n"
            << "+\n"
            << std::string(seq.size(), (char)default_qual);
        add_basecall_fastq(st, gr, oss.str());
        reload();
    }

    //
    // Access Basecall model
    //
    bool
    have_basecall_model(unsigned st, std::string const & gr = std::string()) const
    {
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        return Base::dataset_exists(basecall_model_path(gr_1d, st));
    }
    std::string
    get_basecall_model_file(unsigned st, std::string const & gr = std::string()) const
    {
        std::string res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        Base::read(basecall_model_file_path(gr_1d, st), res);
        return res;
    }
    void
    add_basecall_model_file(unsigned st, std::string const & gr, std::string const & file_name) const
    {
        Base::write_attribute(basecall_model_file_path(gr, st), file_name);
    }
    Basecall_Model_Params
    get_basecall_model_params(unsigned st, std::string const & gr = std::string()) const
    {
        Basecall_Model_Params res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
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
    void
    add_basecall_model_params(unsigned st, std::string const & gr, T const & params) const
    {
        std::string path = basecall_model_path(gr, st);
        Base::write(path + "/scale", false, params.scale);
        Base::write(path + "/shift", false, params.shift);
        Base::write(path + "/drift", false, params.drift);
        Base::write(path + "/var", false, params.var);
        Base::write(path + "/scale_sd", false, params.scale_sd);
        Base::write(path + "/var_sd", false, params.var_sd);
    }
    std::vector< Basecall_Model_State >
    get_basecall_model(unsigned st, std::string const & gr = std::string()) const
    {
        std::vector< Basecall_Model_State > res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        hdf5_tools::Compound_Map m;
        m.add_member("kmer", &Basecall_Model_State::kmer);
        m.add_member("level_mean", &Basecall_Model_State::level_mean);
        m.add_member("level_stdv", &Basecall_Model_State::level_stdv);
        m.add_member("sd_mean", &Basecall_Model_State::sd_mean);
        m.add_member("sd_stdv", &Basecall_Model_State::sd_stdv);
        Base::read(basecall_model_path(gr_1d, st), res, m);
        return res;
    }
    template < typename T >
    void add_basecall_model(unsigned st, std::string const & gr, std::vector< T > const & m)
    {
        hdf5_tools::Compound_Map cm;
        cm.add_member("kmer", &T::kmer);
        cm.add_member("level_mean", &T::level_mean);
        cm.add_member("level_stdv", &T::level_stdv);
        cm.add_member("sd_mean", &T::sd_mean);
        cm.add_member("sd_stdv", &T::sd_stdv);
        auto && gr_1d = get_basecall_1d_group(gr);
        Base::write_dataset(basecall_model_path(gr_1d, st), m, cm);
        reload();
    }

    //
    // Access Basecall events
    //
    bool
    have_basecall_events(unsigned st, std::string const & gr = std::string()) const
    {
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        return (have_basecall_events_unpack(st, gr_1d)
                or have_basecall_events_pack(st, gr_1d));
    }
    bool
    have_basecall_events_unpack(unsigned st, std::string const & gr) const
    {
        return Base::dataset_exists(basecall_events_path(gr, st));
    }
    bool
    have_basecall_events_pack(unsigned st, std::string const & gr) const
    {
        return Base::group_exists(basecall_events_pack_path(gr, st));
    }
    Basecall_Events_Params
    get_basecall_events_params(unsigned st, std::string const & gr = std::string()) const
    {
        Basecall_Events_Params res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        if (have_basecall_events_unpack(st, gr_1d))
        {
            auto path = basecall_events_path(gr_1d, st);
            auto params = get_attr_map(path);
            if (params.count("start_time")) std::istringstream(params["start_time"]) >> res.start_time;
            else res.start_time = 0.0;
            if (params.count("duration")) std::istringstream(params["duration"]) >> res.duration;
            else res.duration = 0.0;
        }
        else if (have_basecall_events_pack(st, gr_1d))
        {
            auto path = basecall_events_pack_path(gr_1d, st);
            Basecall_Events_Pack ev_pack;
            Base::read(path + "/start_time", ev_pack.start_time);
            Base::read(path + "/duration", ev_pack.duration);
            res.start_time = time_to_float(ev_pack.start_time, _channel_id_params.sampling_rate);
            res.duration = time_to_float(ev_pack.duration, _channel_id_params.sampling_rate);
        }
        return res;
    }
    void
    add_basecall_events_params(unsigned st, std::string const & gr,
                               Basecall_Events_Params const & bce_params) const
    {
        auto path = basecall_events_path(gr, st);
        if (not Base::dataset_exists(path))
        {
            throw std::runtime_error("basecall events must be added before their params");
        }
        Base::write_attribute(path + "/start_time", bce_params.start_time);
        Base::write_attribute(path + "/duration", bce_params.duration);
    }
    std::vector< Basecall_Event >
    get_basecall_events(unsigned st, std::string const & gr = std::string()) const
    {
        std::vector< Basecall_Event > res;
        auto && gr_1d = fill_basecall_1d_group(st, gr);
        if (have_basecall_events_unpack(st, gr_1d))
        {
            hdf5_tools::Compound_Map m;
            m.add_member("mean", &Basecall_Event::mean);
            m.add_member("start", &Basecall_Event::start);
            m.add_member("stdv", &Basecall_Event::stdv);
            m.add_member("length", &Basecall_Event::length);
            m.add_member("p_model_state", &Basecall_Event::p_model_state);
            m.add_member("model_state", &Basecall_Event::model_state);
            m.add_member("move", &Basecall_Event::move);
            Base::read(basecall_events_path(gr_1d, st), res, m);
        }
        else if (have_basecall_events_pack(st, gr_1d))
        {
            auto ev_pack = get_basecall_events_pack(st, gr_1d);
            if (not have_basecall_fastq(st, gr_1d))
            {
                throw std::runtime_error("missing fastq for basecall events unpacking");
            }
            auto fq = get_basecall_fastq(st, gr_1d);
            if (not have_eventdetection_events(ev_pack.ed_gr))
            {
                throw std::runtime_error("missing evendetection events for basecall events unpacking");
            }
            auto ed = get_eventdetection_events(ev_pack.ed_gr);
            res = unpack_ev(ev_pack, fq, ed, _channel_id_params.sampling_rate);
        }
        return res;
    }
    template < typename T >
    void
    add_basecall_events(unsigned st, std::string const & gr, std::vector< T > const & ev)
    {
        hdf5_tools::Compound_Map cm;
        cm.add_member("mean", &T::mean);
        cm.add_member("start", &T::start);
        cm.add_member("stdv", &T::stdv);
        cm.add_member("length", &T::length);
        cm.add_member("p_model_state", &T::p_model_state);
        cm.add_member("model_state", &T::model_state);
        cm.add_member("move", &T::move);
        Base::write_dataset(basecall_events_path(gr, st), ev, cm);
        reload();
    }
    Basecall_Events_Pack
    get_basecall_events_pack(unsigned st, std::string const & gr) const
    {
        auto p = basecall_events_pack_path(gr, st);
        Basecall_Events_Pack ev_pack;
        Base::read(p + "/Skip", ev_pack.skip);
        ev_pack.skip_params = get_attr_map(p + "/Skip");
        Base::read(p + "/Move", ev_pack.move);
        ev_pack.move_params = get_attr_map(p + "/Move");
        Base::read(p + "/P_Model_State", ev_pack.p_model_state);
        ev_pack.p_model_state_params = get_attr_map(p + "/P_Model_State");
        Base::read(p + "/ed_gr", ev_pack.ed_gr);
        Base::read(p + "/start_time", ev_pack.start_time);
        Base::read(p + "/duration", ev_pack.duration);
        Base::read(p + "/state_size", ev_pack.state_size);
        Base::read(p + "/p_model_state_bits", ev_pack.p_model_state_bits);
        return ev_pack;
    }
    void
    add_basecall_events(unsigned st, std::string const & gr, Basecall_Events_Pack const & ev_pack)
    {
        auto p = basecall_events_pack_path(gr, st);
        Base::write_dataset(p + "/Skip", ev_pack.skip);
        add_attr_map(p + "/Skip", ev_pack.skip_params);
        Base::write_dataset(p + "/Move", ev_pack.move);
        add_attr_map(p + "/Move", ev_pack.move_params);
        Base::write_dataset(p + "/P_Model_State", ev_pack.p_model_state);
        add_attr_map(p + "/P_Model_State", ev_pack.p_model_state_params);
        Base::write_attribute(p + "/ed_gr", ev_pack.ed_gr);
        Base::write_attribute(p + "/start_time", ev_pack.start_time);
        Base::write_attribute(p + "/duration", ev_pack.duration);
        Base::write_attribute(p + "/state_size", ev_pack.state_size);
        Base::write_attribute(p + "/p_model_state_bits", ev_pack.p_model_state_bits);
        reload();
    }

    //
    // Access Basecall alignment
    //
    bool
    have_basecall_alignment(std::string const & gr = std::string()) const
    {
        auto && _gr = fill_basecall_group(2, gr);
        return have_basecall_alignment_unpack(_gr) or have_basecall_alignment_pack(_gr);
    }
    bool
    have_basecall_alignment_unpack(std::string const & gr) const
    {
        return Base::dataset_exists(basecall_alignment_path(gr));
    }
    bool
    have_basecall_alignment_pack(std::string const & gr) const
    {
        return Base::group_exists(basecall_alignment_pack_path(gr));
    }
    std::vector< Basecall_Alignment_Entry >
    get_basecall_alignment(std::string const & gr = std::string()) const
    {
        std::vector< Basecall_Alignment_Entry > al;
        auto && _gr = fill_basecall_group(2, gr);
        if (have_basecall_alignment_unpack(_gr))
        {
            al = get_basecall_alignment_unpack(_gr);
        }
        else if (have_basecall_alignment_pack(_gr)
                 and have_basecall_seq(2, _gr))
        {
            auto al_pack = get_basecall_alignment_pack(_gr);
            auto seq = get_basecall_seq(2, _gr);
            al = unpack_al(al_pack, seq);
        }
        return al;
    }
    std::vector< Basecall_Alignment_Entry >
    get_basecall_alignment_unpack(std::string const & gr) const
    {
        std::vector< Basecall_Alignment_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("template", &Basecall_Alignment_Entry::template_index);
        m.add_member("complement", &Basecall_Alignment_Entry::complement_index);
        m.add_member("kmer", &Basecall_Alignment_Entry::kmer);
        Base::read(basecall_alignment_path(gr), res, m);
        return res;
    }
    void
    add_basecall_alignment(std::string const & gr, std::vector< Basecall_Alignment_Entry > const & al)
    {
        hdf5_tools::Compound_Map m;
        m.add_member("template", &Basecall_Alignment_Entry::template_index);
        m.add_member("complement", &Basecall_Alignment_Entry::complement_index);
        m.add_member("kmer", &Basecall_Alignment_Entry::kmer);
        Base::write_dataset(basecall_alignment_path(gr), al, m);
        reload();
    }
    Basecall_Alignment_Pack
    get_basecall_alignment_pack(std::string const & gr) const
    {
        Basecall_Alignment_Pack al_pack;
        auto p = basecall_alignment_pack_path(gr);
        Base::read(p + "/Template_Step", al_pack.template_step);
        al_pack.template_step_params = get_attr_map(p + "/Template_Step");
        Base::read(p + "/Complement_Step", al_pack.complement_step);
        al_pack.complement_step_params = get_attr_map(p + "/Complement_Step");
        Base::read(p + "/Move", al_pack.move);
        al_pack.move_params = get_attr_map(p + "/Move");
        Base::read(p + "/template_index_start", al_pack.template_index_start);
        Base::read(p + "/complement_index_start", al_pack.complement_index_start);
        Base::read(p + "/kmer_size", al_pack.kmer_size);
        return al_pack;
    }
    void
    add_basecall_alignment(std::string const & gr, Basecall_Alignment_Pack const & al_pack)
    {
        auto p = basecall_alignment_pack_path(gr);
        Base::write_dataset(p + "/Template_Step", al_pack.template_step);
        add_attr_map(p + "/Template_Step", al_pack.template_step_params);
        Base::write_dataset(p + "/Complement_Step", al_pack.complement_step);
        add_attr_map(p + "/Complement_Step", al_pack.complement_step_params);
        Base::write_dataset(p + "/Move", al_pack.move);
        add_attr_map(p + "/Move", al_pack.move_params);
        Base::write_attribute(p + "/template_index_start", al_pack.template_index_start);
        Base::write_attribute(p + "/complement_index_start", al_pack.complement_index_start);
        Base::write_attribute(p + "/kmer_size", al_pack.kmer_size);
        reload();
    }

    //
    // Packers & Unpackers
    //
    static Raw_Samples_Pack
    pack_rw(std::vector< Raw_Int_Sample > const & rsi)
    {
        Raw_Samples_Pack rsp;
        std::tie(rsp.signal, rsp.signal_params) = rw_coder().encode(rsi, true);
        return rsp;
    }
    static std::vector< Raw_Int_Sample >
    unpack_rw(Raw_Samples_Pack const & rsp)
    {
        return rw_coder().decode< Raw_Int_Sample >(rsp.signal, rsp.signal_params);
    }
    static EventDetection_Events_Pack
    pack_ed(std::vector< EventDetection_Event > const & ed,
            EventDetection_Events_Params const & ed_params)
    {
        EventDetection_Events_Pack ed_pack;
        std::vector< long long > skip;
        std::vector< long long > len;
        long long last_end = ed_params.start_time;
        for (unsigned i = 0; i < ed.size(); ++i)
        {
            skip.push_back(ed[i].start - last_end);
            len.push_back(ed[i].length);
            last_end = ed[i].start + ed[i].length;
        }
        std::tie(ed_pack.skip, ed_pack.skip_params) = ed_skip_coder().encode(skip, false);
        std::tie(ed_pack.len, ed_pack.len_params) = ed_len_coder().encode(len, false);
        return ed_pack;
    }
    static std::vector< EventDetection_Event >
    unpack_ed(EventDetection_Events_Pack const & ed_pack,
              EventDetection_Events_Params const & ed_params,
              std::vector< Raw_Sample > const & rs,
              Raw_Samples_Params const & rs_params)
    {
        auto skip = ed_skip_coder().decode< long long >(ed_pack.skip, ed_pack.skip_params);
        auto len = ed_len_coder().decode< long long >(ed_pack.len, ed_pack.len_params);
        if (skip.size() != len.size())
        {
            throw std::runtime_error("unpack_ed failure: skip and length of different size");
        }
        std::vector< EventDetection_Event > ed(skip.size());
        long long last_end = ed_params.start_time;
        long long off_by_one = ed_params.start_time == rs_params.start_time; // hack
        for (unsigned i = 0; i < skip.size(); ++i)
        {
            ed[i].start = last_end + skip[i];
            ed[i].length = len[i];
            last_end = ed[i].start + ed[i].length;
            // use rs to reconstruct mean and stdv
            long long rs_start_idx = ed[i].start - rs_params.start_time + off_by_one;
            if (rs_start_idx < 0 or rs_start_idx + ed[i].length > (long long)rs.size() + off_by_one)
            {
                throw std::runtime_error("unpack_ed failure: bad rs_start_idx");
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
    static Basecall_Fastq_Pack
    pack_fq(std::string const & fq, unsigned qv_bits = 5)
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
        std::tie(fq_pack.bp, fq_pack.bp_params) = fq_bp_coder().encode(bp, false);
        std::tie(fq_pack.qv, fq_pack.qv_params) = fq_qv_coder().encode(qv, false);
        return fq_pack;
    }
    static std::string
    unpack_fq(Basecall_Fastq_Pack const & fq_pack)
    {
        std::string res;
        res += "@";
        res += fq_pack.read_name;
        res += "\n";
        auto bp = fq_bp_coder().decode< std::int8_t >(fq_pack.bp, fq_pack.bp_params);
        for (auto c : bp) res += c;
        res += "\n+\n";
        auto qv = fq_qv_coder().decode< std::uint8_t >(fq_pack.qv, fq_pack.qv_params);
        for (auto c : qv) res += (char)33 + c;
        res += "\n";
        return res;
    }
    static Basecall_Events_Pack
    pack_ev(std::vector< Basecall_Event > const & ev,
            std::string const & fq,
            Basecall_Events_Params const & ev_params,
            std::vector< EventDetection_Event > const & ed,
            std::string const & ed_gr,
            double sampling_rate,
            unsigned p_model_state_bits)
    {
        Basecall_Events_Pack ev_pack;
        ev_pack.ed_gr = ed_gr;
        ev_pack.start_time = time_to_int(ev_params.start_time, sampling_rate);
        ev_pack.duration = time_to_int(ev_params.duration, sampling_rate);
        ev_pack.state_size = ev[0].get_model_state().size();
        ev_pack.p_model_state_bits = p_model_state_bits;
        auto fqa = split_fq(fq);
        std::vector< long long > skip;
        std::vector< std::uint8_t > mv;
        std::vector< std::uint16_t > p_model_state;
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
                throw std::runtime_error("pack_ev failed: no matching ed event");
            }
            skip.push_back(j - last_j - 1);
            // move
            if (ev[i].move < 0 or ev[i].move > std::numeric_limits< uint8_t >::max())
            {
                throw std::runtime_error("pack_ev failed: invalid move");
            }
            mv.push_back(ev[i].move);
            // state
            auto s = ev[i].get_model_state();
            if (s.size() != ev_pack.state_size)
            {
                throw std::runtime_error("pack_ev failed: different state sizes");
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
            // p_model_state
            std::uint16_t p_model_state_val = ev[i].p_model_state * (1u << p_model_state_bits);
            if (p_model_state_val >= (1u << p_model_state_bits)) p_model_state_val = (1u << p_model_state_bits) - 1;
            p_model_state.push_back(p_model_state_val);
        }
        if (bases != fqa[1])
        {
            throw std::runtime_error("pack_ev failed: sequence of states does not match fastq seq");
        }

        std::tie(ev_pack.skip, ev_pack.skip_params) = ev_skip_coder().encode(skip, false);
        std::tie(ev_pack.move, ev_pack.move_params) = ev_move_coder().encode(mv, false);
        std::tie(ev_pack.p_model_state, ev_pack.p_model_state_params) = bit_packer().encode(p_model_state, p_model_state_bits);
        return ev_pack;
    } // pack_ev()
    static std::vector< Basecall_Event >
    unpack_ev(Basecall_Events_Pack const & ev_pack,
              std::string const & fq,
              std::vector< EventDetection_Event > const & ed,
              double sampling_rate)
    {
        std::vector< Basecall_Event > res;
        auto skip = ev_skip_coder().decode< long long >(ev_pack.skip, ev_pack.skip_params);
        auto mv = ev_move_coder().decode< std::uint8_t >(ev_pack.move, ev_pack.move_params);
        auto p_model_state = bit_packer().decode< std::uint16_t >(ev_pack.p_model_state, ev_pack.p_model_state_params);
        auto fqa = split_fq(fq);
        auto const & bases = fqa[1];
        if (skip.size() != mv.size() or p_model_state.size() != mv.size())
        {
            throw std::runtime_error("unpack_ev failed: skip and move have different sizes");
        }
        res.resize(skip.size());
        long long j = 0;
        std::string s;
        unsigned bases_pos = 0;
        unsigned p_model_state_bits;
        std::istringstream(ev_pack.p_model_state_params.at("num_bits")) >> p_model_state_bits;
        long long unsigned max_p_model_state_int = 1llu << p_model_state_bits;
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
            res[i].p_model_state = (double)p_model_state[i] / max_p_model_state_int;
        }
        return res;
    } // unpack_ev()
    static Basecall_Alignment_Pack
    pack_al(std::vector< Basecall_Alignment_Entry > const & al,
            std::string const & seq)
    {
        Basecall_Alignment_Pack al_pack;
        std::array< std::vector< uint8_t > , 2 > step_v;
        std::vector< int8_t > mv;
        step_v[0].reserve(al.size());
        step_v[1].reserve(al.size());
        mv.reserve(al.size());
        std::array< int, 2 > start_index = {{ -1, -1 }};
        std::array< int, 2 > next_index = {{ -1, -1 }};
        std::array< int, 2 > delta = {{ 1, -1 }};
        auto get_idx = [&] (unsigned i, unsigned k) {
            return k == 0? al[i].template_index : al[i].complement_index;
        };
        unsigned pos = 0;
        for (unsigned i = 0; i < al.size(); ++i)
        {
            for (unsigned k = 0; k < 2; ++k)
            {
                auto idx = get_idx(i, k);
                if (idx >= 0)
                {
                    if (start_index[k] < 0)
                    {
                        start_index[k] = idx;
                        next_index[k] = idx;
                    }
                    if (idx != next_index[k])
                    {
                        throw std::runtime_error("pack_al failed: unexpected index");
                    }
                    step_v[k].push_back(1);
                    next_index[k] += delta[k];
                }
                else // idx < 0
                {
                    step_v[k].push_back(0);
                }
            }
            // compute move
            auto kmer = al[i].get_kmer();
            size_t next_pos = seq.find(kmer, pos);
            if (next_pos == std::string::npos)
            {
                throw std::runtime_error("pack_al failed: cannot find kmer in 2d seq");
            }
            if (next_pos - pos > std::numeric_limits< int8_t >::max())
            {
                throw std::runtime_error("pack_al failed: move too large");
            }
            mv.push_back(next_pos - pos);
            pos = next_pos;
        }
        if (start_index[0] < 0)
        {
            throw std::runtime_error("pack_al failed: no template events");
        }
        if (start_index[1] < 0)
        {
            throw std::runtime_error("pack_al failed: no complement events");
        }
        al_pack.template_index_start = start_index[0];
        al_pack.complement_index_start = start_index[1];
        al_pack.kmer_size = al[0].get_kmer().size();
        std::tie(al_pack.template_step, al_pack.template_step_params) = bit_packer().encode(step_v[0], 1);
        std::tie(al_pack.complement_step, al_pack.complement_step_params) = bit_packer().encode(step_v[1], 1);
        std::tie(al_pack.move, al_pack.move_params) = ev_move_coder().encode(mv, false);
        return al_pack;
    } // pack_al()
    static std::vector< Basecall_Alignment_Entry >
    unpack_al(Basecall_Alignment_Pack const & al_pack,
              std::string const & seq)
    {
        std::vector< Basecall_Alignment_Entry > al;
        std::array< std::vector< uint8_t >, 2 > step_v =
            {{ bit_packer().decode< uint8_t >(al_pack.template_step, al_pack.template_step_params),
               bit_packer().decode< uint8_t >(al_pack.complement_step, al_pack.complement_step_params) }};
        auto mv = ev_move_coder().decode< int8_t >(al_pack.move, al_pack.move_params);
        if (step_v[1].size() != step_v[0].size()
            or mv.size() != step_v[0].size())
        {
            throw std::runtime_error("unpack_al failed: mismatching size");
        }
        al.resize(step_v[0].size());
        std::array< unsigned, 2 > crt_index = {{ al_pack.template_index_start, al_pack.complement_index_start }};
        std::array< int, 2 > delta = {{ 1, -1 }};
        auto pos = 0;
        auto set_idx = [&] (unsigned i, unsigned k, int val) {
            if (k == 0)
            {
                al[i].template_index = val;
            }
            else
            {
                al[i].complement_index = val;
            }
        };
        for (unsigned i = 0; i < step_v[0].size(); ++i)
        {
            for (unsigned k = 0; k < 2; ++k)
            {
                if (step_v[k][i] > 0)
                {
                    set_idx(i, k, crt_index[k]);
                    crt_index[k] += delta[k];
                }
                else
                {
                    set_idx(i, k, -1);
                }
            }
            // set kmer
            pos += mv[i];
            std::copy(seq.begin() + pos, seq.begin() + pos + al_pack.kmer_size, al[i].kmer.begin());
            if (al_pack.kmer_size < MAX_K_LEN) al[i].kmer[al_pack.kmer_size] = 0;
        }
        return al;
    } // unpack_al()

    //
    // Static helpers
    //
    static long long
    time_to_int(double tf, double sampling_rate)
    {
        return tf * sampling_rate + .5;
    }
    static double
    time_to_float(long long ti, double sampling_rate)
    {
        return (long double)ti / sampling_rate;
    }
    static float
    raw_int_sample_to_float(int si, Channel_Id_Params const & channel_id_params)
    {
        return ((float)si + channel_id_params.offset)
            * channel_id_params.range / channel_id_params.digitisation;
    }
    static std::string
    fq2seq(std::string const & fq)
    {
        return split_fq(fq)[1];
    }
    static std::array< std::string, 4 >
    split_fq(std::string const & fq)
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

private:
    //
    // Cached file data
    //
    Channel_Id_Params _channel_id_params;
    std::vector< std::string > _raw_samples_read_names;
    std::vector< std::string > _eventdetection_groups;
    std::map< std::string, std::vector< std::string > > _eventdetection_read_names;
    std::vector< std::string > _basecall_groups;
    std::array< std::vector< std::string >, 3 > _basecall_strand_groups;
    std::map< std::string, std::string > _basecall_1d_group;
    std::map< std::string, std::string > _basecall_eventdetection_group;

    //
    // Cache updaters
    //
    void
    reload()
    {
        load_channel_id_params();
        load_raw_samples_read_names();
        load_eventdetection_groups();
        load_basecall_groups();
    }
    void
    load_channel_id_params()
    {
        if (not Base::group_exists(channel_id_path())) return;
        Base::read(channel_id_path() + "/channel_number", _channel_id_params.channel_number);
        Base::read(channel_id_path() + "/digitisation", _channel_id_params.digitisation);
        Base::read(channel_id_path() + "/offset", _channel_id_params.offset);
        Base::read(channel_id_path() + "/range", _channel_id_params.range);
        Base::read(channel_id_path() + "/sampling_rate", _channel_id_params.sampling_rate);
    }
    void
    load_raw_samples_read_names()
    {
        _raw_samples_read_names.clear();
        if (not Base::group_exists(raw_samples_root_path())) return;
        auto rn_l = Base::list_group(raw_samples_root_path());
        for (auto const & rn : rn_l)
        {
            if (have_raw_samples_unpack(rn)
                or have_raw_samples_pack(rn))
            {
                _raw_samples_read_names.push_back(rn);
            }
        }
    }
    void
    load_eventdetection_groups()
    {
        _eventdetection_groups.clear();
        _eventdetection_read_names.clear();
        if (not Base::group_exists(eventdetection_root_path())) return;
        auto ed_gr_prefix = eventdetection_group_prefix();
        auto gr_l = Base::list_group(eventdetection_root_path());
        for (auto const & g : gr_l)
        {
            if (g.substr(0, ed_gr_prefix.size()) != ed_gr_prefix) continue;
            std::string gr = g.substr(ed_gr_prefix.size());
            _eventdetection_groups.push_back(gr);
            _eventdetection_read_names[gr] = detect_eventdetection_read_names(gr);
        }
    }
    std::vector< std::string >
    detect_eventdetection_read_names(std::string const & gr) const
    {
        std::vector< std::string > res;
        std::string p = eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr + "/Reads";
        if (not Base::group_exists(p)) return res;
        auto rn_l = Base::list_group(p);
        for (auto const & rn : rn_l)
        {
            if (have_eventdetection_events_unpack(gr, rn)
                or have_eventdetection_events_pack(gr, rn))
            {
                res.push_back(rn);
            }
        }
        return res;
    }
    void
    load_basecall_groups()
    {
        _basecall_groups.clear();
        std::for_each(
            _basecall_strand_groups.begin(), _basecall_strand_groups.end(),
            [] (decltype(_basecall_strand_groups)::value_type & v) {
                v.clear();
            });
        _basecall_1d_group.clear();
        _basecall_eventdetection_group.clear();
        if (not Base::group_exists(basecall_root_path())) return;
        auto bc_gr_prefix = basecall_group_prefix();
        auto gr_l = Base::list_group(basecall_root_path());
        for (auto const & g : gr_l)
        {
            if (g.substr(0, bc_gr_prefix.size()) != bc_gr_prefix) continue;
            std::string gr = g.substr(bc_gr_prefix.size());
            _basecall_groups.push_back(gr);
            bool have_1d_subgroups = false;
            for (unsigned st = 0; st < 3; ++st)
            {
                if (Base::group_exists(basecall_strand_group_path(gr, st)))
                {
                    _basecall_strand_groups[st].push_back(gr);
                    if (st < 2)
                    {
                        have_1d_subgroups = true;
                        _basecall_eventdetection_group[gr] = detect_basecall_eventdetection_group(gr);
                    }
                }
            }
            _basecall_1d_group[gr] = (have_1d_subgroups
                                      ? gr
                                      : detect_basecall_1d_group(gr));
        }
    }
    std::string
    detect_basecall_1d_group(std::string const & gr) const
    {
        std::string path = basecall_group_path(gr) + "/basecall_1d";
        if (Base::attribute_exists(path))
        {
            std::string tmp;
            Base::read(path, tmp);
            auto pref = basecall_root_path().substr(1) + "/" + basecall_group_prefix();
            if (tmp.size() >= pref.size()
                and tmp.substr(0, pref.size()) == pref)
            {
                auto gr_1d = tmp.substr(pref.size());
                if (have_basecall_group(gr_1d))
                {
                    return gr_1d;
                }
            }
        }
        return gr;
    }
    std::string
    detect_basecall_eventdetection_group(std::string const & gr) const
    {
        auto bc_params = get_basecall_params(gr);
        if (bc_params.count("event_detection"))
        {
            auto && tmp = bc_params.at("event_detection");
            auto pref = eventdetection_root_path().substr(1) + "/" + eventdetection_group_prefix();
            if (tmp.size() >= pref.size()
                and tmp.substr(0, pref.size()) == pref)
            {
                auto ed_gr = tmp.substr(pref.size());
                if (have_eventdetection_group(ed_gr))
                {
                    return ed_gr;
                }
            }
        }
        return "";
    }

    //
    // Functions that fill in empty arguments with default values
    //
    std::string const &
    fill_raw_samples_read_name(std::string const & rn) const
    {
        return (not rn.empty() or _raw_samples_read_names.empty()
                ? rn
                : _raw_samples_read_names.front());
    }
    std::string const &
    fill_eventdetection_group(std::string const & gr) const
    {
        return (not gr.empty() or _eventdetection_groups.empty()
                ? gr
                : _eventdetection_groups.front());
    }
    std::string const &
    fill_eventdetection_read_name(std::string const & gr, std::string const & rn) const
    {
        return (not rn.empty()
                or _eventdetection_read_names.count(gr) == 0
                or _eventdetection_read_names.at(gr).empty()
                ? rn
                : _eventdetection_read_names.at(gr).front());
    }
    std::string const &
    fill_basecall_group(unsigned st, std::string const & gr) const
    {
        return (not gr.empty()
                or _basecall_strand_groups.at(st).empty()
                ? gr
                : _basecall_strand_groups.at(st).front());
    }
    std::string const &
    fill_basecall_1d_group(unsigned st, std::string const & gr) const
    {
        auto && _gr = fill_basecall_group(st, gr);
        return get_basecall_1d_group(_gr);
    }

    //
    // Fast5 internal paths
    //
    static std::string file_version_path() { return "/file_version"; }
    static std::string channel_id_path()   { return "/UniqueGlobalKey/channel_id"; }
    static std::string tracking_id_path()  { return "/UniqueGlobalKey/tracking_id"; }
    static std::string sequences_path()    { return "/Sequences/Meta"; }
    static std::string raw_samples_root_path() { return "/Raw/Reads"; }
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
        return raw_samples_path(rn) + "_Pack";
    }
    static std::string eventdetection_root_path() { return "/Analyses"; }
    static std::string eventdetection_group_prefix() { return "EventDetection_"; }
    static std::string eventdetection_group_path(std::string const & gr)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + gr;
    }
    static std::string eventdetection_events_params_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_group_path(gr) + "/Reads/" + rn;
    }
    static std::string eventdetection_events_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_group_path(gr) + "/Reads/" + rn + "/Events";
    }
    static std::string eventdetection_events_pack_path(std::string const & gr, std::string const & rn)
    {
        return eventdetection_events_path(gr, rn) + "_Pack";
    }
    static std::string basecall_root_path() { return "/Analyses"; }
    static std::string basecall_group_prefix() { return "Basecall_"; }
    static std::string strand_name(unsigned st)
    {
        static const std::array< std::string, 3 > _strand_name =
            {{ "template", "complement", "2D" }};
        return _strand_name.at(st);
    }
    static std::string basecall_strand_subgroup(unsigned st)
    {
        return std::string("BaseCalled_") + strand_name(st);
    }
    static std::string basecall_group_path(std::string const & gr)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + gr;
    }
    static std::string basecall_strand_group_path(std::string const & gr, unsigned st)
    {
        return basecall_group_path(gr) + "/" + basecall_strand_subgroup(st);
    }
    static std::string basecall_log_path(std::string const & gr)
    {
        return basecall_group_path(gr) + "/Log";
    }
    static std::string basecall_fastq_path(std::string const & gr, unsigned st)
    {
        return basecall_strand_group_path(gr, st) + "/Fastq";
    }
    static std::string basecall_fastq_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_fastq_path(gr, st) + "_Pack";
    }
    static std::string basecall_model_path(std::string const & gr, unsigned st)
    {
        return basecall_strand_group_path(gr, st) + "/Model";
    }
    static std::string basecall_model_file_path(std::string const & gr, unsigned st)
    {
        return basecall_group_path(gr) + "/Summary/basecall_1d_" + strand_name(st) + "/model_file";
    }
    static std::string basecall_events_path(std::string const & gr, unsigned st)
    {
        return basecall_strand_group_path(gr, st) + "/Events";
    }
    static std::string basecall_events_pack_path(std::string const & gr, unsigned st)
    {
        return basecall_events_path(gr, st) + "_Pack";
    }
    static std::string basecall_alignment_path(std::string const & gr)
    {
        return basecall_strand_group_path(gr, 2) + "/Alignment";
    }
    static std::string basecall_alignment_pack_path(std::string const & gr)
    {
        return basecall_alignment_path(gr) + "_Pack";
    }

    //
    // Packers
    //
    static Huffman_Packer const & rw_coder()      { return Huffman_Packer::get_coder("fast5_rw_1"); }
    static Huffman_Packer const & ed_skip_coder() { return Huffman_Packer::get_coder("fast5_ed_skip_1"); }
    static Huffman_Packer const & ed_len_coder()  { return Huffman_Packer::get_coder("fast5_ed_len_1"); }
    static Huffman_Packer const & fq_bp_coder()   { return Huffman_Packer::get_coder("fast5_fq_bp_1"); }
    static Huffman_Packer const & fq_qv_coder()   { return Huffman_Packer::get_coder("fast5_fq_qv_1"); }
    static Huffman_Packer const & ev_skip_coder() { return Huffman_Packer::get_coder("fast5_ev_skip_1"); }
    static Huffman_Packer const & ev_move_coder() { return Huffman_Packer::get_coder("fast5_ev_move_1"); }
    static Bit_Packer     const & bit_packer()    { return Bit_Packer::get_packer(); }
}; // class File

} // namespace fast5

#endif
