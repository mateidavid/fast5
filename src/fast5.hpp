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
    inline static std::string array_to_string(const std::array< char, MAX_K_LEN >& a)
    {
        return std::string(a.begin(), std::find(a.begin(), a.end(), '\0'));
    }
}

namespace fast5
{

struct Channel_Id_Parameters
{
    std::string channel_number;
    double digitisation;
    double offset;
    double range;
    double sampling_rate;
}; // struct Channel_Id_Parameters

typedef std::map< std::string, std::string > Tracking_Id_Parameters;

typedef std::map< std::string, std::string > Sequences_Parameters;

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
    typedef std::vector< std::uint8_t > signal_type;
    typedef std::map< std::string, std::string > signal_param_type;
    signal_type signal;
    signal_param_type signal_param;
}; // struct Raw_Samples_Pack

struct EventDetection_Event_Entry
{
    double mean;
    double stdv;
    long long start;
    long long length;
    friend bool operator == (const EventDetection_Event_Entry& lhs, const EventDetection_Event_Entry& rhs)
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
    typedef std::vector< std::uint8_t > skip_type;
    typedef std::map< std::string, std::string > skip_param_type;
    skip_type skip;
    skip_param_type skip_param;
    typedef std::vector< std::uint8_t > len_type;
    typedef std::map< std::string, std::string > len_param_type;
    len_type len;
    len_param_type len_param;
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
    friend bool operator == (const Model_Entry& lhs, const Model_Entry& rhs)
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
    friend bool operator == (const Event_Entry& lhs, const Event_Entry& rhs)
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
    friend bool operator == (const Event_Alignment_Entry& lhs, const Event_Alignment_Entry& rhs)
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
    File(const std::string& file_name, bool rw = false) { open(file_name, rw); }

    bool is_open() const { return static_cast< const Base* >(this)->is_open(); }
    bool is_rw() const { return static_cast< const Base* >(this)->is_rw(); }
    const std::string& file_name() const { return static_cast< const Base* >(this)->file_name(); }
    void create(const std::string& file_name, bool truncate = false) { static_cast< Base* >(this)->create(file_name, truncate); }
    void close() { static_cast< Base* >(this)->close(); }

    void open(const std::string& file_name, bool rw = false)
    {
        Base::open(file_name, rw);
        if (is_open())
        {
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
        return Base::group_exists(channel_id_path());
    }
    /**
     * Extract "/UniqueGlobalKey/channel_id" attributes.
     */
    Channel_Id_Parameters get_channel_id_params() const
    {
        Channel_Id_Parameters res;
        Base::read(channel_id_path() + "/channel_number", res.channel_number);
        Base::read(channel_id_path() + "/digitisation", res.digitisation);
        Base::read(channel_id_path() + "/offset", res.offset);
        Base::read(channel_id_path() + "/range", res.range);
        Base::read(channel_id_path() + "/sampling_rate", res.sampling_rate);
        return res;
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
        auto channel_id_params = get_channel_id_params();
        return channel_id_params.sampling_rate;
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
    const std::vector< std::string >& get_raw_samples_read_name_list() const
    {
        return _raw_samples_read_name_list;
    }
    /**
     * Check if raw samples exist.
     * If _rn non-empty, check if raw samples exist for given read.
     */
    bool have_raw_samples(const std::string& _rn = std::string()) const
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
    Raw_Samples_Parameters get_raw_samples_params(const std::string& _rn = std::string()) const
    {
        Raw_Samples_Parameters res;
        const std::string& rn = not _rn.empty()? _rn : get_raw_samples_read_name_list().front();
        std::string p = raw_samples_params_path(rn);
        Base::read(p + "/read_id", res.read_id);
        Base::read(p + "/read_number", res.read_number);
        Base::read(p + "/start_mux", res.start_mux);
        Base::read(p + "/start_time", res.start_time);
        Base::read(p + "/duration", res.duration);
        return res;
    }
    /**
     * Get raw samples for given read name as ints (default: first read name).
     */
    std::vector< Raw_Samples_Int_Entry > get_raw_samples_int(const std::string& _rn = std::string()) const
    {
        // get raw samples
        std::vector< Raw_Samples_Int_Entry > res;
        const std::string& rn = not _rn.empty()? _rn : get_raw_samples_read_name_list().front();
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
    std::vector< Raw_Samples_Entry > get_raw_samples(const std::string& _rn = std::string()) const
    {
        // get raw samples
        auto raw_samples_int = get_raw_samples_int(_rn);
        // get scaling parameters
        auto channel_id_params = get_channel_id_params();
        // decode levels
        std::vector< Raw_Samples_Entry > res;
        res.reserve(raw_samples_int.size());
        for (auto int_level : raw_samples_int)
        {
            res.push_back((static_cast< float >(int_level) + channel_id_params.offset)
                          * channel_id_params.range / channel_id_params.digitisation);
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
    const std::vector< std::string >& get_eventdetection_group_list() const
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
    std::vector< std::string > get_eventdetection_read_name_list(const std::string& _ed_gr = std::string()) const
    {
        const std::string& ed_gr = not _ed_gr.empty()? _ed_gr : get_eventdetection_group_list().front();
        return detect_eventdetection_read_name_list(ed_gr);
    }
    /**
     * Check if EventDetection events exist.
     * If _ed_gr given: check if events exist for given group; else: check first EventDetection group.
     * If _rn given: check if events exist for given group and read name.
     */
    bool have_eventdetection_events(
        const std::string& _ed_gr = std::string(),
        const std::string& _rn = std::string()) const
    {
        std::string ed_gr;
        if (_ed_gr.empty())
        {
            auto ed_gr_l = get_eventdetection_group_list();
            if (ed_gr_l.empty()) return false;
            ed_gr = ed_gr_l.front();
        }
        else
        {
            ed_gr = _ed_gr;
        }
        auto rn_l = get_eventdetection_read_name_list(ed_gr);
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
    bool have_eventdetection_events_unpack(std::string const & ed_gr, std::string const & rn) const
    {
        return Base::dataset_exists(eventdetection_events_path(ed_gr, rn));
    }
    bool have_eventdetection_events_pack(std::string const & ed_gr, std::string const & rn) const
    {
        return Base::group_exists(eventdetection_events_pack_path(ed_gr, rn));
    }
    /**
     * Get EventDetection params for given EventDetection group (default: first EventDetection group).
     */
    std::map< std::string, std::string > get_eventdetection_params(const std::string& _ed_gr = std::string()) const
    {
        const std::string& ed_gr = not _ed_gr.empty()? _ed_gr : get_eventdetection_group_list().front();
        return get_attr_map(eventdetection_params_path(ed_gr));
    }
    /**
     * Get EventDetection event params for given EventDetection group, and given read name
     * (default: first EventDetection group, and first read name in it).
     */
    EventDetection_Event_Parameters get_eventdetection_event_params(
        const std::string& _ed_gr = std::string(), const std::string& _rn = std::string()) const
    {
        EventDetection_Event_Parameters res;
        const std::string& ed_gr = not _ed_gr.empty()? _ed_gr : get_eventdetection_group_list().front();
        const std::string rn = not _rn.empty()? _rn : get_eventdetection_read_name_list(ed_gr).front();
        auto p = eventdetection_event_params_path(ed_gr, rn);
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
            res.abasic_found = 0;
        }
        return res;
    }
    /**
     * Get EventDetection events for given EventDetection group, and given read name.
     */
    std::vector< EventDetection_Event_Entry > get_eventdetection_events(
        const std::string& _ed_gr = std::string(), const std::string& _rn = std::string()) const
    {
        std::vector< EventDetection_Event_Entry > res;
        const std::string& ed_gr = not _ed_gr.empty()? _ed_gr : get_eventdetection_group_list().front();
        const std::string rn = not _rn.empty()? _rn : get_eventdetection_read_name_list(ed_gr).front();
        if (have_eventdetection_events_unpack(ed_gr, rn))
        {
            auto p = eventdetection_events_path(ed_gr, rn);
            auto struct_member_names = Base::get_struct_members(p);
            assert(struct_member_names.size() >= 4);
            bool have_stdv = false;
            bool have_variance = false;
            for (const auto& s : struct_member_names)
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
            auto ed_pack = get_eventdetection_events_pack(ed_gr, rn);
            auto ed_param = get_eventdetection_event_params(ed_gr, rn);
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
    const std::vector< std::string >& get_basecall_group_list() const
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
    const std::vector< std::string >& get_basecall_strand_group_list(unsigned st) const
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
    std::map< std::string, std::string > get_basecall_params(const std::string& bc_gr) const
    {
        return get_attr_map(basecall_root_path() + "/" + basecall_group_prefix() + bc_gr);
    }
    /**
     * Check if Basecall log exists for given Basecall group.
     */
    bool have_basecall_log(const std::string& bc_gr) const
    {
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/Log";
        return Base::exists(path);
    }
    /**
     * Get Basecall log for given Basecall group.
     */
    std::string get_basecall_log(const std::string& bc_gr) const
    {
        std::string res;
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/Log";
        Base::read(path, res);
        return res;
    }
    /**
     * Check if Basecall fastq exists for given Basecall group and given strand.
     */
    bool have_basecall_fastq(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        if (_bc_gr.empty() and get_basecall_strand_group_list(st).empty()) return false;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        return Base::dataset_exists(basecall_fastq_path(bc_gr, st));
    }
    /**
     * Get Basecall fastq for given Basecall group and given strand.
     */
    std::string get_basecall_fastq(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        std::string res;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        Base::read(basecall_fastq_path(bc_gr, st), res);
        return res;
    }
    /**
     * Add Basecall fastq
     */
    void add_basecall_fastq(unsigned st, const std::string& bc_gr, const std::string& fq) const
    {
        Base::write(basecall_fastq_path(bc_gr, st), true, fq);
    }
    /**
     * Check if Basecall seq exists for given Basecall group and given strand.
     */
    bool have_basecall_seq(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        return have_basecall_fastq(st, _bc_gr);
    }
    /**
     * Get Basecall sequence for given Basecall group and given strand.
     */
    std::string get_basecall_seq(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        return fq2seq(get_basecall_fastq(st, _bc_gr));
    }
    /**
     * Add Basecall seq
     */
    void add_basecall_seq(unsigned st, const std::string& bc_gr,
                          const std::string& name, const std::string& seq, int default_qual = 33) const
    {
        std::ostringstream oss;
        oss << '@' << name << std::endl
            << seq << std::endl
            << '+' << std::endl
            << std::string(seq.size(), static_cast< char >(default_qual));
        add_basecall_fastq(st, bc_gr, oss.str());
    }
    /**
     * Check if Basecall model exist for given Basecall group and given strand.
     */
    bool have_basecall_model(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        if (_bc_gr.empty() and get_basecall_strand_group_list(st).empty()) return false;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        return Base::dataset_exists(basecall_model_path(bc_gr_1d, st));
    }
    /**
     * Get Basecall model file name for given Basecall group and given strand.
     */
    std::string get_basecall_model_file(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        std::string res;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        assert(Base::exists(basecall_model_file_path(bc_gr_1d, st)));
        Base::read(basecall_model_file_path(bc_gr_1d, st), res);
        return res;
    }
    void add_basecall_model_file(unsigned st, const std::string& bc_gr, const std::string& file_name) const
    {
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        std::string path = basecall_model_file_path(bc_gr_1d, st);
        Base::write(path, false, file_name);
    }
    /**
     * Get Basecall model parameters for given Basecall group and given strand.
     */
    Model_Parameters get_basecall_model_params(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        Model_Parameters res;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        std::string path = basecall_model_path(bc_gr_1d, st);
        Base::read(path + "/scale", res.scale);
        Base::read(path + "/shift", res.shift);
        Base::read(path + "/drift", res.drift);
        Base::read(path + "/var", res.var);
        Base::read(path + "/scale_sd", res.scale_sd);
        Base::read(path + "/var_sd", res.var_sd);
        return res;
    }
    template < typename T >
    void add_basecall_model_params(unsigned st, const std::string& bc_gr, const T& params) const
    {
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        std::string path = basecall_model_path(bc_gr_1d, st);
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
    std::vector< Model_Entry > get_basecall_model(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        std::vector< Model_Entry > res;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        hdf5_tools::Compound_Map m;
        m.add_member("kmer", &Model_Entry::kmer);
        m.add_member("level_mean", &Model_Entry::level_mean);
        m.add_member("level_stdv", &Model_Entry::level_stdv);
        m.add_member("sd_mean", &Model_Entry::sd_mean);
        m.add_member("sd_stdv", &Model_Entry::sd_stdv);
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        Base::read(basecall_model_path(bc_gr_1d, st), res, m);
        return res;
    }
    /**
     * Add Basecall model
     */
    template < typename T >
    void add_basecall_model(unsigned st, const std::string& bc_gr, const std::vector< T >& m) const
    {
        hdf5_tools::Compound_Map cm;
        cm.add_member("kmer", &T::kmer);
        cm.add_member("level_mean", &T::level_mean);
        cm.add_member("level_stdv", &T::level_stdv);
        cm.add_member("sd_mean", &T::sd_mean);
        cm.add_member("sd_stdv", &T::sd_stdv);
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        Base::write(basecall_model_path(bc_gr_1d, st), true, m, cm);
    }
    /**
     * Check if Basecall events exist for given Basecall group and given strand.
     */
    bool have_basecall_events(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        if (_bc_gr.empty() and get_basecall_strand_group_list(st).empty()) return false;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        return Base::dataset_exists(basecall_events_path(bc_gr_1d, st));
    }
    /**
     * Get Basecall events for given Basecall group and given strand.
     */
    std::vector< Event_Entry > get_basecall_events(unsigned st, const std::string& _bc_gr = std::string()) const
    {
        std::vector< Event_Entry > res;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(st).front();
        hdf5_tools::Compound_Map m;
        m.add_member("mean", &Event_Entry::mean);
        m.add_member("start", &Event_Entry::start);
        m.add_member("stdv", &Event_Entry::stdv);
        m.add_member("length", &Event_Entry::length);
        m.add_member("p_model_state", &Event_Entry::p_model_state);
        m.add_member("model_state", &Event_Entry::model_state);
        m.add_member("move", &Event_Entry::move);
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        Base::read(basecall_events_path(bc_gr_1d, st), res, m);
        return res;
    }
    /**
     * Add Basecall events
     */
    template < typename T >
    void add_basecall_events(unsigned st, const std::string& bc_gr, const std::vector< T >& ev) const
    {
        hdf5_tools::Compound_Map cm;
        cm.add_member("mean", &T::mean);
        cm.add_member("start", &T::start);
        cm.add_member("stdv", &T::stdv);
        cm.add_member("length", &T::length);
        cm.add_member("p_model_state", &T::p_model_state);
        cm.add_member("model_state", &T::model_state);
        cm.add_member("move", &T::move);
        auto bc_gr_1d = get_basecall_group_1d(bc_gr);
        Base::write(basecall_events_path(bc_gr_1d, st), true, ev, cm);
    }
    /**
     * Check if Basecall event alignment exist for given Basecall group.
     */
    bool have_basecall_event_alignment(const std::string& _bc_gr = std::string()) const
    {
        if (_bc_gr.empty() and get_basecall_strand_group_list(2).empty()) return false;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(2).front();
        return Base::dataset_exists(basecall_event_alignment_path(bc_gr));
    }
    /**
     * Get Basecall events for given Basecall group.
     */
    std::vector< Event_Alignment_Entry > get_basecall_event_alignment(const std::string& _bc_gr = std::string()) const
    {
        std::vector< Event_Alignment_Entry > res;
        const std::string& bc_gr = not _bc_gr.empty()? _bc_gr : get_basecall_strand_group_list(2).front();
        hdf5_tools::Compound_Map m;
        m.add_member("template", &Event_Alignment_Entry::template_index);
        m.add_member("complement", &Event_Alignment_Entry::complement_index);
        m.add_member("kmer", &Event_Alignment_Entry::kmer);
        Base::read(basecall_event_alignment_path(bc_gr), res, m);
        return res;
    }

    /**
     * Get basecall group holding 1d calls.
     */
    std::string get_basecall_group_1d(const std::string& bc_gr) const
    {
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/basecall_1d";
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
        return bc_gr;
    }
    /**
     * Get EventDetection group for given Basecall group, if available.
     */
    std::string get_basecall_eventdetection_group(const std::string& bc_gr) const
    {
        std::string path = basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/event_detection";
        if (Base::attribute_exists(path))
        {
            std::string tmp;
            Base::read(path, tmp);
            auto pos = tmp.find(eventdetection_group_prefix());
            if (pos != std::string::npos)
            {
                pos += eventdetection_group_prefix().size();
                auto end_pos = tmp.find("/", pos);
                if (end_pos == std::string::npos)
                {
                    end_pos = tmp.size();
                }
                return tmp.substr(pos, end_pos - pos);
            }
        }
        return std::string();
    }

    static std::string fq2seq(const std::string& fq)
    {
        return split_fq(fq)[1];
    }
    static std::array< std::string, 4 > split_fq(const std::string& fq)
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
    static void copy_attributes(File const & src_f, File const & dst_f)
    {
        assert(src_f.is_open());
        assert(dst_f.is_open());
        assert(dst_f.is_rw());
        rec_copy_attributes(src_f, dst_f, "");
    }

    /**
     * Pack raw samples
     */
    static Raw_Samples_Pack pack_rw(std::vector< Raw_Samples_Int_Entry > const & rsi)
    {
        Raw_Samples_Pack rsp;
        std::tie(rsp.signal, rsp.signal_param) = rw_coder().encode(rsi);
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
        std::tie(ed_pack.skip, ed_pack.skip_param) = ed_skip_coder().encode(skip);
        std::tie(ed_pack.len, ed_pack.len_param) = ed_len_coder().encode(len);
        return ed_pack;
    }
    /**
     * Unpack eventdetection events
     */
    static std::vector< EventDetection_Event_Entry > unpack_ed(
        EventDetection_Events_Pack const & ed_pack,
        EventDetection_Event_Parameters const & ed_par,
        std::vector< Raw_Samples_Entry > const & rs,
        Raw_Samples_Parameters const & rs_par)
    {
        auto skip = ed_skip_coder().decode< long long >(ed_pack.skip, ed_pack.skip_param);
        auto len = ed_len_coder().decode< long long >(ed_pack.len, ed_pack.len_param);
        if (skip.size() != len.size())
        {
            throw std::invalid_argument("unpack_ed failure: skip and length of different size");
        }
        std::vector< EventDetection_Event_Entry > ed(skip.size());
        long long last_end = ed_par.start_time;
        for (unsigned i = 0; i < skip.size(); ++i)
        {
            ed[i].start = last_end + skip[i];
            ed[i].length = len[i];
            last_end = ed[i].start + ed[i].length;
            // use rs to reconstruct mean and stdv
            long long rs_start_idx = ed[i].start - rs_par.start_time;
            if (rs_start_idx < 0 or rs_start_idx + ed[i].length > (long long)rs.size())
            {
                throw std::invalid_argument("unpack_ed failure: bad rs_start_idx");
            }
            double s = 0.0;
            double s2 = 0.0;
            unsigned n = ed[i].length;
            for (unsigned j = 0; j < n; ++j)
            {
                s += rs[i];
                s2 += rs[i] * rs[i];
            }
            ed[i].mean = s / n;
            ed[i].stdv = n > 1
                ? std::sqrt(s2/(n - 1) - s*s/(n * (n - 1)))
                : 0;
        }
        return ed;
    }

private:
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
        for (const auto& g : g_list)
        {
            if (g.size() <= eventdetection_group_prefix().size()) continue;
            auto p = std::mismatch(eventdetection_group_prefix().begin(),
                                   eventdetection_group_prefix().end(),
                                   g.begin());
            if (p.first != eventdetection_group_prefix().end()) continue;
            _eventdetection_group_list.emplace_back(p.second, g.end());
        }
    }

    std::vector< std::string > detect_eventdetection_read_name_list(const std::string& ed_gr) const
    {
        std::vector< std::string > res;
        std::string p = eventdetection_root_path() + "/" + eventdetection_group_prefix() + ed_gr + "/Reads";
        if (not Base::group_exists(p)) return res;
        auto rn_list = Base::list_group(p);
        for (const auto& rn : rn_list)
        {
            if (Base::dataset_exists(eventdetection_events_path(ed_gr, rn))
                or Base::group_exists(eventdetection_events_pack_path(ed_gr, rn)))
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
        for (const auto& g : g_list)
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

    std::map< std::string, std::string > get_attr_map(const std::string& path) const
    {
        std::map< std::string, std::string > res;
        auto a_list = Base::get_attr_list(path);
        for (const auto& a : a_list)
        {
            std::string tmp;
            Base::read(path + "/" + a, tmp);
            res[a] = tmp;
        }
        return res;
    }

    void add_attr_map(std::string const & path, std::map< std::string, std::string > const & attr_m) const
    {
        for (auto const & p : attr_m)
        {
            Base::write_attribute(path + "/" + p.first, p.second);
        }
    }

    // list of read names for which we have raw samples
    std::vector< std::string > _raw_samples_read_name_list;

    // list of EventDetection groups
    std::vector< std::string > _eventdetection_group_list;

    // list of Basecall groups
    std::vector< std::string > _basecall_group_list;

    // list of per-strand Basecall groups; 0/1/2 = template/complement/2d
    std::array< std::vector< std::string >, 3 > _basecall_strand_group_list;

    // static paths
    static const std::string& file_version_path()
    {
        static const std::string _file_version_path = "/file_version";
        return _file_version_path;
    }

    static const std::string& channel_id_path()
    {
        static const std::string _channel_id_path = "/UniqueGlobalKey/channel_id";
        return _channel_id_path;
    }
    static const std::string& tracking_id_path()
    {
        static const std::string _tracking_id_path = "/UniqueGlobalKey/tracking_id";
        return _tracking_id_path;
    }
    static const std::string& raw_samples_root_path()
    {
        static const std::string _raw_samples_root_path = "/Raw/Reads";
        return _raw_samples_root_path;
    }
    static std::string raw_samples_params_path(const std::string& rn)
    {
        return raw_samples_root_path() + "/" + rn;
    }
    static std::string raw_samples_path(const std::string& rn)
    {
        return raw_samples_root_path() + "/" + rn + "/Signal";
    }
    static std::string raw_samples_pack_path(const std::string& rn)
    {
        return raw_samples_root_path() + "/" + rn + "/Signal_Pack";
    }
    static const std::string& sequences_path()
    {
        static const std::string _sequences_path = "/Sequences/Meta";
        return _sequences_path;
    }
    static const std::string& eventdetection_root_path()
    {
        static const std::string _eventdetection_root_path = "/Analyses";
        return _eventdetection_root_path;
    }
    static const std::string& eventdetection_group_prefix()
    {
        static const std::string _eventdetection_group_prefix = "EventDetection_";
        return _eventdetection_group_prefix;
    }
    static std::string eventdetection_params_path(const std::string& ed_gr)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + ed_gr;
    }
    static std::string eventdetection_event_params_path(const std::string& ed_gr, const std::string& rn)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + ed_gr + "/Reads/" + rn;
    }
    static std::string eventdetection_events_path(const std::string& ed_gr, const std::string& rn)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + ed_gr + "/Reads/" + rn + "/Events";
    }
    static std::string eventdetection_events_pack_path(const std::string& ed_gr, const std::string& rn)
    {
        return eventdetection_root_path() + "/" + eventdetection_group_prefix() + ed_gr + "/Reads/" + rn + "/Events_Pack";
    }

    static const std::string& basecall_root_path()
    {
        static const std::string _basecall_root_path = "/Analyses";
        return _basecall_root_path;
    }
    static const std::string& basecall_group_prefix()
    {
        static const std::string _basecall_group_prefix = "Basecall_";
        return _basecall_group_prefix;
    }
    static const std::string& basecall_strand_subgroup(unsigned st)
    {
        static const std::array< std::string, 3 > _basecall_strand_subgroup =
            {{ "BaseCalled_template", "BaseCalled_complement", "BaseCalled_2D" }};
        return _basecall_strand_subgroup[st];
    }
    static std::string basecall_fastq_path(const std::string& bc_gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/"
            + basecall_strand_subgroup(st) + "/Fastq";
    }
    static std::string basecall_model_path(const std::string& bc_gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/"
            + basecall_strand_subgroup(st) + "/Model";
    }
    static std::string basecall_model_file_path(const std::string& bc_gr, unsigned st)
    {
        assert(st < 2);
        return basecall_root_path() + "/" + basecall_group_prefix() + bc_gr
            + "/Summary/basecall_1d_" + (st == 0? "template" : "complement") + "/model_file";
    }
    static std::string basecall_events_path(const std::string& bc_gr, unsigned st)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/"
            + basecall_strand_subgroup(st) + "/Events";
    }
    static std::string basecall_event_alignment_path(const std::string& bc_gr)
    {
        return basecall_root_path() + "/" + basecall_group_prefix() + bc_gr + "/"
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
            else
            {
                std::clog << "skipping dataset: " << path + "/" + sg << std::endl;
            }
        }
    } // rec_copy_attributes()

    static fast5_pack::Huffman_Diff_Coder & rw_coder()
    {
        static fast5_pack::Huffman_Diff_Coder _rw_coder;
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
    
    static fast5_pack::Huffman_Diff_Coder & ed_skip_coder()
    {
        static fast5_pack::Huffman_Diff_Coder _ed_skip_coder;
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

    static fast5_pack::Huffman_Diff_Coder & ed_len_coder()
    {
        static fast5_pack::Huffman_Diff_Coder _ed_len_coder;
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

}; // class File

} // namespace fast5

#endif
