#ifndef __FAST5_HPP
#define __FAST5_HPP

#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>

#include "hdf5_tools.hpp"
#define MAX_K_LEN 8

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

struct Raw_Samples_Parameters
{
    std::string read_id;
    long read_number;
    long start_mux;
    long start_time;
    long duration;
}; // struct Raw_Samples_Parameters

struct EventDetection_Event_Entry
{
    double mean;
    double stdv;
    long long start;
    long long length;
}; // struct EventDetection_Event

struct EventDetection_Event_Parameters
{
    unsigned abasic_found;
    /*
    unsigned abasic_event_index;
    double abasic_peak_height;
    unsigned hairpin_found;
    unsigned hairpin_event_index;
    double hairpin_peak_height;
    double hairpin_polyt_level;
    */
    long duration;
    double median_before;
    std::string read_id;
    long read_number;
    long scaling_used;
    long start_mux;
    long start_time;
}; // struct EventDetection_Event_Parameters

//
// This struct represents the expected signal measured
// given the kmer sequence that is in the pore when the
// the observations are made. A pore model consists
// of 1024 of these entries (one per 5-mer) and global
// shift/scaling parameters.
//
struct Model_Entry
{
    char kmer[MAX_K_LEN];
    long long variant;
    double level_mean;
    double level_stdv;
    double sd_mean;
    double sd_stdv;
    double weight;
}; // struct Model_Entry

//
// This struct represents the global transformations
// that must be applied to each Model_Entry
//
struct Model_Parameters
{
    double drift;
    double scale;
    double scale_sd;
    double shift;
    double var;
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
    double start;
    double stdv;
    double length;
    char model_state[MAX_K_LEN];
    double model_level;
    long long move;
    double p_model_state;
    char mp_state[MAX_K_LEN];
    double p_mp_state;
    double p_A;
    double p_C;
    double p_G;
    double p_T;
}; // struct Event_Entry

//
// This struct represents a template-to-complement
// match that is emitted by ONT's 2D basecaller
//
struct Event_Alignment_Entry
{
    long long template_index;
    long long complement_index;
    char kmer[MAX_K_LEN];
}; // struct Event_Alignment_Entry


class File
    : private hdf5_tools::File_Reader
{
private:
    typedef hdf5_tools::File_Reader Base;
public:
    using Base::is_open;
    using Base::file_name;
    using Base::close;
    using Base::get_object_count;
    using Base::is_valid_file;

    File() = default;
    File(const std::string& file_name) { open(file_name); }

    void open(const std::string& file_name)
    {
        Base::open(file_name);
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
        Base::read< std::string >(file_version_path(), res);
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
        Base::read< decltype(res.channel_number) >(channel_id_path() + "/channel_number", res.channel_number);
        Base::read< decltype(res.digitisation) >(channel_id_path() + "/digitisation", res.digitisation);
        Base::read< decltype(res.offset) >(channel_id_path() + "/offset", res.offset);
        Base::read< decltype(res.range) >(channel_id_path() + "/range", res.range);
        Base::read< decltype(res.sampling_rate) >(channel_id_path() + "/sampling_rate", res.sampling_rate);
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
     */
    bool have_raw_samples() const
    {
        return have_channel_id_params() and not get_raw_samples_read_name_list().empty();
    }
    /**
     * Get raw samples attributes for first read name.
     */
    Raw_Samples_Parameters get_raw_samples_params() const
    {
        return get_raw_samples_params(get_raw_samples_read_name_list().front());
    }
    /**
     * Get raw samples attributes for given read name.
     */
    Raw_Samples_Parameters get_raw_samples_params(const std::string& rn) const
    {
        Raw_Samples_Parameters res;
        std::string p = raw_samples_params_path(rn);
        Base::read< decltype(res.read_id) >(p + "/read_id", res.read_id);
        Base::read< decltype(res.read_number) >(p + "/read_number", res.read_number);
        Base::read< decltype(res.start_mux) >(p + "/start_mux", res.start_mux);
        Base::read< decltype(res.start_time) >(p + "/start_time", res.start_time);
        Base::read< decltype(res.duration) >(p + "/duration", res.duration);
        return res;
    }
    /**
     * Get raw samples for first read name.
     */
    std::vector< Raw_Samples_Entry > get_raw_samples() const
    {
        return get_raw_samples(get_raw_samples_read_name_list().front());
    }
    /**
     * Get raw samples for given read name.
     */
    std::vector< Raw_Samples_Entry > get_raw_samples(const std::string& rn) const
    {
        // get raw samples
        std::vector< uint16_t > raw_samples;
        Base::read< uint16_t >(raw_samples_path(rn), raw_samples);
        // get scaling parameters
        auto channel_id_params = get_channel_id_params();
        // decode levels
        std::vector< Raw_Samples_Entry > res;
        res.reserve(raw_samples.size());
        for (auto int_level : raw_samples)
        {
            res.push_back((static_cast< float >(int_level) + channel_id_params.offset)
                          * channel_id_params.range / channel_id_params.digitisation);
        }
        return res;
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
     * Get list of reads for first EventDetection group.
     */
    std::vector< std::string > get_eventdetection_read_name_list() const
    {
        return get_eventdetection_read_name_list(get_eventdetection_group_list().front());
    }
    /**
     * Get list of reads for given EventDetection group.
     */
    std::vector< std::string > get_eventdetection_read_name_list(const std::string& ed_gr) const
    {
        return detect_eventdetection_read_name_list(ed_gr);
    }
    /**
     * Check if EventDetection events exist for first EventDetection group.
     */
    bool have_eventdetection_events() const
    {
        return have_eventdetection_events(get_eventdetection_group_list().front());
    }
    /**
     * Check if EventDetection events exist for given EventDetection group.
     */
    bool have_eventdetection_events(const std::string& ed_gr) const
    {
        return not get_eventdetection_read_name_list(ed_gr).empty();
    }
    /**
     * Get EventDetection params for first EventDetection group.
     */
    std::map< std::string, std::string > get_eventdetection_params() const
    {
        return get_eventdetection_params(get_eventdetection_group_list().front());
    }
    /**
     * Get EventDetection params for given EventDetection group.
     */
    std::map< std::string, std::string > get_eventdetection_params(const std::string& ed_gr) const
    {
        return get_attr_map(eventdetection_params_path(ed_gr));
    }
    /**
     * Get EventDetection event params for first EventDetection group, and first read name.
     */
    EventDetection_Event_Parameters get_eventdetection_event_params() const
    {
        return get_eventdetection_event_params(get_eventdetection_group_list().front());
    }
    /**
     * Get EventDetection event params for given EventDetection group, and first read name.
     */
    EventDetection_Event_Parameters get_eventdetection_event_params(const std::string& ed_gr) const
    {
        return get_eventdetection_event_params(ed_gr, get_eventdetection_read_name_list(ed_gr).front());
    }
    /**
     * Get EventDetection event params for given EventDetection group, and given read name.
     */
    EventDetection_Event_Parameters get_eventdetection_event_params(const std::string& ed_gr, const std::string& rn) const
    {
        EventDetection_Event_Parameters res;
        auto p = eventdetection_event_params_path(ed_gr, rn);
        Base::read< decltype(res.abasic_found) >(p + "/abasic_found", res.abasic_found);
        /*
        if (res.abasic_found)
        {
            Base::read< decltype(res.abasic_event_index) >(p + "/abasic_event_index", res.abasic_event_index);
            Base::read< decltype(res.abasic_peak_height) >(p + "/abasic_peak_height", res.abasic_peak_height);
        }
        Base::read< decltype(res.hairpin_found) >(p + "/hairpin_found", res.hairpin_found);
        if (res.hairpin_found)
        {
            Base::read< decltype(res.hairpin_event_index) >(p + "/hairpin_event_index", res.hairpin_event_index);
            Base::read< decltype(res.hairpin_peak_height) >(p + "/hairpin_peak_height", res.hairpin_peak_height);
            Base::read< decltype(res.hairpin_polyt_level) >(p + "/hairpin_polyt_level", res.hairpin_polyt_level);
        }
        */
        Base::read< decltype(res.duration) >(p + "/duration", res.duration);
        Base::read< decltype(res.median_before) >(p + "/median_before", res.median_before);
        if (Base::exists(p + "/read_id"))
        {
            Base::read< decltype(res.read_id) >(p + "/read_id", res.read_id);
        }
        Base::read< decltype(res.read_number) >(p + "/read_number", res.read_number);
        Base::read< decltype(res.scaling_used) >(p + "/scaling_used", res.scaling_used);
        Base::read< decltype(res.start_mux) >(p + "/start_mux", res.start_mux);
        Base::read< decltype(res.start_time) >(p + "/start_time", res.start_time);
        return res;
    }
    /**
     * Get EventDetection events for first EventDetection group, and first read name.
     */
    std::vector< EventDetection_Event_Entry > get_eventdetection_events() const
    {
        return get_eventdetection_events(get_eventdetection_group_list().front());
    }
    /**
     * Get EventDetection events for given EventDetection group, and first read name.
     */
    std::vector< EventDetection_Event_Entry > get_eventdetection_events(const std::string& ed_gr) const
    {
        return get_eventdetection_events(ed_gr, get_eventdetection_read_name_list(ed_gr).front());
    }
    /**
     * Get EventDetection events for given EventDetection group, and given read name.
     */
    std::vector< EventDetection_Event_Entry > get_eventdetection_events(const std::string& ed_gr, const std::string& rn) const
    {
        std::vector< EventDetection_Event_Entry > res;
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
        Base::read< EventDetection_Event_Entry >(p, res, &m);
        if (not have_stdv)
        {
            // have read variances
            for (auto& e : res)
            {
                e.stdv = std::sqrt(e.stdv);
            }
        }
        return res;
    } // get_eventdetection_events()

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
    const std::vector< std::string >& get_basecall_group_list(unsigned st) const
    {
        return _basecall_strand_group_list[st];
    }
    /**
     * Check if any Basecall groups exist for given strand.
     */
    bool have_basecall_groups(unsigned st) const
    {
        return not get_basecall_group_list(st).empty();
    }
    /**
     * Get Basecall group params for first Basecall group for given strand.
     */
    std::map< std::string, std::string > get_basecall_params(unsigned st) const
    {
        return get_basecall_params(get_basecall_group_list(st).front());
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
        Base::read< std::string >(path, res);
        return res;
    }
    /**
     * Check if Basecall fastq exists for first Basecall group for given strand.
     */
    bool have_basecall_fastq(unsigned st) const
    {
        return not get_basecall_group_list(st).empty()
            and have_basecall_fastq(get_basecall_group_list(st).front(), st);
    }
    /**
     * Check if Basecall fastq exists for given Basecall group and given strand.
     */
    bool have_basecall_fastq(const std::string& bc_gr, unsigned st) const
    {
        return Base::dataset_exists(basecall_fastq_path(bc_gr, st));
    }
    /**
     * Get Basecall fastq for first Basecall group for given strand.
     */
    std::string get_basecall_fastq(unsigned st) const
    {
        return get_basecall_fastq(get_basecall_group_list(st).front(), st);
    }
    /**
     * Get Basecall fastq for given Basecall group and given strand.
     */
    std::string get_basecall_fastq(const std::string& bc_gr, unsigned st) const
    {
        std::string res;
        Base::read< std::string >(basecall_fastq_path(bc_gr, st), res);
        return res;
    }
    /**
     * Check if Basecall seq exists for first Basecall group for given strand.
     */
    bool have_basecall_seq(unsigned st) const
    {
        return have_basecall_fastq(st);
    }
    /**
     * Check if Basecall seq exists for given Basecall group and given strand.
     */
    bool have_basecall_seq(const std::string& bc_gr, unsigned st) const
    {
        return have_basecall_fastq(bc_gr, st);
    }
    /**
     * Get Basecall sequence for first Basecall group for given strand.
     */
    std::string get_basecall_seq(unsigned st) const
    {
        return get_basecall_seq(get_basecall_group_list(st).front(), st);
    }
    /**
     * Get Basecall sequence for given Basecall group and given strand.
     */
    std::string get_basecall_seq(const std::string& bc_gr, unsigned st) const
    {
        return fq2seq(get_basecall_fastq(bc_gr, st));
    }
    /**
     * Check if Basecall model exist for first Basecall group for given strand.
     */
    bool have_basecall_model(unsigned st) const
    {
        return not get_basecall_group_list(st).empty()
            and have_basecall_model(get_basecall_group_list(st).front(), st);
    }
    /**
     * Check if Basecall model exist for given Basecall group and given strand.
     */
    bool have_basecall_model(const std::string& bc_gr, unsigned st) const
    {
        return Base::dataset_exists(basecall_model_path(bc_gr, st));
    }
    /**
     * Get Basecall model file name for first Basecall group for given strand.
     */
    std::string get_basecall_model_file(unsigned st) const
    {
        return get_basecall_model_file(get_basecall_group_list(st).front(), st);
    }
    /**
     * Get Basecall model file name for given Basecall group and given strand.
     */
    std::string get_basecall_model_file(const std::string& bc_gr, unsigned st) const
    {
        std::string res;
        assert(Base::exists(basecall_model_file_path(bc_gr, st)));
        Base::read< std::string >(basecall_model_file_path(bc_gr, st), res);
        return res;
    }
    /**
     * Get Basecall model parameters for first Basecall group for given strand.
     */
    Model_Parameters get_basecall_model_params(unsigned st) const
    {
        return get_basecall_model_params(get_basecall_group_list(st).front(), st);
    }
    /**
     * Get Basecall model parameters for given Basecall group and given strand.
     */
    Model_Parameters get_basecall_model_params(const std::string& bc_gr, unsigned st) const
    {
        Model_Parameters res;
        std::string p = basecall_model_path(bc_gr, st);
        Base::read< double >(p + "/drift", res.drift);
        Base::read< double >(p + "/scale", res.scale);
        Base::read< double >(p + "/scale_sd", res.scale_sd);
        Base::read< double >(p + "/shift", res.shift);
        Base::read< double >(p + "/var", res.var);
        Base::read< double >(p + "/var_sd", res.var_sd);
        return res;
    }
    /**
     * Get Basecall model for first Basecall group for given strand.
     */
    std::vector< Model_Entry > get_basecall_model(unsigned st) const
    {
        return get_basecall_model(get_basecall_group_list(st).front(), st);
    }
    /**
     * Get Basecall model for given Basecall group and given strand.
     */
    std::vector< Model_Entry > get_basecall_model(const std::string& bc_gr, unsigned st) const
    {
        std::vector< Model_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("kmer", &Model_Entry::kmer);
        m.add_member("level_mean", &Model_Entry::level_mean);
        m.add_member("level_stdv", &Model_Entry::level_stdv);
        m.add_member("sd_mean", &Model_Entry::sd_mean);
        m.add_member("sd_stdv", &Model_Entry::sd_stdv);
        Base::read< Model_Entry >(basecall_model_path(bc_gr, st), res, &m);
        return res;
    }
    /**
     * Check if Basecall events exist for first Basecall group for given strand.
     */
    bool have_basecall_events(unsigned st) const
    {
        return not get_basecall_group_list(st).empty()
            and have_basecall_events(get_basecall_group_list(st).front(), st);
    }
    /**
     * Check if Basecall events exist for given Basecall group and given strand.
     */
    bool have_basecall_events(const std::string& bc_gr, unsigned st) const
    {
        return Base::dataset_exists(basecall_events_path(bc_gr, st));
    }
    /**
     * Get Basecall events for first Basecall group for given strand.
     */
    std::vector< Event_Entry > get_basecall_events(unsigned st) const
    {
        return get_basecall_events(get_basecall_group_list(st).front(), st);
    }
    /**
     * Get Basecall events for given Basecall group and given strand.
     */
    std::vector< Event_Entry > get_basecall_events(const std::string& bc_gr, unsigned st) const
    {
        std::vector< Event_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("mean", &Event_Entry::mean);
        m.add_member("start", &Event_Entry::start);
        m.add_member("stdv", &Event_Entry::stdv);
        m.add_member("length", &Event_Entry::length);
        m.add_member("model_state", &Event_Entry::model_state);
        m.add_member("move", &Event_Entry::move);
        Base::read< Event_Entry >(basecall_events_path(bc_gr, st), res, &m);
        return res;
    }
    /**
     * Check if Basecall event alignment exist for first Basecall group with 2d calls.
     */
    bool have_basecall_event_alignment() const
    {
        return not get_basecall_group_list(2).empty()
            and have_basecall_event_alignment(get_basecall_group_list(2).front());
    }
    /**
     * Check if Basecall event alignment exist for given Basecall group.
     */
    bool have_basecall_event_alignment(const std::string& bc_gr) const
    {
        return Base::dataset_exists(basecall_event_alignment_path(bc_gr));
    }
    /**
     * Get Basecall event alignment for first Basecall group with 2d calls.
     */
    std::vector< Event_Alignment_Entry > get_basecall_event_alignment() const
    {
        return get_basecall_event_alignment(get_basecall_group_list(2).front());
    }
    /**
     * Get Basecall events for given Basecall group.
     */
    std::vector< Event_Alignment_Entry > get_basecall_event_alignment(const std::string& bc_gr) const
    {
        std::vector< Event_Alignment_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("template", &Event_Alignment_Entry::template_index);
        m.add_member("complement", &Event_Alignment_Entry::complement_index);
        m.add_member("kmer", &Event_Alignment_Entry::kmer);
        Base::read< Event_Alignment_Entry >(basecall_event_alignment_path(bc_gr), res, &m);
        return res;
    }

    static std::string fq2seq(const std::string& fq)
    {
        size_t nl1_pos = fq.find_first_of('\n');
        if (nl1_pos == std::string::npos) return std::string();
        size_t nl2_pos = fq.find_first_of('\n', nl1_pos + 1);
        if (nl2_pos == std::string::npos) return std::string();
        return fq.substr(nl1_pos + 1, nl2_pos - nl1_pos - 1);
    }

private:
    void detect_raw_samples_read_name_list()
    {
        if (not Base::group_exists(raw_samples_root_path())) return;
        auto rn_list = Base::list_group(raw_samples_root_path());
        for (const auto& rn : rn_list)
        {
            if (not Base::dataset_exists(raw_samples_path(rn))) continue;
            _raw_samples_read_name_list.push_back(rn);
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
            if (not Base::dataset_exists(p + "/" + rn + "/Events")) continue;
            res.push_back(rn);
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
            Base::read< std::string >(path + "/" + a, tmp);
            res[a] = tmp;
        }
        return res;
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
}; // class File

} // namespace fast5

#endif
