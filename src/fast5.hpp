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

#include "hdf5_tools.hpp"
#define MAX_K_LEN 8

namespace fast5
{

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

struct Channel_Id_Parameters
{
    std::string channel_number;
    double digitisation;
    double offset;
    double range;
    double sampling_rate;
}; // struct Channel_Id_Parameters

typedef uint16_t Raw_Sample_Type;

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

    File() = default;
    File(const std::string& file_name) { open(file_name); }

    void open(const std::string& file_name)
    {
        Base::open(file_name);
        // detect basecall version
        if (is_open())
        {
            if (have_basecall_version())
            {
                _basecall_version = basecall_version();
            }
        }
    }

    std::string file_version() const
    {
        double v;
        assert(Base::exists("/file_version"));
        Base::read< double >("/file_version", v);
        // convert it to string
        std::ostringstream os;
        os << v;
        return os.str();
    }

    bool have_basecall_version() const
    {
        return Base::exists(get_bc_2d_root() + "/chimaera version");
    }

    std::string basecall_version() const
    {
        std::string res;
        std::string path = get_bc_2d_root() + "/chimaera version";
        assert(Base::exists(path));
        Base::read< std::string >(path, res);
        return res;
    }

    std::string eventdetection_version() const
    {
        std::string res;
        // only support eventdetection group 000 for now
        std::string path = "/Analyses/EventDetection_000/version";
        assert(Base::exists(path));
        Base::read< std::string >(path, res);
        return res;
    }

    std::string get_log() const
    {
        std::string res;
        std::string path = get_bc_2d_root() + "/Log";
        assert(Base::exists(path));
        Base::read< std::string >(path, res);
        return res;
    }

    bool have_channel_id_params() const
    {
        return Base::group_exists("/UniqueGlobalKey/channel_id");
    }

    Channel_Id_Parameters get_channel_id_params() const
    {
        Channel_Id_Parameters res;
        std::string path("/UniqueGlobalKey/channel_id");
        Base::read< decltype(res.channel_number) >(path + "/channel_number", res.channel_number);
        Base::read< decltype(res.digitisation) >(path + "/digitisation", res.digitisation);
        Base::read< decltype(res.offset) >(path + "/offset", res.offset);
        Base::read< decltype(res.range) >(path + "/range", res.range);
        Base::read< decltype(res.sampling_rate) >(path + "/sampling_rate", res.sampling_rate);
        return res;
    }

    double get_sampling_rate() const
    {
        assert(have_sampling_rate());
        double res;
        Base::read< double >("/UniqueGlobalKey/channel_id/sampling_rate", res);
        return res;
    }

    bool have_sampling_rate() const
    {
        return Base::exists("/UniqueGlobalKey/channel_id/sampling_rate");
    }

    std::string get_model_file(size_t i) const
    {
        std::string res;
        assert(Base::exists(model_file_path(i)));
        Base::read< std::string >(model_file_path(i), res);
        return res;
    }

    std::string sequences_version() const
    {
        std::vector< std::string > tmp;
        assert(Base::exists("/Sequences/Meta/version"));
        Base::read< std::string >("/Sequences/Meta/version", tmp);
        std::string res;
        for (const auto& s: tmp)
        {
            res += s;
        }
        return res;
    }

    bool have_basecalled_group() const
    {
        return Base::group_exists(get_bc_2d_root());
    }

    bool have_eventdetection_group() const
    {
        return Base::group_exists(get_eventdetection_root());
    }

    bool have_sequences_group() const
    {
        return Base::group_exists("/Sequences");
    }

    bool have_basecalled_2D() const
    {
        return Base::exists(get_bc_2d_root() + "/BaseCalled_2D/Fastq");
    }

    std::string basecalled_1D(size_t i) const
    {
        std::string res;
        Base::read< std::string >(fastq_path(i), res);
        return fastq2sequence(res);
    }

    std::string basecalled_2D() const
    {
        std::string res;
        Base::read< std::string >(get_bc_2d_root() + "/BaseCalled_2D/Fastq", res);
        return fastq2sequence(res);
    }

    std::vector< Event_Alignment_Entry > get_event_alignments() const
    {
        std::vector< Event_Alignment_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("template", &Event_Alignment_Entry::template_index);
        m.add_member("complement", &Event_Alignment_Entry::complement_index);
        m.add_member("kmer", &Event_Alignment_Entry::kmer);
        Base::read< Event_Alignment_Entry >(get_bc_2d_root() + "/BaseCalled_2D/Alignment", res, &m);
        return res;
    }

    bool have_model(size_t i) const
    {
        return Base::exists(model_path(i));
    }
    bool have_events(size_t i) const
    {
        return Base::exists(events_path(i));
    }

    std::vector< Model_Entry > get_model(size_t i) const
    {
        std::vector< Model_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("kmer", &Model_Entry::kmer);
        m.add_member("level_mean", &Model_Entry::level_mean);
        m.add_member("level_stdv", &Model_Entry::level_stdv);
        m.add_member("sd_mean", &Model_Entry::sd_mean);
        m.add_member("sd_stdv", &Model_Entry::sd_stdv);
        Base::read< Model_Entry >(model_path(i), res, &m);
        return res;
    }

    Model_Parameters get_model_parameters(size_t i) const
    {
        Model_Parameters res;
        std::string path = model_path(i);
        Base::read< double >(path + "/drift", res.drift);
        Base::read< double >(path + "/scale", res.scale);
        Base::read< double >(path + "/scale_sd", res.scale_sd);
        Base::read< double >(path + "/shift", res.shift);
        Base::read< double >(path + "/var", res.var);
        Base::read< double >(path + "/var_sd", res.var_sd);
        return res;
    }

    std::vector< Event_Entry > get_events(size_t i) const
    {
        std::vector< Event_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("mean", &Event_Entry::mean);
        m.add_member("start", &Event_Entry::start);
        m.add_member("stdv", &Event_Entry::stdv);
        m.add_member("length", &Event_Entry::length);
        m.add_member("model_state", &Event_Entry::model_state);
        m.add_member("move", &Event_Entry::move);
        Base::read< Event_Entry >(events_path(i), res, &m);
        return res;
    }

    bool have_eventdetection_events() const
    {
        if (not Base::group_exists(get_eventdetection_root() + "/Reads")) return false;
        auto reads = Base::list_group(get_eventdetection_root() + "/Reads");
        return not reads.empty();
    }

    std::string get_eventdetection_read_name() const
    {
        auto reads = Base::list_group(get_eventdetection_root() + "/Reads");
        assert(not reads.empty());
        return reads[0];
    }

    std::vector< EventDetection_Event_Entry > get_eventdetection_events() const
    {
        std::vector< EventDetection_Event_Entry > res;
        auto path = get_eventdetection_root() + "/Reads/" + get_eventdetection_read_name() + "/Events";
        auto struct_member_names = Base::get_struct_members(path);
        assert(struct_member_names.size() >= 4);
        bool have_stdv = false;
        bool have_variance = false;
        for (const auto& s : struct_member_names)
        {
            if (s == "stdv") have_stdv = true;
            else if (s == "variance") have_variance = true;
        }
        assert(have_stdv or have_variance);
        hdf5_tools::Compound_Map m;
        m.add_member("mean", &EventDetection_Event_Entry::mean);
        m.add_member("start", &EventDetection_Event_Entry::start);
        m.add_member("length", &EventDetection_Event_Entry::length);
        if (have_stdv)
        {
            m.add_member("stdv", &EventDetection_Event_Entry::stdv);
        }
        else
        {
            m.add_member("variance", &EventDetection_Event_Entry::stdv);
        }
        Base::read< EventDetection_Event_Entry >(path, res, &m);
        if (not have_stdv)
        {
            // have read variances
            for (auto& e : res)
            {
                e.stdv = std::sqrt(e.stdv);
            }
        }
        return res;
    }

    EventDetection_Event_Parameters get_eventdetection_event_parameters() const
    {
        EventDetection_Event_Parameters res;
        auto path = get_eventdetection_root() + "/Reads/" + get_eventdetection_read_name();

        Base::read< decltype(res.abasic_found) >(path + "/abasic_found", res.abasic_found);
        /*
        if (res.abasic_found)
        {
            Base::read< decltype(res.abasic_event_index) >(path + "/abasic_event_index", res.abasic_event_index);
            Base::read< decltype(res.abasic_peak_height) >(path + "/abasic_peak_height", res.abasic_peak_height);
        }
        Base::read< decltype(res.hairpin_found) >(path + "/hairpin_found", res.hairpin_found);
        if (res.hairpin_found)
        {
            Base::read< decltype(res.hairpin_event_index) >(path + "/hairpin_event_index", res.hairpin_event_index);
            Base::read< decltype(res.hairpin_peak_height) >(path + "/hairpin_peak_height", res.hairpin_peak_height);
            Base::read< decltype(res.hairpin_polyt_level) >(path + "/hairpin_polyt_level", res.hairpin_polyt_level);
        }
        */
        Base::read< decltype(res.duration) >(path + "/duration", res.duration);
        Base::read< decltype(res.median_before) >(path + "/median_before", res.median_before);
        if (Base::exists(path + "/read_id"))
        {
            Base::read< decltype(res.read_id) >(path + "/read_id", res.read_id);
            // hack: read_id is currently wrong
            res.read_id.resize(res.read_id.size() - 1);
        }
        Base::read< decltype(res.read_number) >(path + "/read_number", res.read_number);
        Base::read< decltype(res.scaling_used) >(path + "/scaling_used", res.scaling_used);
        Base::read< decltype(res.start_mux) >(path + "/start_mux", res.start_mux);
        Base::read< decltype(res.start_time) >(path + "/start_time", res.start_time);
        return res;
    }

    //
    // Access to raw samples.
    //
    bool have_raw_samples_read_name() const
    {
        if (not Base::group_exists("/Raw/Reads")) return false;
        auto reads = Base::list_group("/Raw/Reads");
        return not reads.empty();
    }

    std::string get_raw_samples_read_name() const
    {
        assert(have_raw_samples_read_name());
        auto reads = Base::list_group("/Raw/Reads");
        return reads[0];
    }

    bool have_raw_samples() const
    {
        return (have_channel_id_params()
                and have_raw_samples_read_name()
                and Base::dataset_exists(std::string("/Raw/Reads/") + get_raw_samples_read_name() + "/Signal"));
    }

    std::vector< float > get_raw_samples() const
    {
        assert(have_raw_samples());
        // get raw samples
        std::vector< Raw_Sample_Type > raw_samples;
        auto path = std::string("/Raw/Reads/") + get_raw_samples_read_name() + "/Signal";
        Base::read< Raw_Sample_Type >(path, raw_samples);
        // get scaling parameters
        auto params = get_channel_id_params();
        // decode levels
        std::vector< float > res;
        res.reserve(raw_samples.size());
        for (auto int_level : raw_samples)
        {
            res.push_back((static_cast< float >(int_level) + params.offset) * params.range / params.digitisation);
        }
        return res;
    }

    void set_eventdetection_group_id(size_t i)
    {
        assert(i <= 999);
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(3) << i;
        _eventdetection_group_id = ss.str();
    }

    void set_basecalled_group_id(size_t i)
    {
        assert(i <= 999);
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(3) << i;
        _basecalled_group_id = ss.str();
    }

    static bool is_valid_file(const std::string& file_name)
    {
        if (not Base::is_valid_file(file_name))
        {
            return false;
        }
        File f(file_name);
        return f.have_eventdetection_group();
    }

private:

    std::string fastq2sequence(std::string& str) const
    {
        // Split the FASTQ record on newlines
        size_t nl1 = str.find_first_of('\n');
        size_t nl2 = str.find_first_of('\n', nl1 + 1);

        if(nl1 == std::string::npos || nl2 == std::string::npos)
            return "";
        else
            return str.substr(nl1 + 1, nl2 - nl1 - 1);
    }

    std::string get_eventdetection_root() const
    {
        return "/Analyses/EventDetection_" + _eventdetection_group_id;
    }

    // Returns the root path of the form:
    // /Analyses/Basecall_2D_ddd/ where ddd is the group
    std::string get_bc_2d_root() const
    {
        return "/Analyses/Basecall_2D_" + _basecalled_group_id;
    }
    // /Analyses/Basecall_1D_ddd/ where ddd is the group
    std::string get_bc_1d_root() const
    {
        return "/Analyses/Basecall_1D_" + _basecalled_group_id;
    }

    std::string model_path(size_t i) const
    {
        static std::vector< std::string > _model_path =
            { "/BaseCalled_template/Model",
              "/BaseCalled_complement/Model" };
        if (_basecall_version < std::string("1.16"))
        {
            return get_bc_2d_root() + _model_path.at(i);
        }
        else
        {
            return get_bc_1d_root() + _model_path.at(i);
        }
    }

    std::string events_path(size_t i) const
    {
        static std::vector< std::string > _events_path =
            { "/BaseCalled_template/Events",
              "/BaseCalled_complement/Events" };
        if (_basecall_version < std::string("1.16"))
        {
            return get_bc_2d_root() + _events_path.at(i);
        }
        else
        {
            return get_bc_1d_root() + _events_path.at(i);
        }
    }

    std::string model_file_path(size_t i) const
    {
        static std::vector< std::string > _model_file_path =
            { "/Summary/basecall_1d_template/model_file",
              "/Summary/basecall_1d_complement/model_file" };
        if (_basecall_version < std::string("1.16"))
        {
            return get_bc_2d_root() + _model_file_path.at(i);
        }
        else
        {
            return get_bc_1d_root() + _model_file_path.at(i);
        }
    }

    std::string fastq_path(size_t i) const
    {
        static std::vector< std::string > _fastq_path =
            { "/BaseCalled_template/Fastq",
              "/BaseCalled_complement/Fastq" };

        if (_basecall_version < std::string("1.16"))
        {
            return get_bc_2d_root() + _fastq_path.at(i);
        }
        else
        {
            return get_bc_1d_root() + _fastq_path.at(i);
        }
    }

    // default to 000 event detection group
    std::string _eventdetection_group_id = "000";

    // default to using the 000 analysis group
    std::string _basecalled_group_id = "000";

    // basecall (chimera) version
    std::string _basecall_version;

}; // class File

} // namespace fast5

#endif
