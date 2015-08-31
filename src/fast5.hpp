#ifndef __FAST5_HPP
#define __FAST5_HPP

#include <cassert>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include "hdf5_tools.hpp"

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
    char kmer[6];
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
    char model_state[6];
    double model_level;
    long long move;
    double p_model_state;
    char mp_state[6];
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
    char kmer[6];
}; // struct Event_Alignment_Entry

class File
    : private hdf5_tools::File_Reader
{
private:
    typedef hdf5_tools::File_Reader Base;
public:
    using Base::Base;

    using Base::is_open;
    using Base::file_name;
    using Base::open;
    using Base::close;

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

    std::string basecall_version() const
    {
        std::string res;
        std::string path = get_bc_2d_root() + "/version";
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

    double get_sampling_rate() const
    {
        assert(have_sampling_rate());

        auto lg = get_log();
        auto idx = lg.find("Sampling rate is");

        std::string line;
        std::stringstream ss1(lg.substr(idx));
        std::getline(ss1,line,'\n');

        std::stringstream ss2(line);

        std::string token;
        std::getline(ss2,token,' ');    //Sampling
        std::getline(ss2,token,' ');    //rate
        std::getline(ss2,token,' ');    //is
        std::getline(ss2,token,' ');    //Hz value

        return std::atof(token.c_str());
    }

    bool have_sampling_rate() const
    {
        auto lg = get_log();
        auto idx = lg.find("Sampling rate is");
        return idx != std::string::npos;
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

    bool have_basecalled_2D() const
    {
        return Base::exists(get_bc_2d_root() + "/BaseCalled_2D/Fastq");
    }

    std::string basecalled_2D() const
    {
        std::string res;
        Base::read< std::string >(get_bc_2d_root() + "/BaseCalled_2D/Fastq", res);
        
        // Split the FASTQ record on newlines
        size_t nl1 = res.find_first_of('\n');
        size_t nl2 = res.find_first_of('\n', nl1 + 1);

        if(nl1 == std::string::npos || nl2 == std::string::npos)
            return "";
        else
            return res.substr(nl1 + 1, nl2 - nl1 - 1);
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
        Base::read< Event_Entry >(events_path(i), res, &m);
        return res;
    }

    void set_basecalled_group_id(size_t i)
    {
        assert(i <= 999);
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(3) << i;
        _basecalled_group_id = ss.str();
    }


private:
    
    // Returns the root path of the form:
    // Analyses/Basecall_2D_ddd/ where ddd is the group
    std::string get_bc_2d_root() const
    {
        return "/Analyses/Basecall_2D_" + _basecalled_group_id;
    }

    std::string model_path(size_t i) const
    {
        static std::vector< std::string > _model_path =
            { "/BaseCalled_template/Model",
              "/BaseCalled_complement/Model" };
        return get_bc_2d_root() + _model_path.at(i);
    }

    std::string events_path(size_t i) const
    {
        static std::vector< std::string > _events_path =
            { "/BaseCalled_template/Events",
              "/BaseCalled_complement/Events" };
        return get_bc_2d_root() + _events_path.at(i);
    }

    std::string model_file_path(size_t i) const
    {
        static std::vector< std::string > _model_file_path =
            { "/Summary/basecall_1d_template/model_file",
              "/Summary/basecall_1d_complement/model_file" };
        return get_bc_2d_root() + _model_file_path.at(i);
    }

    // default to using the 000 analysis group
    std::string _basecalled_group_id = "000";

}; // class File

} // namespace fast5

#endif
