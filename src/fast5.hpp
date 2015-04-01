#ifndef __FAST5_HPP
#define __FAST5_HPP

#include <cassert>
#include <exception>
#include <iostream>
#include <sstream>
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
        assert(Base::exists("/Analyses/Basecall_2D_000/version"));
        Base::read< std::string >("/Analyses/Basecall_2D_000/version", res);
        return res;
    }

    std::string eventdetection_version() const
    {
        std::string res;
        assert(Base::exists("/Analyses/EventDetection_000/version"));
        Base::read< std::string >("/Analyses/EventDetection_000/version", res);
        return res;
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
        return Base::exists("/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq");
    }

    std::string basecalled_2D() const
    {
        std::string res;
        Base::read< std::string >("/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq", res);
        
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
        Base::read< Event_Alignment_Entry >("/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment", res, &m);
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

private:
    static const std::string& model_path(size_t i)
    {
        static std::vector< std::string > _model_path =
            { "/Analyses/Basecall_2D_000/BaseCalled_template/Model",
              "/Analyses/Basecall_2D_000/BaseCalled_complement/Model" };
        return _model_path.at(i);
    }

    static const std::string& events_path(size_t i)
    {
        static std::vector< std::string > _events_path =
            { "/Analyses/Basecall_2D_000/BaseCalled_template/Events",
              "/Analyses/Basecall_2D_000/BaseCalled_complement/Events" };
        return _events_path.at(i);
    }

    static const std::string& model_file_path(size_t i)
    {
        static std::vector< std::string > _model_file_path =
            { "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/model_file",
              "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/model_file" };
        return _model_file_path.at(i);
    }

}; // class File

} // namespace fast5

#endif
