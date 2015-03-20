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
using namespace hdf5_tools;

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

struct Event_Alignment_Entry
{
    long long template_index;
    long long complement_index;
    char kmer[6];
}; // struct Event_Alignment_Entry

struct Model_Parameters
{
    double drift;
    double scale;
    double scale_sd;
    double shift;
    double var;
    double var_sd;
}; // struct Model_Parameters

class File
{
public:
    File() : _file_id(0) {}
    File(const std::string& file_name) : _file_id(0) { open(file_name); }
    File(const File&) = delete;
    File(File&&) = default;
    File& operator = (const File&) = delete;
    File& operator = (File&&) = default;
    ~File()
    {
        if (is_open())
        {
            close();
        }
    }

    bool is_open() const { return _file_id > 0; }

    void open(const std::string& file_name)
    {
        assert(not is_open());
        _file_name = file_name;
        _file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (not is_open()) throw Exception(_file_name + ": error in H5Fopen");
    }

    void close()
    {
        assert(is_open());
        assert(H5Fget_obj_count(_file_id, H5F_OBJ_ALL) == 1);
        int status = H5Fclose(_file_id);
        if (status < 0) throw Exception(_file_name + ": error in H5Fclose");
        _file_id = 0;
    }

    std::string file_version() const
    {
        double v;
        hdf5_tools::Reader< double >()(_file_id, "/file_version", v);
        // convert it to string
        std::ostringstream os;
        os << v;
        return os.str();
    }

    std::string basecall_version() const
    {
        std::string res;
        hdf5_tools::Reader< std::string >()(_file_id, "/Analyses/Basecall_2D_000/version", res);
        return res;
    }

    std::string eventdetection_version() const
    {
        std::string res;
        hdf5_tools::Reader< std::string >()(_file_id, "/Analyses/EventDetection_000/version", res);
        return res;
    }

    std::string sequences_version() const
    {
        std::vector< std::string > tmp;
        hdf5_tools::Reader< std::string >()(_file_id, "/Sequences/Meta/version", tmp);
        std::string res;
        for (const auto& s: tmp)
        {
            res += s;
        }
        return res;
    }

    std::string basecalled_2D() const
    {
        std::string res;
        hdf5_tools::Reader< std::string >()(_file_id, "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq", res);
        
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
        hdf5_tools::Reader< Event_Alignment_Entry >()(_file_id, "/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment", res, &m);
        return res;
    }

    bool have_model(size_t i) const
    {
        return hdf5_tools::addr_exists(_file_id, model_path(i));
    }
    bool have_events(size_t i) const
    {
        return hdf5_tools::addr_exists(_file_id, events_path(i));
    }

    std::vector< Model_Entry > get_model(size_t i) const
    {
        std::vector< Model_Entry > res;
        hdf5_tools::Compound_Map m;
        m.add_member("kmer", &Model_Entry::kmer);
        m.add_member("level_mean", &Model_Entry::level_mean);
        m.add_member("level_stdv", &Model_Entry::level_stdv);
        hdf5_tools::Reader< Model_Entry >()(_file_id, model_path(i), res, &m);
        return res;
    }

    Model_Parameters get_model_parameters(size_t i) const
    {
        Model_Parameters res;
        std::string path = model_path(i);
        hdf5_tools::Reader< double >()(_file_id, path + "/drift", res.drift);
        hdf5_tools::Reader< double >()(_file_id, path + "/scale", res.scale);
        hdf5_tools::Reader< double >()(_file_id, path + "/scale_sd", res.scale_sd);
        hdf5_tools::Reader< double >()(_file_id, path + "/shift", res.shift);
        hdf5_tools::Reader< double >()(_file_id, path + "/var", res.var);
        hdf5_tools::Reader< double >()(_file_id, path + "/var_sd", res.var_sd);
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
        hdf5_tools::Reader< Event_Entry >()(_file_id, events_path(i), res, &m);
        return res;
    }

private:
    std::string _file_name;
    hid_t _file_id;

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

}; // class File

} // namespace fast5

#endif
