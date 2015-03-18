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

private:
    std::string _file_name;
    hid_t _file_id;
}; // class File

} // namespace fast5

#endif
