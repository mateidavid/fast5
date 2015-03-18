#ifndef __HDF5_TOOLS_HPP
#define __HDF5_TOOLS_HPP

#include <cassert>
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

namespace hdf5
{
#include <hdf5.h>
}

namespace hdf5_tools
{
using namespace hdf5;

class Exception
    : public std::exception
{
public:
    Exception(const std::string& msg) : _msg(msg) {}
    const char* what() const noexcept { return _msg.c_str(); }
private:
    std::string _msg;
}; // class Exception

namespace detail
{

/// Template meta-function for deducing memory type to be used during read.
/// `cls`: memory type class
/// `id()`: for H5T_INTEGER and H5T_FLOAT, the memory type to use
template < typename T > struct get_mem_type           { static const hid_t cls = H5T_NO_CLASS; static hid_t id() { return -1;                } };
// signed:
template <> struct get_mem_type< char >               { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_CHAR;    } };
template <> struct get_mem_type< short >              { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_SHORT;   } };
template <> struct get_mem_type< int >                { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_INT;     } };
template <> struct get_mem_type< long >               { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_LONG;    } };
template <> struct get_mem_type< long long >          { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_LLONG;   } };
// unsigned:
template <> struct get_mem_type< unsigned char >      { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_UCHAR;   } };
template <> struct get_mem_type< unsigned short >     { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_USHORT;  } };
template <> struct get_mem_type< unsigned >           { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_UINT;    } };
template <> struct get_mem_type< unsigned long >      { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_ULONG;   } };
template <> struct get_mem_type< unsigned long long > { static const hid_t cls = H5T_INTEGER; static hid_t id() { return H5T_NATIVE_ULLONG;  } };
// float:
template <> struct get_mem_type< float >              { static const hid_t cls = H5T_FLOAT;   static hid_t id() { return H5T_NATIVE_FLOAT;   } };
template <> struct get_mem_type< double >             { static const hid_t cls = H5T_FLOAT;   static hid_t id() { return H5T_NATIVE_DOUBLE;  } };
template <> struct get_mem_type< long double >        { static const hid_t cls = H5T_FLOAT;   static hid_t id() { return H5T_NATIVE_LDOUBLE; } };
// string:
template <> struct get_mem_type< std::string >        { static const hid_t cls = H5T_STRING;  static hid_t id() { return H5T_C_S1;           } };

std::pair< std::string, std::string > get_path_name(const std::string& full_name)
{
    auto last_slash_pos = full_name.find_last_of('/');
    std::string path = last_slash_pos != std::string::npos? full_name.substr(0, last_slash_pos + 1) : std::string();
    std::string name = last_slash_pos != std::string::npos? full_name.substr(last_slash_pos + 1) : full_name;
    return std::make_pair(path, name);
} // get_path_name

// template specialization for reading non-strings
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Extent_Reader
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      hid_t mem_type_id, hid_t obj_id, hid_t,
                      const std::string&, std::function< hid_t(hid_t) >,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        int status = read_fcn(obj_id, mem_type_id, static_cast< void* >(dest.data()));
        if (status < 0) throw Exception(loc_full_name + ": error in " + read_fcn_name);
    }
};

// template specialization for reading strings
template < typename Out_Data_Storage >
struct Extent_Reader< std::string, Out_Data_Storage >
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      hid_t, hid_t obj_id, hid_t obj_space_id,
                      const std::string& get_type_fcn_name, std::function< hid_t(hid_t) > get_type_fcn,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        int status;
        int file_type_id = get_type_fcn(obj_id);
        if (file_type_id < 0) throw Exception(loc_full_name + ": error in " + get_type_fcn_name);
        int is_vlen_str = H5Tis_variable_str(file_type_id);
        if (is_vlen_str < 0) throw Exception(loc_full_name + ": error in H5Tis_variable_str");
        hid_t mem_type_id = H5Tcopy(H5T_C_S1);
        if (mem_type_id < 0) throw Exception(loc_full_name + ": error in H5Tcopy");
        if (is_vlen_str) // stored as variable-length string
        {
            // compute mem_type
            status = H5Tset_size(mem_type_id, H5T_VARIABLE);
            if (status < 0) throw Exception(loc_full_name + ": error in H5Tset_size(variable)");
            // prepare buffer to receive data
            std::vector< char* > char_p_buff(dest.size(), nullptr);
            // perform the read
            status = read_fcn(obj_id, mem_type_id, static_cast< void* >(char_p_buff.data()));
            if (status < 0) throw Exception(loc_full_name + ": error in " + read_fcn_name);
            // transfer strings to destination
            for (size_t i = 0; i < dest.size(); ++i)
            {
                if (not char_p_buff[i]) throw Exception(loc_full_name + ": " + read_fcn_name + " did not fill buffer");
                dest[i] = char_p_buff[i];
            }
            // reclaim memory allocated by libhdf5
            status = H5Dvlen_reclaim(mem_type_id, obj_space_id, H5P_DEFAULT, char_p_buff.data());
            if (status < 0) throw Exception(loc_full_name + ": error in H5Dvlen_reclaim");
        }
        else // stored as fixed-length string
        {
            // compute mem_type
            size_t sz = H5Tget_size(file_type_id);
            if (sz == 0) throw Exception(loc_full_name + ": H5Tget_size returned 0; is this an error?!");
            status = H5Tset_size(mem_type_id, sz + 1);
            if (status < 0) throw Exception(loc_full_name + ": error in H5Tset_size(fixed)");
            // prepare buffer to receieve data
            std::vector< char > char_buff(dest.size() * (sz + 1));
            // perform the read
            status = read_fcn(obj_id, mem_type_id, static_cast< void* >(char_buff.data()));
            if (status < 0) throw Exception(loc_full_name + ": error in " + read_fcn_name);
            // transfer strings to destination
            for (size_t i = 0; i < dest.size(); ++i)
            {
                dest[i] = std::string(&char_buff[i * (sz + 1)], sz);
            }
        }
        status = H5Tclose(mem_type_id);
        if (status < 0) throw Exception(loc_full_name + ": error in H5Tclose(mem_type_id)");
        status = H5Tclose(file_type_id);
        if (status < 0) throw Exception(loc_full_name + ": error in H5Tclose(file_type_id)");
    }
};

template < typename, typename, bool >
struct Object_Reader_impl;

// template specialization for reading scalars
template < typename Out_Data_Type >
struct Object_Reader_impl< Out_Data_Type, Out_Data_Type, true >
{
    void operator () (const std::string& loc_full_name, Out_Data_Type& dest,
                      hid_t mem_type_id, hid_t obj_id, hid_t obj_space_id,
                      const std::string& get_type_fcn_name, std::function< hid_t(hid_t) > get_type_fcn,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        H5S_class_t obj_class_t = H5Sget_simple_extent_type(obj_space_id);
        if (obj_class_t == H5S_NO_CLASS) throw Exception(loc_full_name + ": error in H5Sget_simple_extent_type");
        if (obj_class_t != H5S_SCALAR)
            throw Exception(loc_full_name + ": reading as scalar, but dataspace not H5S_SCALAR");
        std::vector< Out_Data_Type > tmp(1);
        Extent_Reader< Out_Data_Type, std::vector< Out_Data_Type > >()(
            loc_full_name, tmp, mem_type_id, obj_id, obj_space_id,
            get_type_fcn_name, get_type_fcn,
            read_fcn_name, read_fcn);
        dest = std::move(tmp[0]);
    }
};

// template specialization for reading vectors
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Object_Reader_impl< Out_Data_Type, Out_Data_Storage, false >
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      hid_t mem_type_id, hid_t obj_id, hid_t obj_space_id,
                      const std::string& get_type_fcn_name, std::function< hid_t(hid_t) > get_type_fcn,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        H5S_class_t obj_class_t = H5Sget_simple_extent_type(obj_space_id);
        if (obj_class_t == H5S_NO_CLASS) throw Exception(loc_full_name + ": error in H5Sget_simple_extent_type");
        if (obj_class_t != H5S_SIMPLE)
            throw Exception(loc_full_name + ": reading as vector, but dataspace not H5S_SIMPLE");
        int status = H5Sget_simple_extent_dims(obj_space_id, nullptr, nullptr);
        if (status < 0) throw Exception(loc_full_name + ": error in H5Sget_simple_extent_dims");
        if (status != 1) throw Exception(loc_full_name + ": expected extent of dimension 1");
        hsize_t sz;
        H5Sget_simple_extent_dims(obj_space_id, &sz, nullptr);
        dest.clear();
        dest.resize(sz);
        Extent_Reader< Out_Data_Type, Out_Data_Storage >()(
            loc_full_name, dest, mem_type_id, obj_id, obj_space_id,
            get_type_fcn_name, get_type_fcn,
            read_fcn_name, read_fcn);
    }
};

// TMF: split scalar & vector reading branches
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Object_Reader
    : public Object_Reader_impl< Out_Data_Type, Out_Data_Storage, std::is_same< Out_Data_Type, Out_Data_Storage >::value > {};

// open object and object space, then delegate
template < typename Out_Data_Type, typename Out_Data_Storage >
void read_obj_helper(const std::string& loc_full_name, Out_Data_Storage& dest, hid_t mem_type_id,
                     const std::string& open_fcn_name, std::function< hid_t(void) > open_fcn,
                     const std::string& close_fcn_name, std::function< herr_t(hid_t) > close_fcn,
                     const std::string& get_space_fcn_name, std::function< hid_t(hid_t) > get_space_fcn,
                     const std::string& get_type_fcn_name, std::function< hid_t(hid_t) > get_type_fcn,
                     const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
{
    int status;
    // open object
    hid_t obj_id = open_fcn();
    if (obj_id < 0) throw Exception(loc_full_name + ": error in " + open_fcn_name);
    // open object space, check reading ode matches storage mode (scalar/vector)
    hid_t obj_space_id = get_space_fcn(obj_id);
    if (obj_space_id < 0) throw Exception(loc_full_name + ": error in " + get_space_fcn_name);
    // read object
    Object_Reader< Out_Data_Type, Out_Data_Storage >()(
        loc_full_name, dest, mem_type_id, obj_id, obj_space_id,
        get_type_fcn_name, get_type_fcn,
        read_fcn_name, read_fcn);
    // close object space & object
    status = H5Sclose(obj_space_id);
    if (status < 0) throw Exception(loc_full_name + ": error in H5Sclose");
    status = close_fcn(obj_id);
    if (status < 0) throw Exception(loc_full_name + ": error in " + close_fcn_name);
}

// determine if address is attribute or dataset, then delegate
template < typename Out_Data_Type, typename Out_Data_Storage >
void read_addr(hid_t root_id, const std::string& loc_full_name,
               Out_Data_Storage& dest, hid_t mem_type_id)
{
    assert(root_id > 0);
    std::string loc_path;
    std::string loc_name;
    std::tie(loc_path, loc_name) = get_path_name(loc_full_name);
    // determine if object is an attribute; otherwise, assume it's a dataset
    int status;
    status = H5Aexists_by_name(root_id, loc_path.c_str(), loc_name.c_str(), H5P_DEFAULT);
    if (status < 0) throw Exception(loc_full_name + ": error in H5Aexists_by_name");
    bool is_attr = status > 0;
    if (is_attr)
    {
        read_obj_helper< Out_Data_Type, Out_Data_Storage >(
            loc_full_name, dest, mem_type_id,
            "H5Aopen_by_name",
            [&root_id, &loc_path, &loc_name] ()
            {
                return H5Aopen_by_name(root_id, loc_path.c_str(), loc_name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
            },
            "H5Aclose", &H5Aclose,
            "H5Aget_space", &H5Aget_space,
            "H5Aget_type", &H5Aget_type,
            "H5Aread",
            [] (hid_t id, hid_t mem_type_id, void* dest_p)
            {
                return H5Aread(id, mem_type_id, dest_p);
            });
    }
    else
    {
        read_obj_helper< Out_Data_Type, Out_Data_Storage >(
            loc_full_name, dest, mem_type_id,
            "H5Dopen",
            [&root_id, &loc_full_name] ()
            {
                return H5Dopen(root_id, loc_full_name.c_str(), H5P_DEFAULT);
            },
            "H5Dclose", &H5Dclose,
            "H5Dget_space", &H5Dget_space,
            "H5Dget_type", &H5Dget_type,
            "H5Dread",
            [] (hid_t id, hid_t mem_type_id, void* dest_p)
            {
                return H5Dread(id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dest_p);
            });
    }
} // read_addr

} // namespace detail

template < typename Out_Data_Type >
struct Reader
{
    template < typename Out_Data_Storage >
    void operator () (hid_t root_id, const std::string& loc_full_name, Out_Data_Storage& dest, hid_t mem_type_id = -1)
    {
        if (detail::get_mem_type< Out_Data_Type >::cls == H5T_INTEGER
            or detail::get_mem_type< Out_Data_Type >::cls == H5T_FLOAT
            or detail::get_mem_type< Out_Data_Type >::cls == H5T_STRING)
        {
            // HDF5 idiosyncracy:
            //   This check cannot be a static because types such as
            //   H5T_NATIVE_INT are not constants, so id() is not a constexpr.
            assert(mem_type_id == -1);
            mem_type_id = detail::get_mem_type< Out_Data_Type >::id();
        }
        detail::read_addr< Out_Data_Type, Out_Data_Storage >(root_id, loc_full_name, dest, mem_type_id);
    }
}; // struct Reader

} // namespace hdf5

#endif
