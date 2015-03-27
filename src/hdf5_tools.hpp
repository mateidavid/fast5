#ifndef __HDF5_TOOLS_HPP
#define __HDF5_TOOLS_HPP

#include <cassert>
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace hdf5
{
#include <hdf5.h>
}

namespace hdf5_tools
{
using namespace hdf5;

/// Exception class thrown by failed hdf5 operations.
class Exception
    : public std::exception
{
public:
    Exception(const std::string& msg) : _msg(msg) {}
    const char* what() const noexcept { return _msg.c_str(); }
private:
    std::string _msg;
}; // class Exception

// Forward declaration
class Compound_Map;

namespace detail
{

/// TempMetaFunc: Given destination type, deduce memory type to be used in hdf5 read operation.
/// Only useful for numeric types.
/// HDF5 idiosyncracy:
///   Types such as H5T_NATIVE_INT are not constants(!?), so id() is not a constexpr.
template < typename T > struct get_mem_type           { static hid_t id() { return -1;                 } };
// signed integral:
template <> struct get_mem_type< char >               { static hid_t id() { return H5T_NATIVE_CHAR;    } };
template <> struct get_mem_type< short >              { static hid_t id() { return H5T_NATIVE_SHORT;   } };
template <> struct get_mem_type< int >                { static hid_t id() { return H5T_NATIVE_INT;     } };
template <> struct get_mem_type< long >               { static hid_t id() { return H5T_NATIVE_LONG;    } };
template <> struct get_mem_type< long long >          { static hid_t id() { return H5T_NATIVE_LLONG;   } };
// unsigned integral:
template <> struct get_mem_type< unsigned char >      { static hid_t id() { return H5T_NATIVE_UCHAR;   } };
template <> struct get_mem_type< unsigned short >     { static hid_t id() { return H5T_NATIVE_USHORT;  } };
template <> struct get_mem_type< unsigned >           { static hid_t id() { return H5T_NATIVE_UINT;    } };
template <> struct get_mem_type< unsigned long >      { static hid_t id() { return H5T_NATIVE_ULONG;   } };
template <> struct get_mem_type< unsigned long long > { static hid_t id() { return H5T_NATIVE_ULLONG;  } };
// float:
template <> struct get_mem_type< float >              { static hid_t id() { return H5T_NATIVE_FLOAT;   } };
template <> struct get_mem_type< double >             { static hid_t id() { return H5T_NATIVE_DOUBLE;  } };
template <> struct get_mem_type< long double >        { static hid_t id() { return H5T_NATIVE_LDOUBLE; } };

/// TempMetaFunc: Given destination type, can we read it
template < typename Out_Data_Type >
struct can_read
{
    static const bool value =
        std::is_integral< Out_Data_Type >::value
        or std::is_floating_point< Out_Data_Type >::value
        or std::is_same< typename std::remove_extent< Out_Data_Type >::type, char >::value
        or std::is_same< Out_Data_Type, std::string >::value
        or std::is_class< Out_Data_Type >:: value;
};

/// TempMetaFunc: Given a destination type, does it need a compound map
template < typename Out_Data_Type >
struct read_as_atomic
{
    static const bool value =
        std::is_integral< Out_Data_Type >::value
        or std::is_floating_point< Out_Data_Type >::value
        or std::is_same< typename std::remove_extent< Out_Data_Type >::type, char >::value
        or std::is_same< Out_Data_Type, std::string >::value;
};

/// Compute offset of a struct member from a member pointer (runtime version).
template < typename T, typename U >
std::size_t offset_of(U T::* mem_ptr)
{
    return reinterpret_cast< std::size_t >(&(((T*)0)->*mem_ptr));
}

/// Description of a member inside a compound
/// Only works with numeric, string, and struct types.
struct Compound_Member_Description
{
public:
    Compound_Member_Description(const std::string& _name, size_t _offset, hid_t _numeric_type_id)
        : name(_name), offset(_offset), numeric_type_id(_numeric_type_id)
    {
        type = numeric;
    }
    Compound_Member_Description(const std::string& _name, size_t _offset, size_t _char_array_size)
        : name(_name), offset(_offset), char_array_size(_char_array_size)
    {
        type = char_array;
    }
    Compound_Member_Description(const std::string& _name, size_t _offset)
        : name(_name), offset(_offset)
    {
        type = string;
    }
    Compound_Member_Description(const std::string& _name, size_t _offset, const Compound_Map* _compound_map_ptr)
        : name(_name), offset(_offset), compound_map_ptr(_compound_map_ptr)
    {
        type = compound;
    }

    bool is_numeric() const { return type == numeric; }
    bool is_char_array() const { return type == char_array; }
    bool is_string() const { return type == string; }
    bool is_compound() const { return type == compound; }

    std::string name;
    size_t offset;
    union
    {
        hid_t numeric_type_id;
        size_t char_array_size;
        const Compound_Map* compound_map_ptr;
    };

private:
    enum member_type
    {
        numeric,
        char_array,
        string,
        compound
    };
    member_type type;
}; // Compound_Member_Description

} // namespace detail

/// A map of struct fields to tags that is used to read compound datatypes
class Compound_Map
{
public:
    Compound_Map() = default;
    Compound_Map(const Compound_Map&) = delete;
    Compound_Map(Compound_Map&&) = default;
    Compound_Map& operator = (const Compound_Map&) = delete;
    Compound_Map& operator = (Compound_Map&&) = default;
    ~Compound_Map() = default;

    template < typename T, typename U >
    void add_member(const std::string& name, U T::* mem_ptr)
    {
        static_assert(std::is_integral< U >::value
                      or std::is_floating_point< U >::value
                      or std::is_same< typename std::remove_extent< U >::type, char >::value
                      or std::is_same< U, std::string >::value,
                      "add_member(name, mem_ptr) overload expects numerical or string types only ");
        if (std::is_integral< U >::value or std::is_floating_point< U >::value)
        {
            _members.emplace_back(name, detail::offset_of(mem_ptr), detail::get_mem_type< U >::id());
        }
        else if (std::is_same< typename std::remove_extent< U >::type, char >::value)
        {
            _members.emplace_back(name, detail::offset_of(mem_ptr), sizeof(U));
        }
        else if (std::is_same< U, std::string >::value)
        {
            _members.emplace_back(name, detail::offset_of(mem_ptr));
        }
    }

    template < typename T, typename U >
    void add_member(const std::string& name, U T::* mem_ptr, const Compound_Map* compound_map_ptr)
    {
        assert(false); // not currently implemented
        static_assert(std::is_class< U >::value,
                      "add_member(name, mem_ptr, compound_map_ptr) overload expects class types only ");
        _members.emplace_back(name, detail::offset_of(mem_ptr), compound_map_ptr);
    }

    const std::vector< detail::Compound_Member_Description >& members() const { return _members; }

private:
    std::vector< detail::Compound_Member_Description > _members;
}; // Compound_Map

namespace detail
{

// TempSpec: reading numerics
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Extent_Atomic_Reader
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      const Compound_Map*, hid_t obj_id, hid_t,
                      const std::string&, std::function< hid_t(hid_t) >,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        hid_t mem_type_id = get_mem_type< Out_Data_Type >::id();
        assert(mem_type_id != -1);
        int status = read_fcn(obj_id, mem_type_id, static_cast< void* >(dest.data()));
        if (status < 0) throw Exception(loc_full_name + ": error in " + read_fcn_name);
    }
}; // struct Extent_Atomic_Reader

// TempSpec: for reading strings
template < typename Out_Data_Storage >
struct Extent_Atomic_Reader< std::string, Out_Data_Storage >
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      const Compound_Map*, hid_t obj_id, hid_t obj_space_id,
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
}; // struct Extent_Atomic_Reader< std::string >

template < typename Out_Data_Type, typename Out_Data_Storage >
struct Extent_Compound_Reader
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      const Compound_Map* compound_map_ptr, hid_t obj_id, hid_t obj_space_id,
                      const std::string& get_type_fcn_name, std::function< hid_t(hid_t) > get_type_fcn,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        int status;
        assert(compound_map_ptr);
        hid_t file_type_id = get_type_fcn(obj_id);
        if (file_type_id < 0) throw Exception(loc_full_name + ": error in " + get_type_fcn_name);
        H5T_class_t file_type_class = H5Tget_class(file_type_id);
        if (file_type_class == H5T_NO_CLASS) throw Exception(loc_full_name + ": error in H5Tget_class(file_type)");
        if (file_type_class != H5T_COMPOUND) throw Exception(loc_full_name + ": expected H5T_COMPOUND datatype");

        // pass 1
        //   read numeric and char_array members only
        hid_t mem_type_id = H5Tcreate(H5T_COMPOUND, sizeof(Out_Data_Type));
        std::vector< hid_t > mem_stype_id_v;
        for (const auto& e : compound_map_ptr->members())
        {
            assert(not e.is_compound()); // not implemented yet
            if (e.is_string()) continue;
            int file_stype_idx = H5Tget_member_index(file_type_id, e.name.c_str());
            if (file_stype_idx < 0) throw Exception(loc_full_name + ": missing member \"" + e.name + "\"");
            hid_t file_stype_id = H5Tget_member_type(file_type_id, file_stype_idx);
            if (file_stype_id < 0) throw Exception(loc_full_name + ": error in H5Tget_member_type");
            H5T_class_t file_stype_class = H5Tget_class(file_stype_id);
            if (file_stype_class == H5T_NO_CLASS) throw Exception(loc_full_name + ": error in H5Tget_class(file_stype)");
            if (e.is_numeric())
            {
                if (file_stype_class != H5T_INTEGER and file_stype_class != H5T_FLOAT)
                    throw Exception(loc_full_name + ": member \"" + e.name + "\" is numeric, but file_stype is not numeric");
                status = H5Tinsert(mem_type_id, e.name.c_str(), e.offset, e.numeric_type_id);
                if (status < 0) throw Exception(loc_full_name + ": error in H5Tinsert(\"" + e.name + "\")");
            }
            if (e.is_char_array())
            {
                if (file_stype_class != H5T_STRING)
                    throw Exception(loc_full_name + ": member \"" + e.name + "\" is char_array, but file_stype is not H5T_STRING");
                status = H5Tis_variable_str(file_stype_id);
                if (status < 0) throw Exception(loc_full_name + ": error in H5Tis_variable_str(\"" + e.name + "\")");
                if (status) throw Exception(loc_full_name + ": member \"" + e.name + "\" is a char_array, but file_stype is a variable len string");
                //size_t file_stype_size = H5Tget_size(file_stype_id);
                //if (file_stype_size == 0) throw Exception(loc_full_name + ": H5Tget_size(\"" + e.name + "\") returned 0");
                hid_t mem_stype_id = H5Tcopy(H5T_C_S1);
                if (mem_stype_id < 0) throw Exception(loc_full_name + ": member \"" + e.name + "\": error in H5Tcopy");
                status = H5Tset_size(mem_stype_id, e.char_array_size);
                if (status < 0) throw Exception(loc_full_name + ": error in H5Tset_size(\"" + e.name + "\")");
                status = H5Tinsert(mem_type_id, e.name.c_str(), e.offset, mem_stype_id);
                if (status < 0) throw Exception(loc_full_name + ": error in H5Tinsert(\"" + e.name + "\")");
                mem_stype_id_v.push_back(mem_stype_id);
            }
            status = H5Tclose(file_stype_id);
            if (status < 0) throw Exception(loc_full_name + ": member \"" + e.name + "\": error in H5Tclose(file_stype)");
        }
        // perform the actual read
        status = read_fcn(obj_id, mem_type_id, static_cast< void* >(dest.data()));
        if (status < 0) throw Exception(loc_full_name + ": pass 1: error in " + read_fcn_name);
        // release the memory types
        for (const auto& mem_stype_id : mem_stype_id_v)
        {
            status = H5Tclose(mem_stype_id);
            if (status < 0) throw Exception(loc_full_name + ": error in H5Tclose(mem_stype)");
        }
        mem_stype_id_v.clear();
        status = H5Tclose(mem_type_id);
        if (status < 0) throw Exception(loc_full_name + ": error in H5Tclose(mem_type)");

        // pass 2
        //   read strings
        for (const auto& e : compound_map_ptr->members())
        {
            assert(not e.is_compound()); // not implemented yet
            if (e.is_numeric() or e.is_char_array()) continue;
            //TODO
            assert(false);
        }
    }
}; //struct Extent_Compound_Reader

// TempSpec: read extent of atomic types
template < typename Out_Data_Type, typename Out_Data_Storage, bool = true >
struct Extent_Reader_as_atomic
    : Extent_Atomic_Reader< Out_Data_Type, Out_Data_Storage >
{};

// TempSpec: read extent of compound types
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Extent_Reader_as_atomic< Out_Data_Type, Out_Data_Storage, false >
    : Extent_Compound_Reader< Out_Data_Type, Out_Data_Storage >
{};

// branch on atomic/compound destination
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Extent_Reader
    : public Extent_Reader_as_atomic< Out_Data_Type, Out_Data_Storage, read_as_atomic< Out_Data_Type >::value >
{};

template < typename, typename, bool >
struct Object_Reader_impl;

// TempSpec: reading scalars
template < typename Out_Data_Type >
struct Object_Reader_impl< Out_Data_Type, Out_Data_Type, true >
{
    void operator () (const std::string& loc_full_name, Out_Data_Type& dest,
                      const Compound_Map* compound_map_ptr, hid_t obj_id, hid_t obj_space_id,
                      const std::string& get_type_fcn_name, std::function< hid_t(hid_t) > get_type_fcn,
                      const std::string& read_fcn_name, std::function< herr_t(hid_t, hid_t, void*) > read_fcn)
    {
        H5S_class_t obj_class_t = H5Sget_simple_extent_type(obj_space_id);
        if (obj_class_t == H5S_NO_CLASS) throw Exception(loc_full_name + ": error in H5Sget_simple_extent_type");
        if (obj_class_t != H5S_SCALAR)
            throw Exception(loc_full_name + ": reading as scalar, but dataspace not H5S_SCALAR");
        std::vector< Out_Data_Type > tmp(1);
        Extent_Reader< Out_Data_Type, std::vector< Out_Data_Type > >()(
            loc_full_name, tmp, compound_map_ptr, obj_id, obj_space_id,
            get_type_fcn_name, get_type_fcn,
            read_fcn_name, read_fcn);
        dest = std::move(tmp[0]);
    }
};

// TempSpec: reading vectors
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Object_Reader_impl< Out_Data_Type, Out_Data_Storage, false >
{
    void operator () (const std::string& loc_full_name, Out_Data_Storage& dest,
                      const Compound_Map* compound_map_ptr, hid_t obj_id, hid_t obj_space_id,
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
            loc_full_name, dest, compound_map_ptr, obj_id, obj_space_id,
            get_type_fcn_name, get_type_fcn,
            read_fcn_name, read_fcn);
    }
};

// TempMetaFunc: split scalar & vector reading branches
template < typename Out_Data_Type, typename Out_Data_Storage >
struct Object_Reader
    : public Object_Reader_impl< Out_Data_Type, Out_Data_Storage, std::is_same< Out_Data_Type, Out_Data_Storage >::value > {};

// open object and object space, then delegate
template < typename Out_Data_Type, typename Out_Data_Storage >
void read_obj_helper(const std::string& loc_full_name, Out_Data_Storage& dest, const Compound_Map* compound_map_ptr,
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
        loc_full_name, dest, compound_map_ptr, obj_id, obj_space_id,
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
void read_addr(hid_t root_id, const std::string& loc_path, const std::string& loc_name,
               Out_Data_Storage& dest, const Compound_Map* compound_map_ptr)
{
    assert(root_id > 0);
    std::string loc_full_name = loc_path + loc_name;
    // determine if object is an attribute; otherwise, assume it's a dataset
    int status;
    status = H5Aexists_by_name(root_id, loc_path.c_str(), loc_name.c_str(), H5P_DEFAULT);
    if (status < 0) throw Exception(loc_full_name + ": error in H5Aexists_by_name");
    bool is_attr = status > 0;
    if (is_attr)
    {
        read_obj_helper< Out_Data_Type, Out_Data_Storage >(
            loc_full_name, dest, compound_map_ptr,
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
            loc_full_name, dest, compound_map_ptr,
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

// TempSpec: for atomic types
template < typename Out_Data_Type, bool = true >
struct Reader_as_atomic
{
    template < typename Out_Data_Storage >
    void operator () (hid_t root_id, const std::string& loc_path, const std::string& loc_name,
                      Out_Data_Storage& dest)
    {
        static_assert(can_read< Out_Data_Type >::value,
                      "Reader_impl<Out_Data_Type,true>: expected a readable destination");
        static_assert(read_as_atomic< Out_Data_Type >::value,
                      "Reader_impl<Out_Data_Type,true>: expected a type readable as atomic");
        read_addr< Out_Data_Type, Out_Data_Storage >(root_id, loc_path, loc_name, dest, nullptr);
    }
};

// TempSpec: for compound types
template < typename Out_Data_Type >
struct Reader_as_atomic< Out_Data_Type, false >
{
    template < typename Out_Data_Storage >
    void operator () (hid_t root_id, const std::string& loc_path, const std::string& loc_name,
                      Out_Data_Storage& dest, const Compound_Map* compound_map_ptr)
    {
        static_assert(can_read< Out_Data_Type >::value,
                      "Reader_impl<Out_Data_Type,false>: expected a readable destination");
        static_assert(not read_as_atomic< Out_Data_Type >::value,
                      "Reader_impl<Out_Data_Type,false>: expected a type readable as compound");
        read_addr< Out_Data_Type, Out_Data_Storage >(root_id, loc_path, loc_name, dest, compound_map_ptr);
    }
};

template < typename Out_Data_Type >
struct Reader : public Reader_as_atomic< Out_Data_Type, read_as_atomic< Out_Data_Type >::value >
{};

} // namespace detail

/// An HDF5 file reader
class File_Reader
{
public:
    File_Reader() : _file_id(0) {}
    File_Reader(const std::string& file_name) : _file_id(0) { open(file_name); }
    File_Reader(const File_Reader&) = delete;
    File_Reader& operator = (const File_Reader&) = delete;
    ~File_Reader() { if (is_open()) close(); }

    bool is_open() const { return _file_id > 0; }
    const std::string& file_name() const { return _file_name; }

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
        _file_name.clear();
    }

    /// Determine if address is an attribute or dataset
    bool exists(const std::string& loc_full_name) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name[0] == '/');
        std::string loc_path;
        std::string loc_name;
        std::tie(loc_path, loc_name) = split_full_name(loc_full_name);
        int status;
        // check all path elements exist, except for what is to the right of the last '/'
        size_t pos = 0;
        while (true)
        {
            ++pos;
            pos = loc_full_name.find('/', pos);
            if (pos == std::string::npos) break;
            std::string tmp = loc_full_name.substr(0, pos);
            status = H5Lexists(_file_id, tmp.c_str(), H5P_DEFAULT);
            if (status < 0) throw Exception(loc_full_name + ": error in H5Lexists");
            if (not status) return false;
            status = H5Oexists_by_name(_file_id, tmp.c_str(), H5P_DEFAULT);
            if (status < 0) throw Exception(loc_full_name + ": error in H5Oexists_by_name");
            if (not status) return false;
        }
        status = H5Aexists_by_name(_file_id, loc_path.c_str(), loc_name.c_str(), H5P_DEFAULT);
        if (status < 0) throw Exception(loc_full_name + ": error in H5Aexists_by_name");
        if (status) return true;
        // not an attribute: try to open as a dataset
        hid_t ds_id = H5Dopen(_file_id, loc_full_name.c_str(), H5P_DEFAULT);
        if (ds_id < 0) return false;
        status = H5Dclose(ds_id);
        if (status < 0) throw Exception(loc_full_name + ": error in H5Dclose");
        return true;
    }
    /// Read attribute or dataset at address
    template < typename Out_Data_Type, typename ...Args >
    void read(const std::string& loc_full_name, Args&& ...args) const
    {
        assert(is_open());
        assert(not loc_full_name.empty() and loc_full_name[0] == '/');
        std::string loc_path;
        std::string loc_name;
        std::tie(loc_path, loc_name) = split_full_name(loc_full_name);
        detail::Reader< Out_Data_Type >()(_file_id, loc_path, loc_name, std::forward< Args >(args)...);
    }

private:
    std::string _file_name;
    hid_t _file_id;

    /// Split a full name into path and name
    static std::pair< std::string, std::string > split_full_name(const std::string& full_name)
    {
        auto last_slash_pos = full_name.find_last_of('/');
        std::string path = last_slash_pos != std::string::npos? full_name.substr(0, last_slash_pos + 1) : std::string();
        std::string name = last_slash_pos != std::string::npos? full_name.substr(last_slash_pos + 1) : full_name;
        return std::make_pair(path, name);
    } // split_full_name
}; // class File_Reader

} // namespace hdf5_tools

#endif
