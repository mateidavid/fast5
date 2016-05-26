#include <cassert>
#include <iostream>
#include <string>

#include "hdf5_tools.hpp"

using namespace std;
using namespace hdf5;

struct A
{
    int val_1;
    float val_2;
    char val_3[6];
    array< char, 6 > val_4;
};

int main(int argc, char* argv[])
{
    if (argc != 2 and argc != 3)
    {
        cerr << "use: " << argv[0] << " [-f] <fast5_file>" << endl;
        return EXIT_FAILURE;
    }
    bool force = string(argv[1]) == "-f";
    string file_name(argv[force? 2 : 1]);
    {
        hdf5_tools::File f;
        //
        // All fast5 operations are performed inside a try-catch block. This should
        // resist various hdf5 errors without leaking memory.
        //
        try
        {
            //
            // create file; without -f, fail if it exist
            //
            f.create(file_name, force);
            assert(f.is_open());
            assert(f.is_rw());
            int val_1 = 42;
            float val_2 = 3.14;
            char val_3[6] = "ACGTA";
            array< char, 6 > val_4 = { "AACCG" };
            string val_5("CCCGG");
            static_assert(hdf5_tools::detail::mem_type_class< void >::value == 0, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_1) >::value == 1, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_2) >::value == 1, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_3) >::value == 2, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_4) >::value == 2, "");
            static_assert(hdf5_tools::detail::mem_type_class< decltype(val_5) >::value == 3, "");
            static_assert(hdf5_tools::detail::mem_type_class< std::true_type >::value == 4, "");
            //
            // write integer
            //
            f.write("/val_1", false, val_1);
            f.write("/val_1_as_64", false, val_1, H5T_STD_I64LE);
            f.write("/val_1_v", false, vector< int >(3, val_1));
            //
            // write float
            //
            f.write("/val_2", false, val_2);
            f.write("/val_2_as_64", false, val_2, H5T_IEEE_F64LE);
            f.write("/val_2_v", false, vector< float >(3, val_2));
            //
            // write fixlen string: char[]
            //
            f.write("/val_3", false, val_3);
            f.write("/val_3_as_len_3", false, val_3, 3);
            f.write("/val_3_as_varlen", false, val_3, -1);
            //
            // write fixlen string: std::array< char >
            //
            f.write("/val_4", false, val_4);
            f.write("/val_4_as_len_3", false, val_4, 3);
            f.write("/val_4_as_varlen", false, val_4, -1);
            f.write("/val_4_v", false, vector < decltype(val_4) >(3, val_4));
            f.write("/val_4_v_as_len_3", false, vector < decltype(val_4) >(3, val_4), 3);
            f.write("/val_4_v_as_varlen", false, vector < decltype(val_4) >(3, val_4), -1);
            //
            // write varlen string
            //
            f.write("/val_5", false, val_5);
            f.write("/val_5_as_len_3", false, val_5, 3);
            f.write("/val_5_as_fixlen", false, val_5, 0);
            f.write("/val_5_v", false, vector< decltype(val_5) >(3, val_5));
            f.write("/val_5_v_as_len_3", false, vector< decltype(val_5) >(3, val_5), 3);
            f.write("/val_5_v_as_fixlen", false, vector< decltype(val_5) >(1, val_5), 0); // only size 1
            //
            // write compound
            //
            A val_6{ 1, 3.14, "ACGTA", "CGTAC" };
            hdf5_tools::Compound_Map cm;
            cm.add_member("val_1", &A::val_1);
            cm.add_member("val_2", &A::val_2);
            cm.add_member("val_3", &A::val_3);
            cm.add_member("val_4", &A::val_4);
            f.write("/val_6", false, val_6, cm);
            f.write("/val_6_v", false, vector< A >(3, val_6), cm);
        }
        catch (hdf5_tools::Exception& e)
        {
            cout << "hdf5 error: " << e.what() << endl;
        }
        //
        // fast5 file is closed by its destructor at the end of this scope
        //
    }
    assert(hdf5_tools::File::get_object_count() == 0);
}
