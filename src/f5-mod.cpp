#include <cassert>
#include <iostream>
#include <string>

#include "fast5.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "use: " << argv[0] << " <fast5_file>" << endl;
        return EXIT_FAILURE;
    }
    string file_name(argv[1]);
    {
        fast5::File f;
        //
        // All fast5 operations are performed inside a try-catch block. This should
        // resist various hdf5 errors without leaking memory.
        //
        try
        {
            //
            // open file in rw mode
            //
            f.open(file_name, true);
            assert(f.is_open());
            assert(f.is_rw());
            //
            // find next available basecall group with given prefix
            //
            string test_bc_grp_prefix = "Test_";
            auto bc_grp_l = f.get_basecall_group_list();
            set< string > test_bc_grp_suffix_s;
            for (const auto& bc_grp : bc_grp_l)
            {
                if (bc_grp.compare(0, test_bc_grp_prefix.size(), test_bc_grp_prefix) == 0)
                {
                    cerr << "found group: " << test_bc_grp_prefix + bc_grp.substr(test_bc_grp_prefix.size()) << endl;
                }
                test_bc_grp_suffix_s.insert(bc_grp.substr(test_bc_grp_prefix.size()));
            }
            string test_bc_grp_suffix;
            for (unsigned i = 0; i < 1000; ++i)
            {
                ostringstream os;
                os << setw(3) << setfill('0') << i;
                if (test_bc_grp_suffix_s.count(os.str()) == 0)
                {
                    test_bc_grp_suffix = os.str();
                    break;
                }
            }
            assert(not test_bc_grp_suffix.empty());
            clog << "using group: " << test_bc_grp_prefix + test_bc_grp_suffix << endl;
            //
            // add basecall seq
            //
            f.add_basecall_seq(test_bc_grp_prefix + test_bc_grp_suffix, 0, "test_name", "ACGT");
            //
            // add basecall events
            //
            vector< fast5::Event_Entry > ev(3, {55.0, 1.0, 0.05, 0.01, .5, .5, .7, .1, .1, .1, 0,
                        array< char, 8 >{"ACGTA"}, array< char, 8 >{"CGTAC"}});
            f.add_basecall_events(test_bc_grp_prefix + test_bc_grp_suffix, 0, ev);
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
