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

    // open the FAST5 file for reading
    if (not hdf5_tools::File_Reader::is_valid_file(file_name))
    {
        cout << "not a hdf file [" << file_name << "]" << endl;
        return EXIT_SUCCESS;
    }
    if (not fast5::File::is_valid_file(file_name))
    {
        cout << "not a fast5 file [" << file_name << "]" << endl;
        return EXIT_SUCCESS;
    }

    {
        fast5::File f;
        //
        // All fast5 operations are performed inside a try-catch block. This should
        // resist various hdf5 errors without leaking memory.
        //
        try
        {
            //
            // open file
            //
            f.open(file_name);
            assert(f.is_open());

            //
            // extract version information for the ONT software used to generate this dataset
            //
            cout << "file_version=" << f.file_version() << endl;
            cout << "sampling_rate=" << f.get_sampling_rate() << endl;

            //
            // inspect sequences group
            //
            bool have_sequences_group = f.have_sequences_group();
            cout << "have_sequences_group=" << have_sequences_group << endl;
            if (have_sequences_group)
            {
                cout << "sequences_version=" << f.sequences_version() << endl;
            }

            //
            // inspect raw samples
            //
            bool have_raw_samples = f.have_raw_samples();
            cout << "have_raw_samples=" << have_raw_samples << endl;
            if (have_raw_samples)
            {
                auto r = f.get_raw_samples();
                cout << "number_raw_samples=" << r.size() << endl;
            }

            //
            // inspect eventdetection group
            //
            bool have_eventdetection_group = f.have_eventdetection_group();
            cout << "have_eventdetection_group=" << have_eventdetection_group << endl;
            if (have_eventdetection_group)
            {
                cout << "eventdetection_version=" << f.eventdetection_version() << endl;
                bool have_eventdetection_events = f.have_eventdetection_events();
                cout << "have_eventdetection_events=" << have_eventdetection_events << endl;
                if (have_eventdetection_events)
                {
                    cout << "eventdetection_read_name=" << f.get_eventdetection_read_name() << endl;
                    auto ed_params = f.get_eventdetection_event_parameters();
                    cout << "eventdetection_event_parameters=[abasic_found=" << ed_params.abasic_found
                        /*
                          << "abasic_event_index: " << ed_params.abasic_event_index << endl
                          << "abasic_peak_height: " << ed_params.abasic_peak_height << endl
                          << "hairpin_found: " << ed_params.hairpin_found << endl
                          << "hairpin_event_index: " << ed_params.hairpin_event_index << endl
                          << "hairpin_peak_height: " << ed_params.hairpin_peak_height << endl
                          << "hairpin_polyt_level: " << ed_params.hairpin_polyt_level << endl
                        */
                         << ", duration=" << ed_params.duration
                         << ", median_before=" << ed_params.median_before
                         << ", read_id=" << ed_params.read_id
                         << ", read_number=" << ed_params.read_number
                         << ", scaling_used=" << ed_params.scaling_used
                         << ", start_mux=" << ed_params.start_mux
                         << ", start_time=" << ed_params.start_time << "]" << endl;
                    auto ed_events = f.get_eventdetection_events();
                    cout << "eventdetection_events_size=" << ed_events.size() << endl;
                    for (const auto& e : ed_events)
                    {
                        cout << "  (mean=" << e.mean
                             << ", stdv=" << e.stdv
                             << ", start=" << e.start
                             << ", length=" << e.length << ")" << endl;
                    }
                } // have_eventdetection_events
            } // have_eventdetection_group

            //
            // inspect basecall group
            //
            bool have_basecall_version = f.have_basecall_version();
            cout << "have_basecall_version=" << have_basecall_version << endl;
            if (have_basecall_version)
            {
                cout << "basecall_version=" << f.basecall_version() << endl;
            }
            bool have_basecalled_group = f.have_basecalled_group();
            cout << "have_basecalled_group=" << have_basecalled_group << endl;
            bool have_basecalled_2D = f.have_basecalled_2D();
            cout << "have_basecalled_2D=" << have_basecalled_2D << endl;
            if (have_basecalled_2D)
            {
                // Extract the alignment between template and complement events
                // which were generated by the 2D basecaller
                auto v = f.get_event_alignments();
                cout << "event_alignment_size=" << v.size() << endl;
                for (const auto& e : v)
                {
                    cout << "  (template=" << e.template_index
                         << ", complement=" << e.complement_index
                         << ", kmer=" << e.kmer << ")" << endl;
                }
            }
            // Iterate over the template/complement strands
            for (size_t i = 0; i < 2; ++i)
            {
                cout << "have_basecalled_1D(" << i << ")=" << f.have_basecalled_1D(i) << endl;
                // Check if a pore model for this strand exists
                bool have_model = f.have_model(i);
                cout << "have_model(" << i << ")=" << have_model << endl;
                if (have_model)
                {
                    // Print the name of ONT's reference model used to basecall
                    cout << "model_file=" << f.get_model_file(i) << endl;
                    // Extract the global scaling parameters for the pore model
                    auto params = f.get_model_parameters(i);
                    cout << "model_parameters(" << i << ")=[drift=" << params.drift
                         << ", scale="     << params.scale
                         << ", scale_sd="  << params.scale_sd
                         << ", shift="     << params.shift
                         << ", var="       << params.var
                         << ", var_sd="    << params.var_sd << "]" << endl;
                    // Extract the expected current levels for each k-mer
                    auto v = f.get_model(i);
                    cout << "model_size(" << i << ")=" << v.size() << endl;
                    for (const auto& e : v)
                    {
                        cout << "  (kmer=" << e.kmer
                             << ", level_mean=" << e.level_mean
                             << ", level_stdv=" << e.level_stdv << ")" << endl;
                    }
                }
                // Check if this strand has event observations
                bool have_events = f.have_events(i);
                cout << "have_events(" << i << ")=" << have_events << endl;
                if (have_events)
                {
                    // Extract each event
                    auto v = f.get_events(i);
                    cout << "events_size(" << i << ")=" << v.size() << endl;
                    for (const auto& e : v)
                    {
                        cout << "  (mean=" << e.mean
                             << ", start=" << e.start
                             << ", stdv=" << e.stdv
                             << ", length=" << e.length << ")" << endl;
                    }
                }
            }
        }
        catch (hdf5_tools::Exception& e)
        {
            cout << "hdf5 error: " << e.what() << endl;
        }

        //
        // fast5 file is closed by its destructor at the end of this scope
        //
    }
    assert(fast5::File::get_object_count() == 0);
}
