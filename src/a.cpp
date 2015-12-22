#include <cassert>
#include <iostream>
#include <string>

#include "fast5.hpp"

using namespace std;


int main(int argc, char* argv[])
{
    assert(argc == 2);
    string file_name(argv[1]);
    //string ds_name(argv[2]);

    // Open the FAST5 file for reading
    if (not hdf5_tools::File_Reader::is_valid_file(file_name))
    {
        cerr << "error: not a hdf file [" << file_name << "]" << endl;
        return EXIT_FAILURE;
    }
    if (not fast5::File::is_valid_file(file_name))
    {
        cerr << "error: not a fast5 file [" << file_name << "]" << endl;
        return EXIT_FAILURE;
    }
    fast5::File* f_p;
    f_p = new fast5::File(file_name);

    // Check that it opened successfully
    assert(f_p->is_open());

    // Extract version information for the ONT software used to generate this dataset
    cout << "file_version=" << f_p->file_version() << endl;

    cout << "sampling rate: " << f_p->get_sampling_rate() << endl;

    if (f_p->have_sequences_group())
    {
        cout << "sequences_version=" << f_p->sequences_version() << endl;
    }

    if (f_p->have_basecalled_group())
    {
        // This function checks to see if 2D basecalls are available
        if(f_p->have_basecalled_2D())
        {
            if (f_p->have_basecall_version())
            {
                cout << "basecall_version=" << f_p->basecall_version() << endl;
            }

            cout << "basecalled_2D=" << f_p->basecalled_2D() << endl;

            // Extract the alignment between template and complement events
            // which were generated by the 2D basecaller
            auto v = f_p->get_event_alignments();
            cout << "event_alignment().size()=" << v.size() << endl;
            for (const auto& e : v)
            {
                cout << "(template=" << e.template_index << ", complement=" << e.complement_index << ", kmer=" << e.kmer << ")" << endl;
            }
        }

        // Iterate over the template/complement strands
        for (size_t i = 0; i < 2; ++i)
        {
            // Check if a pore model for this strand exists
            if (f_p->have_model(i))
            {
                // Print the name of ONT's reference model used to basecall
                cout << "Model file: " << f_p->get_model_file(i) << endl;

                // Extract the global scaling parameters for the pore model
                auto params = f_p->get_model_parameters(i);
                cout << "model drift=" << params.drift <<
                    ", scale="     << params.scale <<
                    ", scale_sd="  << params.scale_sd <<
                    ", shift="     << params.shift <<
                    ", var="       << params.var <<
                    ", var_sd="    << params.var_sd << endl;
            
                // Extract the expected current levels for each k-mer
                auto v = f_p->get_model(i);
                cout << "model(" << i << ").size()=" << v.size() << endl;
                for (const auto& e : v)
                {
                    cout << "(kmer=" << e.kmer << ", level_mean=" << e.level_mean << ", level_stdv=" << e.level_stdv << ")" << endl;
                }
            }

            // Check if this strand has event observations
            if (f_p->have_events(i))
            {
                // Extract each event
                auto v = f_p->get_events(i);
                cout << "events(" << i << ").size()=" << v.size() << endl;
                for (const auto& e : v)
                {
                    cout << "(mean=" << e.mean << ", start=" << e.start << ", stdv=" << e.stdv << ", length=" << e.length << ")" << endl;
                }
            }
        }
    } // if is_basecalled

    if (f_p->have_eventdetection_group())
    {
        cout << "eventdetection_version=" << f_p->eventdetection_version() << endl;
    }
    if (f_p->have_eventdetection_events())
    {
        cout << "EventDetection read name: " << f_p->get_eventdetection_read_name() << endl;
        auto ed_params = f_p->get_eventdetection_event_parameters();
        cout << "EventDetection event parameters:" << endl
             << "abasic_found: " << ed_params.abasic_found << endl
            /*
             << "abasic_event_index: " << ed_params.abasic_event_index << endl
             << "abasic_peak_height: " << ed_params.abasic_peak_height << endl
             << "hairpin_found: " << ed_params.hairpin_found << endl
             << "hairpin_event_index: " << ed_params.hairpin_event_index << endl
             << "hairpin_peak_height: " << ed_params.hairpin_peak_height << endl
             << "hairpin_polyt_level: " << ed_params.hairpin_polyt_level << endl
            */
             << "duration: " << ed_params.duration << endl
             << "median_before: " << ed_params.median_before << endl
             << "read_id: " << ed_params.read_id << endl
             << "read_number: " << ed_params.read_number << endl
             << "scaling_used: " << ed_params.scaling_used << endl
             << "start_mux: " << ed_params.start_mux << endl
             << "start_time: " << ed_params.start_time << endl;
        auto ed_events = f_p->get_eventdetection_events();
        cout << "EventDetection num events: " << ed_events.size() << endl;
        for (const auto& e : ed_events)
        {
            cout << "(mean=" << e.mean << ", stdv=" << e.stdv << ", start=" << e.start << ", length=" << e.length << ")" << endl;
        }
    }

    // Cleanup the file pointer, which closes the file
    delete f_p;
}
