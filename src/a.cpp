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

    fast5::File* f_p;
    f_p = new fast5::File(file_name);
    assert(f_p->is_open());
    cout << "file_version=" << f_p->file_version() << endl;
    cout << "basecall_version=" << f_p->basecall_version() << endl;
    cout << "eventdetection_version=" << f_p->eventdetection_version() << endl;
    cout << "sequences_version=" << f_p->sequences_version() << endl;

    for (size_t i = 0; i < 2; ++i)
    {
        if (f_p->have_model(i))
        {
            auto v = f_p->get_model(i);
            cout << "model(" << i << ").size()=" << v.size() << endl;
            for (const auto& e : v)
            {
                cout << "(kmer=" << e.kmer << ", level_mean=" << e.level_mean << ", level_stdv=" << e.level_stdv << ")" << endl;
            }
        }
        if (f_p->have_events(i))
        {
            auto v = f_p->get_events(i);
            cout << "events(" << i << ").size()=" << v.size() << endl;
            for (const auto& e : v)
            {
                cout << "(mean=" << e.mean << ", start=" << e.start << ", stdv=" << e.stdv << ", length=" << e.length << ")" << endl;
            }
        }
    }

    delete f_p;
}
