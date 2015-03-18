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
    //DataSet* ds_p = new DataSet (f_p->openDataSet(ds_name));

    //delete ds_p;
    delete f_p;
}
