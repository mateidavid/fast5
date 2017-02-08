#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>

#include <tclap/CmdLine.h>
#include "alg.hpp"
#include "logger.hpp"

#include "fast5.hpp"
#include "File_Packer.hpp"

using namespace std;

namespace opts
{
    using namespace TCLAP;
    string description = "Pack an ONT fast5 file.";
    CmdLine cmd_parser(description);
    //
    MultiArg< string > log_level("", "log", "Log level. (default: info)", false, "string", cmd_parser);
    MultiSwitchArg extra_verbosity("v", "", "Increase verbosity", cmd_parser);
    //
    ValueArg< unsigned > p_model_state_bits("", "p-model-state-bits", "P_Model_State bits to keep.", false, fast5::File_Packer::opts::p_model_state_bits(), "int", cmd_parser);
    ValueArg< unsigned > qv_bits("", "qv-bits", "QV bits to keep.", false, fast5::File_Packer::opts::max_qv_bits(), "int", cmd_parser);
    //
    SwitchArg al_drop("", "al-drop", "Drop basecall alignment data.", cmd_parser);
    SwitchArg al_copy("", "al-copy", "Copy basecall alignment data.", cmd_parser);
    SwitchArg al_unpack("", "al-unpack", "Unpack basecall alignment data.", cmd_parser);
    SwitchArg al_pack("", "al-pack", "Pack basecall alignment data.", cmd_parser);
    //
    SwitchArg ev_drop("", "ev-drop", "Drop basecall event data.", cmd_parser);
    SwitchArg ev_copy("", "ev-copy", "Copy basecall event data.", cmd_parser);
    SwitchArg ev_unpack("", "ev-unpack", "Unpack basecall event data.", cmd_parser);
    SwitchArg ev_pack("", "ev-pack", "Pack basecall event data.", cmd_parser);
    //
    SwitchArg fq_drop("", "fq-drop", "Drop basecall fastq data.", cmd_parser);
    SwitchArg fq_copy("", "fq-copy", "Copy basecall fastq data.", cmd_parser);
    SwitchArg fq_unpack("", "fq-unpack", "Unpack basecall fatsq data.", cmd_parser);
    SwitchArg fq_pack("", "fq-pack", "Pack basecall fastq data.", cmd_parser);
    //
    SwitchArg ed_drop("", "ed-drop", "Drop event detection data.", cmd_parser);
    SwitchArg ed_copy("", "ed-copy", "Copy event detection data.", cmd_parser);
    SwitchArg ed_unpack("", "ed-unpack", "Unpack event detection data.", cmd_parser);
    SwitchArg ed_pack("", "ed-pack", "Pack event detection data.", cmd_parser);
    //
    SwitchArg rw_drop("", "rw-drop", "Drop raw samples data.", cmd_parser);
    SwitchArg rw_copy("", "rw-copy", "Copy raw samples data.", cmd_parser);
    SwitchArg rw_unpack("", "rw-unpack", "Unpack raw samples data.", cmd_parser);
    SwitchArg rw_pack("", "rw-pack", "Pack raw samples data.", cmd_parser);
    //
    SwitchArg check("c", "check", "Check packing.", cmd_parser);
    SwitchArg unpack("u", "unpack", "Unpack fast5 file.", cmd_parser);
    SwitchArg pack("p", "pack", "Pack fast5 file (default, if no other pack/unpack/copy options).", cmd_parser);
    SwitchArg force("f", "force", "Overwrite output file if it exists.", cmd_parser);
    //
    UnlabeledValueArg< string > input_fn("input", "Input fast5 file.", true, "", "file", cmd_parser);
    UnlabeledValueArg< string > output_fn("output", "Output fast5 file.", true, "", "file", cmd_parser);
} // opts

void real_main()
{
    fast5::File src_f;
    fast5::File dst_f;
    try
    {
        // open files
        src_f.open(opts::input_fn);
        dst_f.create(opts::output_fn, opts::force);
        assert(src_f.is_open());
        assert(dst_f.is_open());
        assert(dst_f.is_rw());
        // copy all attributes
        fast5::File::copy_attributes(src_f, dst_f, "/UniqueGlobalKey", true);
        set< string > bc_gr_s;
        // process raw samples
        if (opts::rw_pack)
        {
            fast5::File_Packer::pack_rw(src_f, dst_f);
        }
        else if (opts::rw_unpack)
        {
            fast5::File_Packer::unpack_rw(src_f, dst_f);
        }
        else if (opts::rw_copy)
        {
            fast5::File_Packer::copy_rw(src_f, dst_f);
        }
        // process eventdetection events
        if (opts::ed_pack)
        {
            fast5::File_Packer::pack_ed(src_f, dst_f);
        }
        else if (opts::ed_unpack)
        {
            fast5::File_Packer::unpack_ed(src_f, dst_f);
        }
        else if (opts::ed_copy)
        {
            fast5::File_Packer::copy_ed(src_f, dst_f);
        }
        // process basecall fastq
        if (opts::fq_pack)
        {
            fast5::File_Packer::pack_fq(src_f, dst_f, bc_gr_s);
        }
        else if (opts::fq_unpack)
        {
            fast5::File_Packer::unpack_fq(src_f, dst_f, bc_gr_s);
        }
        else if (opts::fq_copy)
        {
            fast5::File_Packer::copy_fq(src_f, dst_f, bc_gr_s);
        }
        // process basecall events
        if (opts::ev_pack)
        {
            fast5::File_Packer::pack_ev(src_f, dst_f, bc_gr_s);
        }
        else if (opts::ev_unpack)
        {
            fast5::File_Packer::unpack_ev(src_f, dst_f, bc_gr_s);
        }
        else if (opts::ev_copy)
        {
            fast5::File_Packer::copy_ev(src_f, dst_f, bc_gr_s);
        }
        // process basecall alignments
        if (opts::al_pack)
        {
            fast5::File_Packer::pack_al(src_f, dst_f, bc_gr_s);
        }
        else if (opts::al_unpack)
        {
            fast5::File_Packer::unpack_al(src_f, dst_f, bc_gr_s);
        }
        else if (opts::al_copy)
        {
            fast5::File_Packer::copy_al(src_f, dst_f, bc_gr_s);
        }
        // copy basecall params
        fast5::File_Packer::copy_basecall_params(src_f, dst_f, bc_gr_s);
        // close files
        src_f.close();
        dst_f.close();
    }
    catch (hdf5_tools::Exception& e)
    {
        cerr << opts::input_fn.get() << ": HDF5 error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
} // real_main()

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    // set log levels
    auto default_level = (int)logger::level::info + opts::extra_verbosity.getValue();
    logger::Logger::set_default_level(default_level);
    logger::Logger::set_levels_from_options(opts::log_level, &clog);
    // print options
    LOG(info) << "program: " << opts::cmd_parser.getProgramName() << endl;
    LOG(info) << "version: " << opts::cmd_parser.getVersion() << endl;
    LOG(info) << "args: " << opts::cmd_parser.getOrigArgv() << endl;
    // what to pack/unpack
    if (opts::pack + opts::unpack > 1)
    {
        LOG(error) << "at most one of --pack/--unpack may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop > 1)
    {
        LOG(error) << "at most one of --rw-pack/--rw-unpack/--rw-copy/--rw-drop may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop > 1)
    {
        LOG(error) << "at most one of --ed-pack/--ed-unpack/--ed-copy/--ed-drop may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::fq_pack + opts::fq_unpack + opts::fq_copy + opts::fq_drop > 1)
    {
        LOG(error) << "at most one of --fq-pack/--fq-unpack/--fq-copy/--fq-drop may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop > 1)
    {
        LOG(error) << "at most one of --ev-pack/--ev-unpack/--ev-copy/--ev-drop may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::al_pack + opts::al_unpack + opts::al_copy + opts::al_drop > 1)
    {
        LOG(error) << "at most one of --al-pack/--al-unpack/--al-copy/--al-drop may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::pack + opts::unpack
        + opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop
        + opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop
        + opts::fq_pack + opts::fq_unpack + opts::fq_copy + opts::fq_drop
        + opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop
        + opts::al_pack + opts::al_unpack + opts::al_copy + opts::al_drop
        == 0)
    {
        opts::pack.set(true);
    }
    if (opts::pack)
    {
        opts::rw_pack.set(true);
        opts::ed_pack.set(true);
        opts::fq_pack.set(true);
        opts::ev_pack.set(true);
        opts::al_pack.set(true);
    }
    else if (opts::unpack)
    {
        opts::rw_unpack.set(true);
        opts::ed_unpack.set(true);
        opts::fq_unpack.set(true);
        opts::ev_unpack.set(true);
        opts::al_unpack.set(true);
    }
    else
    {
        if (opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop == 0) opts::rw_copy.set(true);
        if (opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop == 0) opts::ed_copy.set(true);
        if (opts::fq_pack + opts::fq_unpack + opts::fq_copy + opts::fq_drop == 0) opts::fq_copy.set(true);
        if (opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop == 0) opts::ev_copy.set(true);
        if (opts::al_pack + opts::al_unpack + opts::al_copy + opts::al_drop == 0) opts::al_copy.set(true);
    }
    LOG(info) << "rw: " << (opts::rw_pack? "pack" : opts::rw_unpack? "unpack" : opts::rw_copy? "copy" : "drop") << endl;
    LOG(info) << "ed: " << (opts::ed_pack? "pack" : opts::ed_unpack? "unpack" : opts::ed_copy? "copy" : "drop") << endl;
    LOG(info) << "fq: " << (opts::fq_pack? "pack" : opts::fq_unpack? "unpack" : opts::fq_copy? "copy" : "drop") << endl;
    LOG(info) << "ev: " << (opts::ev_pack? "pack" : opts::ev_unpack? "unpack" : opts::ev_copy? "copy" : "drop") << endl;
    LOG(info) << "al: " << (opts::al_pack? "pack" : opts::al_unpack? "unpack" : opts::al_copy? "copy" : "drop") << endl;
    LOG(info) << "check: " << (opts::check? "yes" : "no") << endl;
    // set File_Packer options
    fast5::File_Packer::opts::check() = opts::check;
    fast5::File_Packer::opts::qv_bits() = opts::qv_bits;
    fast5::File_Packer::opts::p_model_state_bits() = opts::p_model_state_bits;
    real_main();
}
