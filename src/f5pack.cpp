#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>

#include <tclap/CmdLine.h>
#include "alg.hpp"
#include "logger.hpp"

#include "fast5.hpp"

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
    //ValueArg< unsigned > float_prec("", "float-prec", "Float precision.", false, 10, "int", cmd_parser);
    //SwitchArg rw_time("", "rw-time", "Add timepoints to raw data.", cmd_parser);
    //SwitchArg curr_int("", "curr-int", "Dump current data encoded as int (raw samples only).", cmd_parser);
    //SwitchArg time_int("", "time-int", "Dump start/length data encoded as int.", cmd_parser);
    //
    //ValueArg< string > rn("", "rn", "Read name.", false, "", "Read_1015|...", cmd_parser);
    //ValueArg< unsigned > st("", "st", "Strand.", false, 0, "0|1|2", cmd_parser);
    //ValueArg< string > gr("", "gr", "Group name suffix.", false, "", "000|RNN_001|...", cmd_parser);
    //
    //SwitchArg fq("", "fq", "Dump basecall fastq data.", cmd_parser);
    SwitchArg ev_drop("", "ev-drop", "Drop basecall event data.", cmd_parser);
    SwitchArg ev_copy("", "ev-copy", "Copy basecall event data.", cmd_parser);
    SwitchArg ev_unpack("", "ev-unpack", "Unpack basecall event data.", cmd_parser);
    SwitchArg ev_pack("", "ev-pack", "Pack basecall event data.", cmd_parser);
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

template < typename U, typename V >
void print_map(ostream& os, const map< U, V >& m, const string& prefix)
{
    for (const auto& p : m)
    {
        os << prefix << p.first << "=" << p.second << endl;
    }
}

unsigned time_int(double f, fast5::Channel_Id_Parameters const & channel_id_params)
{
    return f * channel_id_params.sampling_rate;
}

void do_pack_rw(fast5::File const & src_f, fast5::File const & dst_f)
{
    auto rn_l = src_f.get_raw_samples_read_name_list();
    for (auto const & rn : rn_l)
    {
        auto rsi = src_f.get_raw_samples_int(rn);
        auto rs_pack = src_f.pack_rw(rsi);
        if (opts::check)
        {
            auto rs_unpack = src_f.unpack_rw(rs_pack);
            if (rs_unpack.size() != rsi.size())
            {
                LOG(error)
                    << "check_failed rs_unpack.size=" << rs_unpack.size()
                    << " rs_orig.size=" << rsi.size() << endl;
                exit(EXIT_FAILURE);
            }
            for (unsigned i = 0; i < rs_unpack.size(); ++i)
            {
                if (rs_unpack[i] != rsi[i])
                {
                    LOG(error)
                        << "check_failed i=" << i
                        << " rs_unpack=" << rs_unpack[i]
                        << " rs_orig=" << rsi[i] << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        dst_f.add_raw_samples_pack(rn, rs_pack);
    }
}
void do_unpack_rw(fast5::File const & src_f, fast5::File const & dst_f)
{
    auto rn_l = src_f.get_raw_samples_read_name_list();
    for (auto const & rn : rn_l)
    {
        auto rsi = src_f.get_raw_samples_int(rn);
        dst_f.add_raw_samples_int(rn, rsi);
    }
}
void do_copy_rw(fast5::File const & src_f, fast5::File const & dst_f)
{
    auto rn_l = src_f.get_raw_samples_read_name_list();
    for (auto const & rn : rn_l)
    {
        if (src_f.have_raw_samples_unpack(rn))
        {
            auto rsi = src_f.get_raw_samples_int(rn);
            dst_f.add_raw_samples_int(rn, rsi);
        }
        else if (src_f.have_raw_samples_pack(rn))
        {
            auto p = src_f.get_raw_samples_pack(rn);
            dst_f.add_raw_samples_pack(rn, p);
        }
    }
}

void do_pack_ed(fast5::File const & src_f, fast5::File const & dst_f)
{
    auto gr_l = src_f.get_eventdetection_group_list();
    for (auto const & gr : gr_l)
    {
        auto rn_l = src_f.get_eventdetection_read_name_list(gr);
        for (auto const & rn : rn_l)
        {
            auto ed = src_f.get_eventdetection_events(gr, rn);
            auto ed_param = src_f.get_eventdetection_event_params(gr, rn);
            auto ed_pack = src_f.pack_ed(ed, ed_param);
            if (opts::check)
            {
                auto rs = src_f.get_raw_samples(rn);
                auto rs_param = src_f.get_raw_samples_params(rn);
                auto ed_unpack = src_f.unpack_ed(ed_pack, ed_param, rs, rs_param);
                if (ed_unpack.size() != ed.size())
                {
                    LOG(error)
                        << "check_failed gr=" << gr
                        << " ed_unpack.size=" << ed_unpack.size()
                        << " ed_orig.size=" << ed.size() << endl;
                    exit(EXIT_FAILURE);
                }
                for (unsigned i = 0; i + 1 < ed_unpack.size(); ++i) // skip last event
                {
                    LOG(debug1)
                        << "gr=" << gr
                        << " i=" << i
                        << " ed_unpack=(" << ed_unpack[i].start
                        << "," << ed_unpack[i].length
                        << "," << ed_unpack[i].mean
                        << "," << ed_unpack[i].stdv
                        << ") ed_orig=(" << ed[i].start
                        << "," << ed[i].length
                        << "," << ed[i].mean
                        << "," << ed[i].stdv
                        << ")" << endl;
                    if (ed_unpack[i].start != ed[i].start
                        or ed_unpack[i].length != ed[i].length
                        or abs(ed_unpack[i].mean - ed[i].mean) > .1
                        or abs(ed_unpack[i].stdv - ed[i].stdv) > .1)
                    {
                        LOG(error)
                            << "check_failed gr=" << gr
                            << " i=" << i
                            << " ed_unpack=(" << ed_unpack[i].start
                            << "," << ed_unpack[i].length
                            << "," << ed_unpack[i].mean
                            << "," << ed_unpack[i].stdv
                            << ") ed_orig=(" << ed[i].start
                            << "," << ed[i].length
                            << "," << ed[i].mean
                            << "," << ed[i].stdv
                            << ")" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            } // if check
            dst_f.add_eventdetection_events_pack(gr, rn, ed_pack);
        } // for rn
    } // for gr
} // do_pack_ed()
void do_unpack_ed(fast5::File const & src_f, fast5::File const & dst_f)
{
    auto gr_l = src_f.get_eventdetection_group_list();
    for (auto const & gr : gr_l)
    {
        auto rn_l = src_f.get_eventdetection_read_name_list(gr);
        for (auto const & rn : rn_l)
        {
            auto ed = src_f.get_eventdetection_events(gr, rn);
            dst_f.add_eventdetection_events(gr, rn, ed);
        }
    }
}
void do_copy_ed(fast5::File const & src_f, fast5::File const & dst_f)
{
    auto gr_l = src_f.get_eventdetection_group_list();
    for (auto const & gr : gr_l)
    {
        auto rn_l = src_f.get_eventdetection_read_name_list(gr);
        for (auto const & rn : rn_l)
        {
            if (src_f.have_eventdetection_events_unpack(gr, rn))
            {
                auto ed = src_f.get_eventdetection_events(gr, rn);
                dst_f.add_eventdetection_events(gr, rn, ed);
            }
            else if (src_f.have_eventdetection_events_pack(gr, rn))
            {
                auto ed_pack = src_f.get_eventdetection_events_pack(gr, rn);
                dst_f.add_eventdetection_events_pack(gr, rn, ed_pack);
            }
        }
    }
}
void do_pack_ev(fast5::File const & src_f, fast5::File const & dst_f)
{
    for (unsigned st = 0; st < 2; ++st)
    {
        auto gr_l = src_f.get_basecall_strand_group_list(st);
        for (auto const & gr : gr_l)
        {
            if (src_f.have_basecall_events_pack(st, gr))
            {
                auto ev_pack = src_f.get_basecall_events_pack(st, gr);
                dst_f.add_basecall_events_pack(st, gr, ev_pack);
            }
            else if (src_f.have_basecall_events_unpack(st, gr))
            {
                auto ev = src_f.get_basecall_events(st, gr);
                auto ev_param = src_f.get_basecall_event_params(st, gr);
                // sampling rate
                auto channel_id_params = src_f.get_channel_id_params();
                // ed group
                auto ed_gr = src_f.get_basecall_eventdetection_group(gr);
                if (ed_gr.empty())
                {
                    LOG(error)
                        << "missing eventdetection events for basecall events: st=" << st << " gr=" << gr << endl;
                    exit(EXIT_FAILURE);
                }
                auto ed = src_f.get_eventdetection_events(ed_gr);
                auto ev_pack = src_f.pack_ev(ev, ev_param, ed, ed_gr, channel_id_params.sampling_rate);
                if (opts::check)
                {
                    auto ev_unpack = src_f.unpack_ev(ev_pack, ed, channel_id_params.sampling_rate);
                    if (ev_unpack.size() != ev.size())
                    {
                        LOG(error)
                            << "check_failed st=" << st
                            << " gr=" << gr
                            << " ev_unpack.size=" << ev_unpack.size()
                            << " ev_orig.size=" << ev.size() << endl;
                        exit(EXIT_FAILURE);
                    }
                    for (unsigned i = 0; i < ev_unpack.size(); ++i)
                    {
                        if (abs(ev_unpack[i].start - ev[i].start) > 1e-3
                            or abs(ev_unpack[i].length - ev[i].length) > 1e-3
                            or abs(ev_unpack[i].mean - ev[i].mean) > 1e-1
                            or abs(ev_unpack[i].stdv - ev[i].stdv) > 1e-1
                            or ev_unpack[i].move != ev[i].move
                            or ev_unpack[i].model_state != ev[i].model_state)
                        {
                            LOG(error)
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " i=" << i
                                << " ev_unpack=(" << ev_unpack[i].start
                                << "," << ev_unpack[i].length
                                << "," << ev_unpack[i].mean
                                << "," << ev_unpack[i].stdv
                                << "," << ev_unpack[i].move
                                << "," << ev_unpack[i].get_model_state()
                                << ") ev_orig=(" << ev[i].start
                                << "," << ev[i].length
                                << "," << ev[i].mean
                                << "," << ev[i].stdv
                                << "," << ev[i].move
                                << "," << ev[i].get_model_state()
                                << ")" << endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                }
                dst_f.add_basecall_events_pack(st, gr, ev_pack);
            }
        }
    }
} // do_pack_ev
void do_unpack_ev(fast5::File const & src_f, fast5::File const & dst_f)
{
    for (unsigned st = 0; st < 2; ++st)
    {
        auto gr_l = src_f.get_basecall_strand_group_list(st);
        for (auto const & gr : gr_l)
        {
            if (not src_f.have_basecall_events(st, gr)) continue;
            auto ev = src_f.get_basecall_events(st, gr);
            auto ev_param = src_f.get_basecall_event_params(st, gr);
            dst_f.add_basecall_events(st, gr, ev);
            dst_f.add_basecall_event_params(st, gr, ev_param);
        }
    }
} // do_unpack_ev
void do_copy_ev(fast5::File const & src_f, fast5::File const & dst_f)
{
    for (unsigned st = 0; st < 2; ++st)
    {
        auto gr_l = src_f.get_basecall_strand_group_list(st);
        for (auto const & gr : gr_l)
        {
            if (src_f.have_basecall_events_unpack(st, gr))
            {
                auto ev = src_f.get_basecall_events(st, gr);
                auto ev_param = src_f.get_basecall_event_params(st, gr);
                dst_f.add_basecall_events(st, gr, ev);
                dst_f.add_basecall_event_params(st, gr, ev_param);
            }
            else if (src_f.have_basecall_events_pack(st, gr))
            {
                auto ev_pack = src_f.get_basecall_events_pack(st, gr);
                dst_f.add_basecall_events_pack(st, gr, ev_pack);
            }
        }
    }
} // do_copy_ev

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
        fast5::File::copy_attributes(src_f, dst_f);
        // process raw samples
        if (opts::rw_pack)
        {
            do_pack_rw(src_f, dst_f);
        }
        else if (opts::rw_unpack)
        {
            do_unpack_rw(src_f, dst_f);
        }
        else if (opts::rw_copy)
        {
            do_copy_rw(src_f, dst_f);
        }
        // process eventdetection events
        if (opts::ed_pack)
        {
            do_pack_ed(src_f, dst_f);
        }
        else if (opts::ed_unpack)
        {
            do_unpack_ed(src_f, dst_f);
        }
        else if (opts::ed_copy)
        {
            do_copy_ed(src_f, dst_f);
        }
        // process basecall events
        if (opts::ev_pack)
        {
            do_pack_ev(src_f, dst_f);
        }
        else if (opts::ev_unpack)
        {
            do_unpack_ev(src_f, dst_f);
        }
        else if (opts::ev_copy)
        {
            do_copy_ev(src_f, dst_f);
        }
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
    if (opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop > 1)
    {
        LOG(error) << "at most one of --ev-pack/--ev-unpack/--ev-copy/--ev-drop may be given" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts::pack + opts::unpack
        + opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop
        + opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop
        + opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop == 0)
    {
        opts::pack.set(true);
    }
    if (opts::pack)
    {
        opts::rw_pack.set(true);
        opts::ed_pack.set(true);
        opts::ev_pack.set(true);
    }
    else if (opts::unpack)
    {
        opts::rw_unpack.set(true);
        opts::ed_unpack.set(true);
        opts::ev_unpack.set(true);
    }
    else
    {
        if (opts::rw_pack + opts::rw_unpack + opts::rw_copy + opts::rw_drop == 0) opts::rw_copy.set(true);
        if (opts::ed_pack + opts::ed_unpack + opts::ed_copy + opts::ed_drop == 0) opts::ed_copy.set(true);
        if (opts::ev_pack + opts::ev_unpack + opts::ev_copy + opts::ev_drop == 0) opts::ev_copy.set(true);
    }
    LOG(info) << "rw: " << (opts::rw_pack? "pack" : opts::rw_unpack? "unpack" : opts::rw_copy? "copy" : "drop") << endl;
    LOG(info) << "ed: " << (opts::ed_pack? "pack" : opts::ed_unpack? "unpack" : opts::ed_copy? "copy" : "drop") << endl;
    LOG(info) << "ev: " << (opts::ev_pack? "pack" : opts::ev_unpack? "unpack" : opts::ev_copy? "copy" : "drop") << endl;
    LOG(info) << "check: " << (opts::check? "yes" : "no") << endl;
    real_main();
}
