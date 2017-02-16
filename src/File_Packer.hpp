#ifndef __FILE_PACKER_HPP
#define __FILE_PACKER_HPP

#include <string>
#include <set>

#include "fast5.hpp"
#include "logger.hpp"

#define STATIC_MEMBER_WRAPPER(_type, _id, _init) \
    static _type & _id() { static _type _ ## _id = _init; return _ ## _id; }

namespace fast5
{

class File_Packer
{
public:
    File_Packer() :
        File_Packer(1)
    {}

    File_Packer(int _policy) :
        File_Packer(_policy, _policy, _policy, _policy, _policy)
    {}

    File_Packer(int _rw_policy, int _ed_policy, int _fq_policy, int _ev_policy, int _al_policy) :
        rw_policy(_rw_policy),
        ed_policy(_ed_policy),
        fq_policy(_fq_policy),
        ev_policy(_ev_policy),
        al_policy(_al_policy),
        check(true),
        force(false),
        qv_bits(max_qv_bits()),
        p_model_state_bits(default_p_model_state_bits())
    {}

    void set_check(bool _check) { check = _check; }
    void set_force(bool _force) { force = _force; }
    void set_qv_bits(unsigned _qv_bits) { qv_bits = _qv_bits; }
    void set_p_model_state_bits(unsigned _p_model_state_bits) { p_model_state_bits = _p_model_state_bits; }

    STATIC_MEMBER_WRAPPER(unsigned const, max_qv_bits, 5)
    STATIC_MEMBER_WRAPPER(unsigned const, max_qv_mask, ((unsigned)1 << max_qv_bits()) - 1)
    STATIC_MEMBER_WRAPPER(unsigned const, default_p_model_state_bits, 2)

    void
    run(std::string const & ifn, std::string const & ofn) const
    {
        File src_f;
        File dst_f;
        try
        {
            // open files
            src_f.open(ifn);
            dst_f.create(ofn, force);
            assert(src_f.is_open());
            assert(dst_f.is_open());
            assert(dst_f.is_rw());
            // copy attributes under / and /UniqueGlobalKey
            copy_attributes(src_f, dst_f, "", false);
            copy_attributes(src_f, dst_f, "/UniqueGlobalKey", true);
            std::set< std::string > bc_gr_s;
            // process raw samples
            if (rw_policy == 1)
            {
                pack_rw(src_f, dst_f);
            }
            else if (rw_policy == 2)
            {
                unpack_rw(src_f, dst_f);
            }
            else if (rw_policy == 3)
            {
                copy_rw(src_f, dst_f);
            }
            // process eventdetection events
            if (ed_policy == 1)
            {
                pack_ed(src_f, dst_f);
            }
            else if (ed_policy == 2)
            {
                unpack_ed(src_f, dst_f);
            }
            else if (ed_policy == 3)
            {
                copy_ed(src_f, dst_f);
            }
            // process basecall fastq
            if (fq_policy == 1)
            {
                pack_fq(src_f, dst_f, bc_gr_s);
            }
            else if (fq_policy == 2)
            {
                unpack_fq(src_f, dst_f, bc_gr_s);
            }
            else if (fq_policy == 3)
            {
                copy_fq(src_f, dst_f, bc_gr_s);
            }
            // process basecall events
            if (ev_policy == 1)
            {
                pack_ev(src_f, dst_f, bc_gr_s);
            }
            else if (ev_policy == 2)
            {
                unpack_ev(src_f, dst_f, bc_gr_s);
            }
            else if (ev_policy == 3)
            {
                copy_ev(src_f, dst_f, bc_gr_s);
            }
            // process basecall alignments
            if (al_policy == 1)
            {
                pack_al(src_f, dst_f, bc_gr_s);
            }
            else if (al_policy == 2)
            {
                unpack_al(src_f, dst_f, bc_gr_s);
            }
            else if (al_policy == 3)
            {
                copy_al(src_f, dst_f, bc_gr_s);
            }
            // copy basecall params
            copy_basecall_params(src_f, dst_f, bc_gr_s);
            // close files
            src_f.close();
            dst_f.close();
        }
        catch (hdf5_tools::Exception& e)
        {
            LOG_EXIT << ifn << ": HDF5 error: " << e.what() << std::endl;
        }
    } // run()

private:
    int rw_policy;
    int ed_policy;
    int fq_policy;
    int ev_policy;
    int al_policy;
    bool check;
    bool force;
    unsigned qv_bits;
    unsigned p_model_state_bits;

    void
    pack_rw(File const & src_f, File & dst_f) const
    {
        auto rn_l = src_f.get_raw_samples_read_name_list();
        for (auto const & rn : rn_l)
        {
            if (src_f.have_raw_samples_pack(rn))
            {
                auto rs_pack = src_f.get_raw_samples_pack(rn);
                dst_f.add_raw_samples(rn, rs_pack);
            }
            else if (src_f.have_raw_samples_unpack(rn))
            {
                auto rsi_ds = src_f.get_raw_int_samples_dataset(rn);
                auto & rsi = rsi_ds.first;
                auto & rs_params = rsi_ds.second;
                auto rs_pack = src_f.pack_rw(rsi_ds);
                if (check)
                {
                    auto rsi_ds_unpack = src_f.unpack_rw(rs_pack);
                    auto & rsi_unpack = rsi_ds_unpack.first;
                    auto & rs_params_unpack = rsi_ds_unpack.second;
                    if (not (rs_params_unpack == rs_params))
                    {
                        LOG_ABORT
                            << "check_failed rs_params_unpack!=rs_params" << std::endl;
                    }
                    if (rsi_unpack.size() != rsi.size())
                    {
                        LOG_ABORT
                            << "check_failed rs_unpack.size=" << rsi_unpack.size()
                            << " rs_orig.size=" << rsi.size() << std::endl;
                    }
                    for (unsigned i = 0; i < rsi_unpack.size(); ++i)
                    {
                        if (rsi_unpack[i] != rsi[i])
                        {
                            LOG_ABORT
                                << "check_failed i=" << i
                                << " rs_unpack=" << rsi_unpack[i]
                                << " rs_orig=" << rsi[i] << std::endl;
                        }
                    }
                }
                dst_f.add_raw_samples(rn, rs_pack);
                LOG(info)
                    << "rn=" << rn
                    << " rs_size=" << rsi.size()
                    << " signal_bits=" << rs_pack.signal_params.at("avg_bits")
                    << std::endl;
            }
        }
    } // pack_rw()

    void
    unpack_rw(File const & src_f, File & dst_f) const
    {
        auto rn_l = src_f.get_raw_samples_read_name_list();
        for (auto const & rn : rn_l)
        {
            auto rsi_ds = src_f.get_raw_int_samples_dataset(rn);
            dst_f.add_raw_samples_dataset(rn, rsi_ds);
        }
    } // unpack_rw()

    void
    copy_rw(File const & src_f, File & dst_f) const
    {
        auto rn_l = src_f.get_raw_samples_read_name_list();
        for (auto const & rn : rn_l)
        {
            if (src_f.have_raw_samples_unpack(rn))
            {
                auto rsi_ds = src_f.get_raw_int_samples_dataset(rn);
                dst_f.add_raw_samples_dataset(rn, rsi_ds);
            }
            else if (src_f.have_raw_samples_pack(rn))
            {
                auto rs_pack = src_f.get_raw_samples_pack(rn);
                dst_f.add_raw_samples(rn, rs_pack);
            }
        }
    } // copy_rw()

    void
    pack_ed(File const & src_f, File & dst_f) const
    {
        auto gr_l = src_f.get_eventdetection_group_list();
        for (auto const & gr : gr_l)
        {
            auto rn_l = src_f.get_eventdetection_read_name_list(gr);
            for (auto const & rn : rn_l)
            {
                auto ed_params = src_f.get_eventdetection_params(gr);
                dst_f.add_eventdetection_params(gr, ed_params);
                if (src_f.have_eventdetection_events_pack(gr, rn))
                {
                    auto ede_pack = src_f.get_eventdetection_events_pack(gr, rn);
                    dst_f.add_eventdetection_events(gr, rn, ede_pack);
                }
                else if (src_f.have_eventdetection_events(gr, rn))
                {
                    auto ede_ds = src_f.get_eventdetection_events_dataset(gr, rn);
                    auto & ede = ede_ds.first;
                    auto & ede_params = ede_ds.second;
                    auto ede_pack = src_f.pack_ed(ede_ds);
                    if (check)
                    {
                        auto rs_ds = src_f.get_raw_samples_dataset(rn);
                        auto ede_ds_unpack = src_f.unpack_ed(ede_pack, rs_ds);
                        auto & ede_unpack = ede_ds_unpack.first;
                        auto & ede_params_unpack = ede_ds_unpack.second;
                        if (not (ede_params_unpack == ede_params))
                        {
                            LOG_ABORT
                                << "check_failed ede_params_unpack!=ede_params" << std::endl;
                        }
                        if (ede_unpack.size() != ede.size())
                        {
                            LOG_ABORT
                                << "check_failed gr=" << gr
                                << " ede_unpack.size=" << ede_unpack.size()
                                << " ede_orig.size=" << ede.size() << std::endl;
                        }
                        for (unsigned i = 0; i + 1 < ede_unpack.size(); ++i) // skip last event
                        {
                            LOG(debug1)
                                << "gr=" << gr
                                << " i=" << i
                                << " ede_unpack=(" << ede_unpack[i].start
                                << "," << ede_unpack[i].length
                                << "," << ede_unpack[i].mean
                                << "," << ede_unpack[i].stdv
                                << ") ed_orig=(" << ede[i].start
                                << "," << ede[i].length
                                << "," << ede[i].mean
                                << "," << ede[i].stdv
                                << ")" << std::endl;
                            if (ede_unpack[i].start != ede[i].start
                                or ede_unpack[i].length != ede[i].length
                                or abs(ede_unpack[i].mean - ede[i].mean) > .1
                                or abs(ede_unpack[i].stdv - ede[i].stdv) > .1)
                            {
                                LOG_ABORT
                                    << "check_failed gr=" << gr
                                    << " i=" << i
                                    << " ede_unpack=(" << ede_unpack[i].start
                                    << "," << ede_unpack[i].length
                                    << "," << ede_unpack[i].mean
                                    << "," << ede_unpack[i].stdv
                                    << ") ed_orig=(" << ede[i].start
                                    << "," << ede[i].length
                                    << "," << ede[i].mean
                                    << "," << ede[i].stdv
                                    << ")" << std::endl;
                            }
                        }
                    } // if check
                    dst_f.add_eventdetection_events(gr, rn, ede_pack);
                    LOG(info)
                        << "gr=" << gr
                        << " rn=" << rn
                        << " ed_size=" << ede.size()
                        << " skip_bits=" << ede_pack.skip_params.at("avg_bits")
                        << " len_bits=" << ede_pack.len_params.at("avg_bits")
                        << std::endl;
                }
            } // for rn
        } // for gr
    } // pack_ed()

    void
    unpack_ed(File const & src_f, File & dst_f) const
    {
        auto gr_l = src_f.get_eventdetection_group_list();
        for (auto const & gr : gr_l)
        {
            auto rn_l = src_f.get_eventdetection_read_name_list(gr);
            for (auto const & rn : rn_l)
            {
                auto ed_params = src_f.get_eventdetection_params(gr);
                dst_f.add_eventdetection_params(gr, ed_params);
                auto ede_ds = src_f.get_eventdetection_events_dataset(gr, rn);
                dst_f.add_eventdetection_events_dataset(gr, rn, ede_ds);
            }
        }
    } // unpack_ed()

    void
    copy_ed(File const & src_f, File & dst_f) const
    {
        auto gr_l = src_f.get_eventdetection_group_list();
        for (auto const & gr : gr_l)
        {
            auto rn_l = src_f.get_eventdetection_read_name_list(gr);
            for (auto const & rn : rn_l)
            {
                auto ed_params = src_f.get_eventdetection_params(gr);
                dst_f.add_eventdetection_params(gr, ed_params);
                if (src_f.have_eventdetection_events_unpack(gr, rn))
                {
                    auto ede_ds = src_f.get_eventdetection_events_dataset(gr, rn);
                    dst_f.add_eventdetection_events_dataset(gr, rn, ede_ds);
                }
                else if (src_f.have_eventdetection_events_pack(gr, rn))
                {
                    auto ede_pack = src_f.get_eventdetection_events_pack(gr, rn);
                    dst_f.add_eventdetection_events(gr, rn, ede_pack);
                }
            }
        }
    } // copy_ed()

    void
    pack_fq(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        for (unsigned st = 0; st < 3; ++st)
        {
            auto gr_l = src_f.get_basecall_strand_group_list(st);
            for (auto const & gr : gr_l)
            {
                if (src_f.have_basecall_fastq_pack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto fq_pack = src_f.get_basecall_fastq_pack(st, gr);
                    dst_f.add_basecall_fastq(st, gr, fq_pack);
                }
                else if (src_f.have_basecall_fastq_unpack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto fq = src_f.get_basecall_fastq(st, gr);
                    auto fqa = src_f.split_fq(fq);
                    auto fq_pack = src_f.pack_fq(fq, qv_bits);
                    if (check)
                    {
                        auto fq_unpack = src_f.unpack_fq(fq_pack);
                        auto fqa_unpack = src_f.split_fq(fq_unpack);
                        if (fqa_unpack[0] != fqa[0])
                        {
                            LOG_ABORT
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " fq_unpack_name=" << fqa_unpack[0]
                                << " fq_orig_name=" << fqa[0] << std::endl;
                        }
                        if (fqa_unpack[1] != fqa[1])
                        {
                            LOG_ABORT
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " fq_unpack_bp=" << fqa_unpack[1]
                                << " fq_orig_bp=" << fqa[1] << std::endl;
                        }
                        if (fqa_unpack[3].size() != fqa[3].size())
                        {
                            LOG_ABORT
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " fq_unpack_qv_size=" << fqa_unpack[3].size()
                                << " fq_orig_qv_size=" << fqa[3].size() << std::endl;
                        }
                        auto qv_mask = max_qv_mask() & (max_qv_mask() << (max_qv_bits() - qv_bits));
                        for (unsigned i = 0; i < fqa_unpack[3].size(); ++i)
                        {
                            if ((std::min<unsigned>(fqa_unpack[3][i] - 33, max_qv_mask()) & qv_mask) !=
                                (std::min<unsigned>(fqa[3][i] - 33, max_qv_mask()) & qv_mask))
                            {
                                LOG_ABORT
                                    << "check_failed st=" << st
                                    << " gr=" << gr
                                    << " i=" << i
                                    << " fq_unpack_qv=" << fqa_unpack[3][i]
                                    << " fq_orig_qv=" << fqa[3][i] << std::endl;
                            }
                        }
                    }
                    dst_f.add_basecall_fastq(st, gr, fq_pack);
                    LOG(info)
                        << "gr=" << gr
                        << " st=" << st
                        << " bp_size=" << fqa[1].size()
                        << " bp_bits=" << fq_pack.bp_params.at("avg_bits")
                        << " qv_bits=" << fq_pack.qv_params.at("avg_bits")
                        << std::endl;
                }
            }
        }
    } // pack_fq()

    void
    unpack_fq(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        for (unsigned st = 0; st < 3; ++st)
        {
            auto gr_l = src_f.get_basecall_strand_group_list(st);
            for (auto const & gr : gr_l)
            {
                if (src_f.have_basecall_fastq(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto fq = src_f.get_basecall_fastq(st, gr);
                    dst_f.add_basecall_fastq(st, gr, fq);
                }
            }
        }
    } // unpack_fq()

    void
    copy_fq(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        for (unsigned st = 0; st < 3; ++st)
        {
            auto gr_l = src_f.get_basecall_strand_group_list(st);
            for (auto const & gr : gr_l)
            {
                if (src_f.have_basecall_fastq_unpack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto fq = src_f.get_basecall_fastq(st, gr);
                    dst_f.add_basecall_fastq(st, gr, fq);
                }
                else if (src_f.have_basecall_fastq_pack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto fq_pack = src_f.get_basecall_fastq_pack(st, gr);
                    dst_f.add_basecall_fastq(st, gr, fq_pack);
                }
            }
        }
    } // copy_fq()

    void
    pack_ev(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        for (unsigned st = 0; st < 2; ++st)
        {
            auto gr_l = src_f.get_basecall_strand_group_list(st);
            for (auto const & gr : gr_l)
            {
                if (src_f.have_basecall_events_pack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto ev_pack = src_f.get_basecall_events_pack(st, gr);
                    dst_f.add_basecall_events(st, gr, ev_pack);
                }
                else if (src_f.have_basecall_events_unpack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto ev_ds = src_f.get_basecall_events_dataset(st, gr);
                    auto & ev = ev_ds.first;
                    auto & ev_params = ev_ds.second;
                    // sampling rate
                    auto cid_params = src_f.get_channel_id_params();
                    // basecall fq
                    if (not src_f.have_basecall_fastq(st, gr))
                    {
                        LOG_ABORT
                            << "missing fastq for basecall events: st=" << st << " gr=" << gr << std::endl;
                    }
                    auto sq = src_f.get_basecall_seq(st, gr);
                    // ed group
                    auto ed_gr = src_f.get_basecall_eventdetection_group(gr);
                    std::vector< EventDetection_Event > ed;
                    if (not ed_gr.empty())
                    {
                        ed = src_f.get_eventdetection_events(ed_gr);
                    }
                    auto ev_pack = src_f.pack_ev(ev_ds, sq, ed, ed_gr,
                                                 cid_params, p_model_state_bits);
                    if (check)
                    {
                        if (ed_gr.empty())
                        {
                            auto rs_ds = src_f.get_raw_samples_dataset("");
                            ed = src_f.unpack_implicit_ed(ev_pack, rs_ds);
                        }
                        auto ev_ds_unpack = src_f.unpack_ev(ev_pack, sq, ed, cid_params);
                        auto & ev_unpack = ev_ds_unpack.first;
                        auto & ev_params_unpack = ev_ds_unpack.second;
                        if (not (ev_params_unpack == ev_params))
                        {
                            LOG_ABORT
                                << "check_failed ev_params_unpack!=ev_params" << std::endl;
                        }
                        if (ev_unpack.size() != ev.size())
                        {
                            LOG_ABORT
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " ev_unpack.size=" << ev_unpack.size()
                                << " ev_orig.size=" << ev.size() << std::endl;
                        }
                        for (unsigned i = 0; i < ev_unpack.size(); ++i)
                        {
                            if (abs(ev_unpack[i].start - ev[i].start) > 1e-3
                                or abs(ev_unpack[i].length - ev[i].length) > 1e-3
                                or abs(ev_unpack[i].mean - ev[i].mean) > 1e-1
                                or abs(ev_unpack[i].stdv - ev[i].stdv) > 1e-1
                                //or ev_unpack[i].move != ev[i].move // allow workaround for invalid moves
                                or ev_unpack[i].model_state != ev[i].model_state)
                            {
                                LOG_ABORT
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
                                    << ")" << std::endl;
                            }
                        }
                    }
                    dst_f.add_basecall_events(st, gr, ev_pack);
                    std::ostringstream oss;
                    if (not ev_pack.rel_skip.empty())
                    {
                        oss
                            << "rel_skip_bits=" << ev_pack.rel_skip_params.at("avg_bits");
                    }
                    else
                    {
                        oss
                            << "skip_bits=" << ev_pack.skip_params.at("avg_bits")
                            << " len_bits=" << ev_pack.len_params.at("avg_bits");
                    }
                    LOG(info)
                        << "gr=" << gr
                        << " st=" << st
                        << " ev_size=" << ev.size()
                        << " " << oss.str()
                        << " move_bits=" << ev_pack.move_params.at("avg_bits")
                        << " p_model_state_bits=" << ev_pack.p_model_state_params.at("num_bits")
                        << std::endl;
                }
            }
        }
    } // pack_ev()

    void
    unpack_ev(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        for (unsigned st = 0; st < 2; ++st)
        {
            auto gr_l = src_f.get_basecall_strand_group_list(st);
            for (auto const & gr : gr_l)
            {
                if (src_f.have_basecall_events(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto ev_ds = src_f.get_basecall_events_dataset(st, gr);
                    dst_f.add_basecall_events_dataset(st, gr, ev_ds);
                }
            }
        }
    } // unpack_ev()

    void
    copy_ev(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        for (unsigned st = 0; st < 2; ++st)
        {
            auto gr_l = src_f.get_basecall_strand_group_list(st);
            for (auto const & gr : gr_l)
            {
                if (src_f.have_basecall_events_unpack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto ev_ds = src_f.get_basecall_events_dataset(st, gr);
                    dst_f.add_basecall_events_dataset(st, gr, ev_ds);
                }
                else if (src_f.have_basecall_events_pack(st, gr))
                {
                    bc_gr_s.insert(gr);
                    auto ev_pack = src_f.get_basecall_events_pack(st, gr);
                    dst_f.add_basecall_events(st, gr, ev_pack);
                }
            }
        }
    } // copy_ev()

    void
    pack_al(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        auto gr_l = src_f.get_basecall_strand_group_list(2);
        for (auto const & gr : gr_l)
        {
            if (src_f.have_basecall_alignment_pack(gr))
            {
                bc_gr_s.insert(gr);
                auto al_pack = src_f.get_basecall_alignment_pack(gr);
                dst_f.add_basecall_alignment(gr, al_pack);
            }
            else if (src_f.have_basecall_alignment_unpack(gr))
            {
                bc_gr_s.insert(gr);
                auto al = src_f.get_basecall_alignment(gr);
                // basecall seq
                if (not src_f.have_basecall_seq(2, gr))
                {
                    LOG_ABORT
                        << "missing fastq for basecall alignment: gr=" << gr << std::endl;
                }
                auto seq = src_f.get_basecall_seq(2, gr);
                auto al_pack = src_f.pack_al(al, seq);
                if (check)
                {
                    auto al_unpack = src_f.unpack_al(al_pack, seq);
                    if (al_unpack.size() != al.size())
                    {
                        LOG_ABORT
                            << "check_failed gr=" << gr
                            << " al_unpack.size=" << al_unpack.size()
                            << " al_orig.size=" << al.size() << std::endl;
                    }
                    for (unsigned i = 0; i < al.size(); ++i)
                    {
                        if (al_unpack[i].template_index != al[i].template_index
                            or al_unpack[i].complement_index != al[i].complement_index
                            or al_unpack[i].get_kmer() != al[i].get_kmer())
                        {
                            LOG_ABORT
                                << "check_failed gr=" << gr
                                << " i=" << i
                                << " al_unpack=(" << al_unpack[i].template_index
                                << "," << al_unpack[i].complement_index
                                << "," << al_unpack[i].get_kmer()
                                << ") al_orig=(" << al[i].template_index
                                << "," << al[i].complement_index
                                << "," << al[i].get_kmer()
                                << ")" << std::endl;
                        }
                    }
                }
                dst_f.add_basecall_alignment(gr, al_pack);
                LOG(info)
                    << "gr=" << gr
                    << " al_size=" << al.size()
                    << " template_step_bits=" << al_pack.template_step_params.at("num_bits")
                    << " complement_step_bits=" << al_pack.complement_step_params.at("num_bits")
                    << " move_bits=" << al_pack.move_params.at("avg_bits")
                    << std::endl;
            }
        }
    } // pack_al()

    void
    unpack_al(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        auto gr_l = src_f.get_basecall_strand_group_list(2);
        for (auto const & gr : gr_l)
        {
            if (src_f.have_basecall_alignment(gr))
            {
                bc_gr_s.insert(gr);
                auto al = src_f.get_basecall_alignment(gr);
                dst_f.add_basecall_alignment(gr, al);
            }
        }
    } // unpack_al()

    void
    copy_al(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s) const
    {
        auto gr_l = src_f.get_basecall_strand_group_list(2);
        for (auto const & gr : gr_l)
        {
            if (src_f.have_basecall_alignment_unpack(gr))
            {
                bc_gr_s.insert(gr);
                auto al = src_f.get_basecall_alignment(gr);
                dst_f.add_basecall_alignment(gr, al);
            }
            else if (src_f.have_basecall_alignment_pack(gr))
            {
                bc_gr_s.insert(gr);
                auto al_pack = src_f.get_basecall_alignment_pack(gr);
                dst_f.add_basecall_alignment(gr, al_pack);
            }
        }
    } // copy_al()

    void
    copy_basecall_params(File const & src_f, File & dst_f, std::set< std::string > const & bc_gr_s) const
    {
        for (auto const & gr : bc_gr_s)
        {
            auto bc_params = src_f.get_basecall_params(gr);
            dst_f.add_basecall_params(gr, bc_params);
        }
    } // copy_basecall_params()

    void
    copy_attributes(File const & src_f, File const & dst_f, std::string const & p, bool recurse = false) const
    {
        File::Base::copy_attributes(src_f, dst_f, p, recurse);
    } // copy_attributes()
}; // class File_Packer

} // namespace fast5

#endif
