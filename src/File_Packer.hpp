#ifndef __FILE_PACKER_HPP
#define __FILE_PACKER_HPP

#include <string>
#include <set>

#include "fast5.hpp"

#include "logger.hpp"


namespace fast5
{

struct File_Packer
{
    struct opts
    {
        static bool & check() { static bool _check = true; return _check; }
        static unsigned const & max_qv_bits() { static unsigned const _max_qv_bits = 5; return _max_qv_bits; }
        static std::uint8_t const & max_qv() { static std::uint8_t const _max_qv = ((std::uint8_t)1 << max_qv_bits()) - 1; return _max_qv; }
        static unsigned & qv_bits() { static unsigned _qv_bits = max_qv_bits(); return _qv_bits; }
        static unsigned & p_model_state_bits() { static unsigned _p_model_state_bits = 2; return _p_model_state_bits; }
    };

    static void
    pack_rw(File const & src_f, File & dst_f)
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
                if (opts::check())
                {
                    auto rsi_ds_unpack = src_f.unpack_rw(rs_pack);
                    auto & rsi_unpack = rsi_ds_unpack.first;
                    auto & rs_params_unpack = rsi_ds_unpack.second;
                    if (not (rs_params_unpack == rs_params))
                    {
                        LOG(error)
                            << "check_failed rs_params_unpack!=rs_params" << std::endl;
                        abort();
                    }
                    if (rsi_unpack.size() != rsi.size())
                    {
                        LOG(error)
                            << "check_failed rs_unpack.size=" << rsi_unpack.size()
                            << " rs_orig.size=" << rsi.size() << std::endl;
                        abort();
                    }
                    for (unsigned i = 0; i < rsi_unpack.size(); ++i)
                    {
                        if (rsi_unpack[i] != rsi[i])
                        {
                            LOG(error)
                                << "check_failed i=" << i
                                << " rs_unpack=" << rsi_unpack[i]
                                << " rs_orig=" << rsi[i] << std::endl;
                            abort();
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

    static void
    unpack_rw(File const & src_f, File & dst_f)
    {
        auto rn_l = src_f.get_raw_samples_read_name_list();
        for (auto const & rn : rn_l)
        {
            auto rsi_ds = src_f.get_raw_int_samples_dataset(rn);
            dst_f.add_raw_samples_dataset(rn, rsi_ds);
        }
    } // unpack_rw()

    static void
    copy_rw(File const & src_f, File & dst_f)
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

    static void
    pack_ed(File const & src_f, File & dst_f)
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
                    if (opts::check())
                    {
                        auto rs_ds = src_f.get_raw_samples_dataset(rn);
                        auto ede_ds_unpack = src_f.unpack_ed(ede_pack, rs_ds);
                        auto & ede_unpack = ede_ds_unpack.first;
                        auto & ede_params_unpack = ede_ds_unpack.second;
                        if (not (ede_params_unpack == ede_params))
                        {
                            LOG(error)
                                << "check_failed ede_params_unpack!=ede_params" << std::endl;
                            abort();
                        }
                        if (ede_unpack.size() != ede.size())
                        {
                            LOG(error)
                                << "check_failed gr=" << gr
                                << " ede_unpack.size=" << ede_unpack.size()
                                << " ede_orig.size=" << ede.size() << std::endl;
                            abort();
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
                                LOG(error)
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
                                abort();
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

    static void
    unpack_ed(File const & src_f, File & dst_f)
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

    static void
    copy_ed(File const & src_f, File & dst_f)
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

    static void
    pack_fq(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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
                    auto fq_pack = src_f.pack_fq(fq, opts::qv_bits());
                    if (opts::check())
                    {
                        auto fq_unpack = src_f.unpack_fq(fq_pack);
                        auto fqa_unpack = src_f.split_fq(fq_unpack);
                        if (fqa_unpack[0] != fqa[0])
                        {
                            LOG(error)
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " fq_unpack_name=" << fqa_unpack[0]
                                << " fq_orig_name=" << fqa[0] << std::endl;
                            abort();
                        }
                        if (fqa_unpack[1] != fqa[1])
                        {
                            LOG(error)
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " fq_unpack_bp=" << fqa_unpack[1]
                                << " fq_orig_bp=" << fqa[1] << std::endl;
                            abort();
                        }
                        if (fqa_unpack[3].size() != fqa[3].size())
                        {
                            LOG(error)
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " fq_unpack_qv_size=" << fqa_unpack[3].size()
                                << " fq_orig_qv_size=" << fqa[3].size() << std::endl;
                            abort();
                        }
                        auto qv_mask = opts::max_qv() & (opts::max_qv() << (opts::max_qv_bits() - opts::qv_bits()));
                        for (unsigned i = 0; i < fqa_unpack[3].size(); ++i)
                        {
                            if ((std::min<std::uint8_t>(fqa_unpack[3][i] - 33, opts::max_qv()) & qv_mask) !=
                                (std::min<std::uint8_t>(fqa[3][i] - 33, opts::max_qv()) & qv_mask))
                            {
                                LOG(error)
                                    << "check_failed st=" << st
                                    << " gr=" << gr
                                    << " i=" << i
                                    << " fq_unpack_qv=" << fqa_unpack[3][i]
                                    << " fq_orig_qv=" << fqa[3][i] << std::endl;
                                abort();
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

    static void
    unpack_fq(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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

    static void
    copy_fq(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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

    static void
    pack_ev(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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
                    auto channel_id_params = src_f.get_channel_id_params();
                    // basecall fq
                    if (not src_f.have_basecall_fastq(st, gr))
                    {
                        LOG(error)
                            << "missing fastq for basecall events: st=" << st << " gr=" << gr << std::endl;
                        abort();
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
                                                 channel_id_params.sampling_rate, opts::p_model_state_bits());
                    if (opts::check())
                    {
                        if (ed_gr.empty())
                        {
                            auto rs_ds = src_f.get_raw_samples_dataset("");
                            ed = src_f.unpack_implicit_ed(ev_pack, rs_ds);
                        }
                        auto ev_ds_unpack = src_f.unpack_ev(ev_pack, sq, ed, channel_id_params.sampling_rate);
                        auto & ev_unpack = ev_ds_unpack.first;
                        auto & ev_params_unpack = ev_ds_unpack.second;
                        if (not (ev_params_unpack == ev_params))
                        {
                            LOG(error)
                                << "check_failed ev_params_unpack!=ev_params" << std::endl;
                            abort();
                        }
                        if (ev_unpack.size() != ev.size())
                        {
                            LOG(error)
                                << "check_failed st=" << st
                                << " gr=" << gr
                                << " ev_unpack.size=" << ev_unpack.size()
                                << " ev_orig.size=" << ev.size() << std::endl;
                            abort();
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
                                    << ")" << std::endl;
                                abort();
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

    static void
    unpack_ev(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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

    static void
    copy_ev(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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

    static void
    pack_al(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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
                    LOG(error)
                        << "missing fastq for basecall alignment: gr=" << gr << std::endl;
                    abort();
                }
                auto seq = src_f.get_basecall_seq(2, gr);
                auto al_pack = src_f.pack_al(al, seq);
                if (opts::check())
                {
                    auto al_unpack = src_f.unpack_al(al_pack, seq);
                    if (al_unpack.size() != al.size())
                    {
                        LOG(error)
                            << "check_failed gr=" << gr
                            << " al_unpack.size=" << al_unpack.size()
                            << " al_orig.size=" << al.size() << std::endl;
                        abort();
                    }
                    for (unsigned i = 0; i < al.size(); ++i)
                    {
                        if (al_unpack[i].template_index != al[i].template_index
                            or al_unpack[i].complement_index != al[i].complement_index
                            or al_unpack[i].get_kmer() != al[i].get_kmer())
                        {
                            LOG(error)
                                << "check_failed gr=" << gr
                                << " i=" << i
                                << " al_unpack=(" << al_unpack[i].template_index
                                << "," << al_unpack[i].complement_index
                                << "," << al_unpack[i].get_kmer()
                                << ") al_orig=(" << al[i].template_index
                                << "," << al[i].complement_index
                                << "," << al[i].get_kmer()
                                << ")" << std::endl;
                            abort();
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

    static void
    unpack_al(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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

    static void
    copy_al(File const & src_f, File & dst_f, std::set< std::string > & bc_gr_s)
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

    static void
    copy_basecall_params(File const & src_f, File & dst_f, std::set< std::string > const & bc_gr_s)
    {
        for (auto const & gr : bc_gr_s)
        {
            auto bc_params = src_f.get_basecall_params(gr);
            dst_f.add_basecall_params(gr, bc_params);
        }
    } // copy_basecall_params()

    static void
    copy_attributes(File const & src_f, File const & dst_f, std::string const & p, bool recurse = false)
    {
        File::Base::copy_attributes(src_f, dst_f, p, recurse);
    } // copy_attributes()
}; // class File_Packer

} // namespace fast5

#endif
