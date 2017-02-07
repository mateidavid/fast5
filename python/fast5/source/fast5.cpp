#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/overloads.hpp>

#include "fast5.hpp"

namespace bp = boost::python;

// member functions with default arguments
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_raw_samples_overloads, have_raw_samples,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_raw_samples_params_overloads, get_raw_samples_params,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_raw_int_samples_overloads, get_raw_int_samples,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_raw_samples_overloads, get_raw_samples,
    0, 1)
//
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_eventdetection_group_overloads, have_eventdetection_group,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_eventdetection_read_name_list_overloads, get_eventdetection_read_name_list,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_eventdetection_events_overloads, have_eventdetection_events,
    0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_eventdetection_params_overloads, get_eventdetection_params,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_eventdetection_events_params_overloads, get_eventdetection_events_params,
    0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_eventdetection_events_overloads, get_eventdetection_events,
    0, 2)
//
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_group_overloads, have_basecall_group,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_strand_group_overloads, have_basecall_strand_group,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_fastq_overlords, have_basecall_fastq,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_fastq_overlords, get_basecall_fastq,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_seq_overlords, have_basecall_seq,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_seq_overlords, get_basecall_seq,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_model_overlords, have_basecall_model,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_model_file_overlords, get_basecall_model_file,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_model_params_overlords, get_basecall_model_params,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_model_overlords, get_basecall_model,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_events_overlords, have_basecall_events,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_events_params_overlords, get_basecall_events_params,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_events_overlords, get_basecall_events,
    1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    have_basecall_alignment_overlords, have_basecall_alignment,
    0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
    get_basecall_alignment_overlords, get_basecall_alignment,
    0, 1)

BOOST_PYTHON_MODULE(fast5)
{
    bp::class_<fast5::Channel_Id_Params>("Channel_Id_Params")
        .def_readwrite("channel_number", &fast5::Channel_Id_Params::channel_number)
        .def_readwrite("digitisation", &fast5::Channel_Id_Params::digitisation)
        .def_readwrite("offset", &fast5::Channel_Id_Params::offset)
        .def_readwrite("range", &fast5::Channel_Id_Params::range)
        .def_readwrite("sampling_rate", &fast5::Channel_Id_Params::sampling_rate)
        ;
    bp::class_<fast5::Raw_Samples_Params>("Raw_Samples_Params")
        .def_readwrite("read_id", &fast5::Raw_Samples_Params::read_id)
        .def_readwrite("read_number", &fast5::Raw_Samples_Params::read_number)
        .def_readwrite("start_mux", &fast5::Raw_Samples_Params::start_mux)
        .def_readwrite("start_time", &fast5::Raw_Samples_Params::start_time)
        .def_readwrite("duration", &fast5::Raw_Samples_Params::duration)
        ;;
    bp::class_<fast5::EventDetection_Events_Params>("EventDetection_Events_Params")
        .def_readwrite("read_id", &fast5::EventDetection_Events_Params::read_id)
        .def_readwrite("read_number", &fast5::EventDetection_Events_Params::read_number)
        .def_readwrite("scaling_used", &fast5::EventDetection_Events_Params::scaling_used)
        .def_readwrite("start_mux", &fast5::EventDetection_Events_Params::start_mux)
        .def_readwrite("start_time", &fast5::EventDetection_Events_Params::start_time)
        .def_readwrite("duration", &fast5::EventDetection_Events_Params::duration)
        .def_readwrite("median_before", &fast5::EventDetection_Events_Params::median_before)
        .def_readwrite("abasic_found", &fast5::EventDetection_Events_Params::abasic_found)
        ;
    bp::class_<fast5::EventDetection_Event>("EventDetection_Event")
        .def_readwrite("mean", &fast5::EventDetection_Event::mean)
        .def_readwrite("stdv", &fast5::EventDetection_Event::stdv)
        .def_readwrite("start", &fast5::EventDetection_Event::start)
        .def_readwrite("length", &fast5::EventDetection_Event::length)
        ;
    bp::class_<fast5::Basecall_Model_State>("Basecall_Model_State")
        .def_readwrite("variant", &fast5::Basecall_Model_State::variant)
        .def_readwrite("level_mean", &fast5::Basecall_Model_State::level_mean)
        .def_readwrite("level_stdv", &fast5::Basecall_Model_State::level_stdv)
        .def_readwrite("sd_mean", &fast5::Basecall_Model_State::sd_mean)
        .def_readwrite("sd_stdv", &fast5::Basecall_Model_State::sd_stdv)
        .def_readwrite("weight", &fast5::Basecall_Model_State::weight)
        .def_readwrite("kmer", &fast5::Basecall_Model_State::kmer)
        ;
    bp::class_<fast5::Basecall_Model_Params>("Basecall_Model_Params")
        .def_readwrite("scale", &fast5::Basecall_Model_Params::scale)
        .def_readwrite("shift", &fast5::Basecall_Model_Params::shift)
        .def_readwrite("drift", &fast5::Basecall_Model_Params::drift)
        .def_readwrite("var", &fast5::Basecall_Model_Params::var)
        .def_readwrite("scale_sd", &fast5::Basecall_Model_Params::scale_sd)
        .def_readwrite("var_sd", &fast5::Basecall_Model_Params::var_sd)
        ;
    bp::class_<fast5::Basecall_Event>("Basecall_Event")
        .def_readwrite("mean", &fast5::Basecall_Event::mean)
        .def_readwrite("stdv", &fast5::Basecall_Event::stdv)
        .def_readwrite("start", &fast5::Basecall_Event::start)
        .def_readwrite("length", &fast5::Basecall_Event::length)
        .def_readwrite("p_model_state", &fast5::Basecall_Event::p_model_state)
        .def_readwrite("p_mp_state", &fast5::Basecall_Event::p_mp_state)
        .def_readwrite("p_A", &fast5::Basecall_Event::p_A)
        .def_readwrite("p_C", &fast5::Basecall_Event::p_C)
        .def_readwrite("p_G", &fast5::Basecall_Event::p_G)
        .def_readwrite("p_T", &fast5::Basecall_Event::p_T)
        .def_readwrite("move", &fast5::Basecall_Event::move)
        .def_readwrite("model_state", &fast5::Basecall_Event::model_state)
        .def_readwrite("mp_state", &fast5::Basecall_Event::mp_state)
        ;;
    bp::class_<fast5::Basecall_Alignment_Entry>("Basecall_Alignment_Entry")
        .def_readwrite("template_index", &fast5::Basecall_Alignment_Entry::template_index)
        .def_readwrite("complement_index", &fast5::Basecall_Alignment_Entry::complement_index)
        .def("get_kmer", &fast5::Basecall_Alignment_Entry::get_kmer)
        ;;

    bp::class_<std::map<std::string, std::string>>("Map_Str_Str")
        .def(bp::map_indexing_suite<std::map<std::string, std::string>>())
        ;
    bp::class_<std::vector<std::string>>("Vec_Str")
        .def(bp::vector_indexing_suite<std::vector<std::string>>())
        ;
    bp::class_<std::vector<fast5::Raw_Int_Sample>>("Vec_Raw_Int_Sample")
        .def(bp::vector_indexing_suite<std::vector<fast5::Raw_Int_Sample>>())
        ;
    bp::class_<std::vector<fast5::Raw_Sample>>("Vec_Raw_Sample")
        .def(bp::vector_indexing_suite<std::vector<fast5::Raw_Sample>>())
        ;
    bp::class_<std::vector<fast5::EventDetection_Event>>("Vec_EventDetection_Event")
        .def(bp::vector_indexing_suite<std::vector<fast5::EventDetection_Event>>())
        ;
    bp::class_<std::vector<fast5::Basecall_Model_State>>("Vec_Basecall_Model_State")
        .def(bp::vector_indexing_suite<std::vector<fast5::Basecall_Model_State>>())
        ;
    bp::class_<std::vector<fast5::Basecall_Event>>("Vec_Basecall_Event")
        .def(bp::vector_indexing_suite<std::vector<fast5::Basecall_Event>>())
        ;
    bp::class_<std::vector<fast5::Basecall_Alignment_Entry>>("Vec_Basecall_Alignment_Entry")
        .def(bp::vector_indexing_suite<std::vector<fast5::Basecall_Alignment_Entry>>())
        ;

    bp::class_<fast5::File, boost::noncopyable>("File")
        .def(bp::init<std::string, bp::optional<bool>>())
        .def("is_open", &fast5::File::is_open)
        .def("is_rw", &fast5::File::is_rw)
        .def("file_name", &fast5::File::file_name, bp::return_value_policy<bp::copy_const_reference>())
        .def("open", &fast5::File::open)
        .def("create", &fast5::File::create)
        .def("close", &fast5::File::close)
        .def("is_valid_file", &hdf5_tools::File::is_valid_file).staticmethod("is_valid_file")
        .def("get_object_count", &hdf5_tools::File::get_object_count).staticmethod("get_object_count")
        //
        .def("file_version", &fast5::File::file_version)
        //
        .def("have_channel_id_params", &fast5::File::have_channel_id_params)
        .def("get_channel_id_params", &fast5::File::get_channel_id_params)
        //
        .def("have_sampling_rate", &fast5::File::have_sampling_rate)
        .def("get_sampling_rate", &fast5::File::get_sampling_rate)
        //
        .def("have_tracking_id_params", &fast5::File::have_tracking_id_params)
        .def("get_tracking_id_params", &fast5::File::get_tracking_id_params)
        //
        .def("have_sequences_params", &fast5::File::have_sequences_params)
        .def("get_sequences_params", &fast5::File::get_sequences_params)
        //
        .def("get_raw_samples_read_name_list",
             &fast5::File::get_raw_samples_read_name_list,
             bp::return_value_policy<bp::copy_const_reference>())
        .def("have_raw_samples",
             &fast5::File::have_raw_samples,
             have_raw_samples_overloads())
        .def("get_raw_samples_params",
             &fast5::File::get_raw_samples_params,
             get_raw_samples_params_overloads())
        .def("get_raw_int_samples",
             &fast5::File::get_raw_int_samples,
             get_raw_int_samples_overloads())
        .def("get_raw_samples",
             &fast5::File::get_raw_samples,
             get_raw_samples_overloads())
        //
        .def("get_eventdetection_group_list",
             &fast5::File::get_eventdetection_group_list,
             bp::return_value_policy<bp::copy_const_reference>())
        .def("have_eventdetection_group",
             &fast5::File::have_eventdetection_group)
        .def("get_eventdetection_read_name_list",
             &fast5::File::get_eventdetection_read_name_list,
             bp::return_value_policy<bp::copy_const_reference>(),
             get_eventdetection_read_name_list_overloads())
        .def("have_eventdetection_events",
             &fast5::File::have_eventdetection_events,
             have_eventdetection_events_overloads())
        .def("get_eventdetection_params",
             &fast5::File::get_eventdetection_params,
             get_eventdetection_params_overloads())
        .def("get_eventdetection_events_params",
             &fast5::File::get_eventdetection_events_params,
             get_eventdetection_events_params_overloads())
        .def("get_eventdetection_events",
             &fast5::File::get_eventdetection_events,
             get_eventdetection_events_overloads())
        //
        .def("get_basecall_group_list",
             &fast5::File::get_basecall_group_list,
             bp::return_value_policy<bp::copy_const_reference>())
        .def("have_basecall_group",
             &fast5::File::have_basecall_group)
        .def("get_basecall_strand_group_list",
             &fast5::File::get_basecall_strand_group_list,
             bp::return_value_policy<bp::copy_const_reference>())
        .def("have_basecall_strand_group",
             &fast5::File::have_basecall_strand_group)
        .def("get_basecall_1d_group",
             &fast5::File::get_basecall_1d_group,
             bp::return_value_policy<bp::copy_const_reference>())
        .def("get_basecall_eventdetection_group",
             &fast5::File::get_basecall_eventdetection_group,
             bp::return_value_policy<bp::copy_const_reference>())
        //
        .def("get_basecall_params",
             &fast5::File::get_basecall_params)
        //
        .def("have_basecall_log",
             &fast5::File::have_basecall_log)
        .def("get_basecall_log",
             &fast5::File::get_basecall_log)
        //
        .def("have_basecall_fastq",
             &fast5::File::have_basecall_fastq,
             have_basecall_fastq_overlords())
        .def("get_basecall_fastq",
             &fast5::File::get_basecall_fastq,
             get_basecall_fastq_overlords())
        .def("have_basecall_seq",
             &fast5::File::have_basecall_seq,
             have_basecall_seq_overlords())
        .def("get_basecall_seq",
             &fast5::File::get_basecall_seq,
             get_basecall_seq_overlords())
        //
        .def("have_basecall_model",
             &fast5::File::have_basecall_model,
             have_basecall_model_overlords())
        .def("get_basecall_model_file",
             &fast5::File::get_basecall_model_file,
             get_basecall_model_file_overlords())
        .def("get_basecall_model_params",
             &fast5::File::get_basecall_model_params,
             get_basecall_model_params_overlords())
        .def("get_basecall_model",
             &fast5::File::get_basecall_model,
             get_basecall_model_overlords())
        //
        .def("have_basecall_events",
             &fast5::File::have_basecall_events,
             have_basecall_events_overlords())
        .def("get_basecall_events_params",
             &fast5::File::get_basecall_events_params,
             get_basecall_events_params_overlords())
        .def("get_basecall_events",
             &fast5::File::get_basecall_events,
             get_basecall_events_overlords())
        //
        .def("have_basecall_alignment",
             &fast5::File::have_basecall_alignment,
             have_basecall_alignment_overlords())
        .def("get_basecall_alignment",
             &fast5::File::get_basecall_alignment,
             get_basecall_alignment_overlords())
        ;
}
