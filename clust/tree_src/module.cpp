#include <pybind11/pybind11.h>
#include "trie.h"

PYBIND11_MODULE(tree, m)
{
    m.doc() = "Python binding for tree search algorithms";
    namespace py = pybind11;
    py::class_<trie_t>(m, "Tree");
    py::class_<result_t>(m, "Result")
        .def(py::init<>())
        .def_readwrite("label", &result_t::label)
        .def_readwrite("distance", &result_t::distance);

    py::class_<config_t>(m, "Config")
        .def(py::init<>())
        .def_readwrite("tree_depth", &config_t::tree_depth)
        .def_readwrite("index_begin", &config_t::index_begin)
        .def_readwrite("main_tree_threshold", &config_t::main_tree_threshold)
        .def_readwrite("sub_tree_threshold", &config_t::sub_tree_threshold)
        .def_readwrite("depth_limit", &config_t::depth_limit);

    py::class_<Process>(m, "Process")
        .def(py::init<config_t *>())
        .def("cluster", &Process::cluster);

    m.def("new_tree", &new_trie, "Construct tree");
    m.def("search", &poucet_search, "Search in the tree");
    m.def("quick_search", &quick_search, "Search in the tree");
    m.def("insert", &insert_string, "Insert a string to the tree");
    m.def("delete", &delete_string, "Delete a string from the tree");
}
