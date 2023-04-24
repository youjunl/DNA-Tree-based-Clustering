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

   m.def("new_tree", &new_trie, "Construct tree");
   m.def("search", &search, "Search in the tree");
   m.def("quick_search", &quick_search, "Search in the tree");
   m.def("insert", &insert_string, "Insert a string to the tree");
   m.def("delete", &delete_string, "Delete a string from the tree");
}
