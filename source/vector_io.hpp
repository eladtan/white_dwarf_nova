#ifndef VECTOR_IO_HPP
#define VECTOR_IO_HPP 1

#include <vector>
#include <string>
#include <cassert>

using std::vector;
using std::string;

vector<double> load_txt(const string& fname);

void save_txt(const string& fname, const vector<double>& v);

#endif // VECTOR_IO_HPP
