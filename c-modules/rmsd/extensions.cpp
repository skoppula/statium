#include <iostream>
#include "strtk.hpp"
 
using namespace std;
 
void query_distance(const char* libfile, const char* pairs, const char* values) {
	//std::string int_string = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15";
	//std::vector<int> int_list;
	//strtk::parse(int_string,",",int_list);
	cout << "Hello " <<  libfile << pairs << values << "!\n";
}
 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(statium_cpp)
{
	def("query_distance", query_distance);
}
