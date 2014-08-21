#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
 
//<deque>, <list>, <map>, <set>, <queue>, <stack>, <vector>
//<vector>: .size(), .push_back(), .pop_back()
//int a[7] vs int* p = new int[7]
using namespace std;

#include <boost/algorithm/string.hpp>


std::string query_distance(const char* libfile, const char* pairs, const char* values) {
	//std::string int_string = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15";
	//std::vector<int> int_list;
	//strtk::parse(int_string,",",int_list);
	
	//std::string original("This,is,a,test,for,splitting.");	
	// std::string original(pairs);
	//std::vector<std::string> split; 
	//boost::split(split, original, boost::algorithm::is_any_of(","));	
	// std::copy(split.begin(), split.end(), std::ostream_iterator<std::string>(std::cout, " "));	
	cout << "Hello " <<  libfile << pairs << values << "!\n";
	// return original;

	std::string line;
	ifstream infile;
	infile.open(libfile);	
	while(!infile.eof()) {
		getline(infile, line);
		std::vector<std::string> split; 
		boost::split(split, line, boost::algorithm::is_any_of(","));	
		std::copy(split.begin(), split.end(), std::ostream_iterator<std::string>(std::cout, " "));	
		cout << "\n";
	}
	infile.close();

	return "out";
}
 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(statium_cpp)
{
	def("query_distance", query_distance);
}
