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

std::vector<double> convertStringVectortoDoubleVector(const std::vector<std::string>& stringVector) {
	std::vector<double> doubleVector(stringVector.size());
	std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const std::string& val)
		{
			return std::stod(val);
		});
	return doubleVector;
}

std::string query_distance(const char* libfile, const char* pairs, const char* values) {
	
	//Splitting template distance values into vector
	std::string dist_str(values);
	std::vector<std::string> dist_split_str; 
	boost::split(dist_split_str, dist_str, boost::algorithm::is_any_of(","));	
	std::vector<double> dists = convertStringVectortoDoubleVector(dist_split_str);
	//Printing out distance values
	std::copy(dists.begin(), dists.end(), std::ostream_iterator<double>(std::cout, " "));	

	std::string line;
	ifstream infile;
	infile.open(libfile);	
	while(!infile.eof()) {
		getline(infile, line);
		std::vector<std::string> split1; 
		boost::split(split1, line, boost::algorithm::is_any_of(";"));
		if(split1.size() > 1) {
			std::string lib_dists_str(split1.at(1));
			std::vector<std::string> split2; 
			boost::split(split2, lib_dists_str, boost::algorithm::is_any_of(","));
			std::vector<double> lib_dists = convertStringVectortoDoubleVector(split2);
			std::copy(lib_dists.begin(), lib_dists.end(), std::ostream_iterator<double>(std::cout, " "));	
		}
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
