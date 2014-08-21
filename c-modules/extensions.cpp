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

//Ignore threshold_adjust and pair_check for now
std::string query_distance(const char* libfile, const char* pairs, const char* values, bool pair_check=false, bool threshold_adjust=true) {
	
	//Splitting template distance values into vector
	string dist_str(values);
	vector<string> dist_split_str; 
	boost::split(dist_split_str, dist_str, boost::algorithm::is_any_of(","));	
	vector<double> dists = convertStringVectortoDoubleVector(dist_split_str);
	//Printing out distance values
	copy(dists.begin(), dists.end(), std::ostream_iterator<double>(cout, " "));	

	char AAs[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	float ip_probs[] = {0.064, 0.018, 0.038, 0.043, 0.065, 0.052, 0.025, 0.084, 0.04, 0.127, 0.024, 0.033, 0.038, 0.029, 0.052, 0.044, 0.05, 0.095, 0.025, 0.055};

	float threshold_counts = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	float thresholds = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}

	int num_dists = dists.size()

	string line;
	ifstream infile;
	infile.open(libfile);	
	while(!infile.eof()) {
		getline(infile, line);
		vector<std::string> split1; 
		boost::split(split1, line, boost::algorithm::is_any_of(";"));
		if(split1.size() > 1) {
			string lib_dists_str(split1.at(1));
			vector<std::string> split2; 
			boost::split(split2, lib_dists_str, boost::algorithm::is_any_of(","));
			vector<double> lib_dists = convertStringVectortoDoubleVector(split2);
			copy(lib_dists.begin(), lib_dists.end(), std::ostream_iterator<double>(cout, " "));	
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
