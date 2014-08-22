#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <boost/algorithm/string.hpp>
 
//<deque>, <list>, <map>, <set>, <queue>, <stack>, <vector>
//<vector>: .size(), .push_back(), .pop_back()
//int a[7] vs int* p = new int[7]

using namespace std;

vector<double> convertStringVectortoDoubleVector(const vector<string>& stringVector) {
	vector<double> doubleVector(stringVector.size());
	transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const string& val)
		{
			return stod(val);
		});
	return doubleVector;
}


string double_array_to_string(double double_array[], int size_of_array) {
	string s = "";
	for (int i = 0; i < size_of_array; i++) {
		string temp;
		ostringstream convert;
		convert << double_array[i];
		s += convert.str() + ", ";
	}
	s.erase(s.end()-2, s.end());
	return s;
}

void print_double_vector(vector<double>& v) {
	copy(v.begin(), v.end(), std::ostream_iterator<double>(cout, " "));	
	cout << "" << endl;
}

void print_int_vector(vector<int>& v) {
	copy(v.begin(), v.end(), std::ostream_iterator<int>(cout, " "));	
	cout << "" << endl;
}

void print_array (int arg[], int length) {
	for (int n=0; n<length; ++n)
		cout << arg[n] << ' ';
	cout << '\n';
}

void print_array (double arg[], int length) {
	for (int n=0; n<length; ++n)
		cout << arg[n] << ' ';
	cout << '\n';
}

//Ignore threshold_adjust and pair_check for now
string query_distance(const char* libfile, const char* values, int threshold_adjust_counts=100, float threshold=0.4) {

	bool threshold_adjust = threshold_adjust_counts == -1 ? false:true;
	
	//Splitting template distance values into vector
	string dist_str(values);
	vector<string> dist_split_str; 
	boost::split(dist_split_str, dist_str, boost::algorithm::is_any_of(","));	
	vector<double> dists = convertStringVectortoDoubleVector(dist_split_str);
	
	int num_dists = dists.size();

	//RMSD values and AA identities for every IP in library file
	vector<int> aas;
	vector<double> rmsds;

	//Keep track of lowest 100 RMSDs
	vector<double> thresholds;

	string line;
	ifstream infile;
	infile.open(libfile);	
	while (!infile.eof()) {
		getline(infile, line);
		vector<std::string> split1; 
		boost::split(split1, line, boost::algorithm::is_any_of(";"));
		if (split1.size() > 1) {
			string lib_dists_str(split1.at(0));
			vector<std::string> split2; 
			boost::split(split2, lib_dists_str, boost::algorithm::is_any_of(","));
			vector<double> lib_dists = convertStringVectortoDoubleVector(split2);

			float rmsd = 0;
			for (int i = 0; i < num_dists; i++) {
				rmsd += (lib_dists.at(i) - dists.at(i)) * (lib_dists.at(i) - dists.at(i));
			}
			rmsd = sqrt(rmsd/(num_dists-1));
			rmsds.push_back(rmsd);

			if (threshold_adjust) {
				if (thresholds.size() < threshold_adjust_counts) {
					thresholds.push_back(rmsd);
					push_heap(thresholds.begin(), thresholds.end());

				} else if (rmsd < thresholds.front()) {
					pop_heap (thresholds.begin(), thresholds.end());
					thresholds.pop_back();

					thresholds.push_back(rmsd);
					push_heap(thresholds.begin(), thresholds.end());
				}
			}

			int aa = stoi(split1.at(1));
			aas.push_back(aa);
		}
		cout << "\n";
	}

	const float ip_probs[20] = {0.064, 0.018, 0.038, 0.043, 0.065, 0.052, 0.025, 0.084, 0.04, 0.127, 0.024, 0.033, 0.038, 0.029, 0.052, 0.044, 0.05, 0.095, 0.025, 0.055};
	infile.close();

	threshold = threshold_adjust ? thresholds.front(): threshold;
	cout << "Using threshold " << threshold << endl;

	int match_counts[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double match_probs[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	//Tally up matches within RMSD threshold
	for (unsigned int i = 0; i < aas.size(); i++)
		if (rmsds[i] < threshold) match_counts[aas[i]] += 1;
	
	int match_total = 0;
	for (int i = 0; i < 20; i++) {
		if (match_counts[i] == 0) match_counts[i] += 1;
		match_total += match_counts[i];
	}

	for (int i = 0; i < 20; i++)
		match_probs[i] = match_counts[i]/match_total;
	
	double energies[20];
	for (int i = 0; i < 20; i++) {
		if (abs(match_probs[i] - ip_probs[i]) < 0.0001 || match_probs[i] == 0.0) {
			energies[i] = 0.0;
		} else {
			energies[i] = -1.0*log10(match_probs[i] / ip_probs[i]);
		}
	}

	return double_array_to_string(energies, 20);
}
 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(statium_cpp)
{
	def("query_distance", query_distance);
}
