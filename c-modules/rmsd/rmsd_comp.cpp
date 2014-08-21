#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main(int argc, char* argv[]) {
    string arg_name;
    vector<float> template_d;
    char* input;
    char* output;
    int rsize;
    int fsize;
    for(int i = 0; i < argc; i++) {
	arg_name = argv[i];
	if (arg_name == "-distances") {
	    for(int j = i + 5; j < argc; j++) {
	        float tmpf = atof(argv[j]);
		template_d.push_back(tmpf);
	    }
	    input = argv[i + 1];
	    output = argv[i + 2];
	    rsize = atoi(argv[i + 3]);
	    fsize = atoi(argv[i + 4]);
	}
    }

    cout << "Allocating memory for " << input << "..." << endl;
    char aa_dict [20] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    float ip_probs [20] = {0.064, 0.018, 0.038, 0.043, 0.065, 0.052, 0.025, 0.084, 0.04, 0.127, 0.024, 0.033, 0.038, 0.029, 0.052, 0.044, 0.05, 0.095, 0.025, 0.055};
    float lib_dists [rsize];
    float threshold_counts [10];
    float thresholds [10];
    float* rmsd_list = new float[fsize];
    int* aa_list = new int[fsize];
    float incr = 0.1;
    for(int i = 0; i < 10; i++) {
	threshold_counts[i] = 0;
	thresholds[i] = incr;
	incr += 0.1;
    }
    cout << "Reading distances for " << input << "..." << endl;
    ifstream file;
    file.open(input);
    string str;
    int count = -1;
    int aa;
    float dist;
    while (0==0)
    {
	if (file.eof()) break;

	if (count % 100000 == 0) cout << count << endl;
	count++;
	for(int i = 0; i < rsize; i++) {
	    file >> dist;
	    lib_dists[i] = dist;
	}
	
	file >> aa;
	float rmsd = 0.0;
	for(int i = 0; i < rsize; i++) {
	    rmsd = rmsd + ((lib_dists[i] - template_d[i]) * (lib_dists[i] - template_d[i]));
	}
	float libN = rsize - 1;
	rmsd = sqrt(rmsd / libN);
	
	float idxf = rmsd * 10.0;
	if (idxf > 10.0) idxf = 10.0;
	int idx = idxf;
        if (rmsd < thresholds[idx]) threshold_counts[idx] += 1;
	rmsd_list[count] = rmsd;
	aa_list[count] = aa;	
    }
    file.close();
    cout << "Calculating threshold for " << input << "..." << endl;
    int running_count = 0;
    float use_threshold;
    for(int i = 0; i < 10; i++) {
	running_count += threshold_counts[i];
	if (running_count > 100) {
	    use_threshold = thresholds[i]; 
	    break;
	}
    }
    cout << use_threshold << endl;;
    vector<float> match_counts;
    vector<float> match_probs;
    for(int i = 0; i < 20; i++) {
	match_counts.push_back(0.0);
	match_probs.push_back(0.0);
    }
    for(int i = 0; i < fsize; i++) {
         if (rmsd_list[i] < use_threshold) {
	     match_counts[aa_list[i]] += 1.0;
	 }
    }

    float match_total = 0.0;
    for(int i = 0; i < 20; i++) {
	if (match_counts[i] == 0.0) match_counts[i] += 1.0;
	match_total += match_counts[i];
    }
    for(int i = 0; i < 20; i++) {
	match_probs[i] = match_counts[i] / match_total;
    }
    cout << "Calculating scores for " << input << "..." << endl;
    ofstream outfile;
    outfile.open(output);
    float e;
    for(int i = 0; i < 20; i++) {
	if (match_probs[i] != ip_probs[i]) e = -1.0*log10(match_probs[i] / ip_probs[i]);
	else e = 0.0;
	if (e < 0.0001 && e > -0.0001) e = 0.0;
	outfile << aa_dict[i] << ' ' << e << endl;
    }
    outfile.close();

    delete[] rmsd_list;
    delete[] aa_list;

    cout << "Done with cpp for " << input << "..." << endl;
    return 0;
}
