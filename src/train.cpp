#include <iostream>
#include <string.h>
#include <fstream>
#include "hmm.h"
using namespace std;


typedef struct{
   int sequence_size;
   int state_num;
   double alpha[100][100];
   double beta[100][100];
   double gamme[100][100];
} HMM_Params;

void get_alpha(HMM *, HMM_Params *, string);

/*-------------------------------------------------------------
typedef struct{
   char *model_name;
   int state_num;								//number of state
   int observ_num;								//number of observation
   double initial[MAX_STATE];					//initial prob.
   double transition[MAX_STATE][MAX_STATE];		//transition prob.
   double observation[MAX_OBSERV][MAX_STATE];	//observation prob.
} HMM;
----------------------------------------------------------------*/


int main(int argc, char **argv){
	if (argc != 5){
		cout << "Execute Error" << endl;
		cout << "You need to use format" << endl;
		cout << "./train iteration model_init.txt seq_model_01.txt model_01.txt" << endl;
		return 0;
	}
	int iteration = stoi(argv[1]);
	const char* initial_model = argv[2];
	const char* observe_sequence = argv[3];
	const char* output_file = argv[4];
	string sequence="";
	HMM hmm;
	HMM_Params hmm_params;
	loadHMM(&hmm, initial_model);
	
	for (int i=0; i<iteration; i++){
		ifstream ifs(observe_sequence, ifstream::in);
		while(getline(ifs, sequence)){
			get_alpha(&hmm, &hmm_params, sequence);
			break;
			//get_beta(sequence)
			//get_gamma()
		}
		ifs.close();
	}
	for (int i=0; i<hmm_params.state_num; i++){
		cout << hmm_params.alpha[i][0] <<endl;
	}
	dumpHMM(open_or_die(output_file, "w"), &hmm);
	return 0;
}

void get_alpha(HMM *hmm, HMM_Params *hmm_params, string sequence){
	hmm_params->state_num = hmm->state_num;
	hmm_params->sequence_size = sequence.size();

	// initialization
	for (int i=0; i<hmm_params->state_num; i++){
		hmm_params->alpha[i][0] = (hmm->initial[i]) * (hmm->observation[sequence[0]-65][i]);
	}
	// induction
	
}


