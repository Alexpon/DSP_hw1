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
   double gamma[100][100];
} HMM_Params;

void get_alpha_beta(HMM *, HMM_Params *, string);

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
			get_alpha_beta(&hmm, &hmm_params, sequence);
			break;
			//get_beta(sequence)
			//get_gamma()
		}
		ifs.close();
	}
	for (int j=0; j<hmm_params.sequence_size; j++){
		for (int i=0; i<hmm_params.state_num; i++){
			cout << hmm_params.beta[i][j] << " ";
		}
		cout << endl;
	}
	dumpHMM(open_or_die(output_file, "w"), &hmm);
	return 0;
}

void get_alpha_beta(HMM *hmm, HMM_Params *hmm_params, string sequence){
	double tmp_alpha, tmp_beta;
	hmm_params->state_num = hmm->state_num;
	hmm_params->sequence_size = sequence.size();

	// initialization
	for (int state=0; state<hmm_params->state_num; state++){
		hmm_params->alpha[state][0] = (hmm->initial[state]) * (hmm->observation[sequence[0]-65][state]);
		hmm_params->beta[state][hmm_params->sequence_size-1] = 1.0;
	}
	
	for (int observ=1; observ<hmm_params->sequence_size; observ++){
		for (int state=0; state<hmm_params->state_num; state++){
			tmp_alpha = 0.0;
			for (int sub_state=0; sub_state<hmm_params->state_num; sub_state++){
				tmp_alpha += hmm_params->alpha[sub_state][observ-1]*(hmm->transition[sub_state][state]);
			}
			hmm_params->alpha[state][observ] = tmp_alpha*(hmm->observation[sequence[observ]-65][state]);
		}
	}
	
	for (int observ=hmm_params->sequence_size-2; observ>=0; observ--){
		for (int state=0; state<hmm_params->state_num; state++){
			tmp_beta = 0.0;
			for (int sub_state=0; sub_state<hmm_params->state_num; sub_state++){
				tmp_beta += hmm->transition[state][sub_state]*hmm->observation[sequence[observ+1]-65][sub_state]*hmm_params->beta[sub_state][observ+1];
			}
			hmm_params->beta[state][observ] = tmp_beta;
		}
	}
}

