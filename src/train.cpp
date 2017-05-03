#include <iostream>
#include <string.h>
#include <fstream>
#include "hmm.h"
using namespace std;

#define STATE_NUM 6
#define OBSERV_NUM 6
#define SEQ_NUM 100

typedef struct{
	int sequence_size;
	int state_num;
	double alpha[STATE_NUM][SEQ_NUM];
	double beta[STATE_NUM][SEQ_NUM];
	double gamma[STATE_NUM][SEQ_NUM];
	double epsilon[SEQ_NUM][STATE_NUM][STATE_NUM];
} HMM_Params;

typedef struct{
	int counter;
	double pi[STATE_NUM];
	double transition_epsilon[STATE_NUM][STATE_NUM];
	double transition_gamma[STATE_NUM][STATE_NUM];
	double observation_numerator[OBSERV_NUM][STATE_NUM];
	double observation_denominator[OBSERV_NUM][STATE_NUM];
} HMM_Cumulate;

void hmm_initial(HMM_Cumulate *);
void calculate_alpha_beta(HMM *, HMM_Params *, string);
void calculate_gamma(HMM_Params *);
void calculate_epsilon(HMM *, HMM_Params *, string);
void cumulate_pi_A_B(HMM_Params *, HMM_Cumulate *, string);
void update_parameter(HMM *, HMM_Cumulate *);

/*-------------------------------------------------------------
typedef struct{
   char *model_name;
   int state_num;								//number of state
   int observ_num;								//number of observation
   double initial[STATE_NUM];					//initial prob.
   double transition[STATE_NUM][STATE_NUM];		//transition prob.
   double observation[OBSERV_NUM][STATE_NUM];	//observation prob.
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
	HMM_Cumulate hmm_cumulate;
	loadHMM(&hmm, initial_model);
	
	for (int i=0; i<iteration; i++){
		cout << "Iter: " << i+1 << endl;
		hmm_initial(&hmm_cumulate);
		ifstream ifs(observe_sequence, ifstream::in);
		while(getline(ifs, sequence)){
			calculate_alpha_beta(&hmm, &hmm_params, sequence);
			calculate_gamma(&hmm_params);
			calculate_epsilon(&hmm, &hmm_params, sequence);
			cumulate_pi_A_B(&hmm_params, &hmm_cumulate, sequence);
		}
		update_parameter(&hmm, &hmm_cumulate);
		ifs.close();
	}
	/*
	for (int j=0; j<hmm_params.sequence_size; j++){
		for (int i=0; i<hmm_params.state_num; i++){
			cout << hmm_params.gamma[i][j] << " ";
		}
		cout << endl;
	}*/
	dumpHMM(open_or_die(output_file, "w"), &hmm);
	return 0;
}

void hmm_initial(HMM_Cumulate *hmm_cumulate){
	hmm_cumulate->counter = 0;
	for (int i=0; i<STATE_NUM; i++){
		for (int j=0; j<STATE_NUM; j++){
			hmm_cumulate->transition_epsilon[i][j] = 0.0;
			hmm_cumulate->transition_gamma[i][j] = 0.0;
		}
		for (int k=0; k<OBSERV_NUM; k++){
			hmm_cumulate->observation_numerator[k][i] = 0.0;
			hmm_cumulate->observation_denominator[k][i] = 0.0;
		}
		hmm_cumulate->pi[i] = 0.0;
	}
}

void calculate_alpha_beta(HMM *hmm, HMM_Params *hmm_params, string sequence){
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


void calculate_gamma(HMM_Params *hmm_params){
	double cumulate_ab=0.0;
	for (int observ=0; observ<hmm_params->sequence_size; observ++){
		cumulate_ab=0.0;
		for (int state=0; state<hmm_params->state_num; state++){
			cumulate_ab += hmm_params->alpha[state][observ]*hmm_params->beta[state][observ];
		}
		for (int state=0; state<hmm_params->state_num; state++){
			hmm_params->gamma[state][observ] = hmm_params->alpha[state][observ]*hmm_params->beta[state][observ]/cumulate_ab;
		}
	}
}

void calculate_epsilon(HMM *hmm, HMM_Params *hmm_params, string sequence){
	for (int observ=0; observ<hmm_params->sequence_size-1; observ++){
		double cumulate = 0.0;
		for (int state_i=0; state_i<hmm_params->state_num; state_i++){
			for (int state_j=0; state_j<hmm_params->state_num; state_j++){
				cumulate += (hmm_params->alpha[state_i][observ]) * 
							(hmm->transition[state_i][state_j]) *
							(hmm->observation[sequence[observ+1]-65][state_j]) * 
							(hmm_params->beta[state_i][observ+1]);
			}
		}
		for (int state_i=0; state_i<hmm_params->state_num; state_i++){
			for (int state_j=0; state_j<hmm_params->state_num; state_j++){
				hmm_params->epsilon[observ][state_i][state_j] = 
							(hmm_params->alpha[state_i][observ]) *
							(hmm->transition[state_i][state_j]) *
							(hmm->observation[sequence[observ+1]-65][state_j]) * 
							(hmm_params->beta[state_i][observ+1]) / 
							cumulate;
			}
		}
	}
}

void cumulate_pi_A_B(HMM_Params *hmm_params, HMM_Cumulate *hmm_cumulate, string sequence){
	double tmp_epsilon=0.0;
	double tmp_gamma=0.0;
	double tmp_numerator=0.0;
	double tmp_denominator=0.0;

	hmm_cumulate->counter++;

	// Pi and A(transition prob)
	for (int state_i=0; state_i<STATE_NUM; state_i++){
		hmm_cumulate->pi[state_i] += hmm_params->gamma[state_i][0];

		tmp_gamma = 0.0;
		for (int t=0; t<hmm_params->sequence_size-1; t++)
			tmp_gamma += hmm_params->gamma[state_i][t];
		for (int state_j=0; state_j<STATE_NUM; state_j++){
			tmp_epsilon = 0.0;
			for (int t=0; t<hmm_params->sequence_size-1; t++)
				tmp_epsilon += hmm_params->epsilon[t][state_i][state_j];
			hmm_cumulate->transition_epsilon[state_i][state_j] += tmp_epsilon;
			hmm_cumulate->transition_gamma[state_i][state_j] += tmp_gamma;
		}
	}

	// B(observation prob)
	for (int state_j=0; state_j<STATE_NUM; state_j++){
		tmp_denominator = 0.0;
		for (int t=0; t<hmm_params->sequence_size; t++)
			tmp_denominator += hmm_params->gamma[state_j][t];
		for (int observ_k=0; observ_k<OBSERV_NUM; observ_k++){
			tmp_numerator = 0.0;
			for (int t=0; t<hmm_params->sequence_size; t++){
				if ((sequence[t]-65)==observ_k)
					tmp_numerator += hmm_params->gamma[state_j][t];
			}
			hmm_cumulate->observation_numerator[observ_k][state_j] += tmp_numerator;
			hmm_cumulate->observation_denominator[observ_k][state_j] += tmp_denominator;
		}
	}
}

void update_parameter(HMM *hmm, HMM_Cumulate *hmm_cumulate){
	for (int state_i=0; state_i<STATE_NUM; state_i++){
		hmm->initial[state_i] = hmm_cumulate->pi[state_i]/hmm_cumulate->counter;
	}
	for (int state_i=0; state_i<STATE_NUM; state_i++){
		for (int state_j=0; state_j<STATE_NUM; state_j++){
			hmm->transition[state_i][state_j] = (hmm_cumulate->transition_epsilon[state_i][state_j])/
									(hmm_cumulate->transition_gamma[state_i][state_j]);
		}
	}
	for (int state_j=0; state_j<STATE_NUM; state_j++){
		for (int observ_k=0; observ_k<OBSERV_NUM; observ_k++){
			hmm->observation[observ_k][state_j] = (hmm_cumulate->observation_numerator[observ_k][state_j])/
													(hmm_cumulate->observation_denominator[observ_k][state_j]);
		}
	}
}