#include <iostream>
#include <string.h>
#include <fstream>
#include "hmm.h"
using namespace std;

#define MAX_PREDICT_SIZE 100000

double viterbi(HMM *, string);
void get_acc(int *, int);

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
	if (argc != 4){
		cout << "Execute Error" << endl;
		cout << "You need to use format" << endl;
		cout << "./test modellist.txt testing_data.txt result.txt" << endl;
		return 0;
	}

	const char* modellist_filename = argv[1];
	const char* testing_data_filename = argv[2];
	const char* result_filename = argv[3];
	double optimal_prob = 0.0;
	double model_prob = 0.0;
	int optimal_model = 1;
	int predict[MAX_PREDICT_SIZE];
	int cnt = 0;
	HMM hmms[5];
	load_models( modellist_filename, hmms, 5);

	string sequence="";
	ifstream ifs(testing_data_filename, ifstream::in);
	ofstream ofs(result_filename, ofstream::out);
	while(getline(ifs, sequence)){
		optimal_prob = 0.0;
		optimal_model = 1;
		for (int i=0; i<5; i++){
			model_prob = viterbi(&hmms[i], sequence);
			if (model_prob>optimal_prob){
				optimal_prob = model_prob;
				optimal_model = i+1;
				predict[cnt] = i+1;
			}
		}
		ofs << "model0" << optimal_model << ".txt " << optimal_prob << endl;
		cnt ++;
	}
	ifs.close();	
	ofs.close();
	
	get_acc(predict, cnt);

	return 0;
}

double viterbi(HMM *hmm, string sequence){
	int sequence_size = sequence.size();
	double delta[hmm->state_num][sequence_size];
	double max_pre_state = 0.0;
	double max_total_prob = 0.0;
	
	for (int state_i=0; state_i<hmm->state_num; state_i++)
		delta[state_i][0] = hmm->initial[state_i]*hmm->observation[sequence[0]-65][state_i];

	// recursion
	for (int observ_k=1; observ_k<sequence_size; observ_k++){
		for (int state_i=0; state_i<hmm->state_num; state_i++){
			max_pre_state = 0.0;
			for (int pre_state=0; pre_state<hmm->state_num; pre_state++){
				if (delta[pre_state][observ_k-1]*hmm->transition[pre_state][state_i] > max_pre_state){
					max_pre_state = delta[pre_state][observ_k-1]*hmm->transition[pre_state][state_i];
				}
			}
			delta[state_i][observ_k] = max_pre_state * hmm->observation[sequence[observ_k]-65][state_i];
		}
	}

	for (int state_i=0; state_i<hmm->state_num; state_i++){
		if (delta[state_i][sequence_size-1]>max_total_prob)
			max_total_prob = delta[state_i][sequence_size-1];
	}
	return max_total_prob;
}

void get_acc(int predict[], int total){
	int cnt = 0;
	int match = 0;
	string line="";

	ifstream ifs("../testing_answer.txt", ifstream::in);
	ofstream ofs("acc.txt", ofstream::out);
	while (getline(ifs, line)){
		if ((line[7]-48)==predict[cnt]){
			match ++;
		}
		cnt ++;
	}
	if (cnt != total)
		cout << "The size of prediction is not equal to answer" << endl;
	else
		ofs << (match/(double)cnt) << endl;
	ofs.close();
}
