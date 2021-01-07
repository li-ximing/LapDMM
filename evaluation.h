/*
* Copyright (C) 2019 by
* 
* 	Ximing Li
*	liximing86@gmail.com
* 	College of Computer Science and Technology, JiLin University
*
*/

#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include <time.h>
#include "gamma.h"
#include "dataset.h"
#include "hungarian.h"
#include <iostream>  
#include <fstream>  
#include <string> 
#include <cstdlib> 
using namespace std;


// Evaluation for topics
class EvaluationTopic {
public:

	string filename;							// the filename of vocabulary
	int topword;								// the number of top words for each topic

	EvaluationTopic(string fileBasis);

	~EvaluationTopic();

	double ClusteringNMI(dataset* trainData, int* z, int K);

	double ClusteringACC(dataset* trainData, int* z, int K);
};

#endif