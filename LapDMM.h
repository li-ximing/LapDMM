/*
* Copyright (C) 2019 by
* 
* 	Ximing Li
*	liximing86@gmail.com
* 	College of Computer Science and Technology, JiLin University
*
*/

#ifndef	_LAPDMM_H_
#define	_LAPDMM_H_

#include "dataset.h"
#include "evaluation.h"

#define BETA 0.01;
#define ALPHA 0.1;

// collapsed variational inference (generalizaed EM) for Dirichlet mixture model with manifold regularization
class LapDMM
{
public:
	string fileBasis;

	dataset * trainData;	// pointer to training dataset object

	EvaluationTopic* evaluation;	// class for evaluation

	// --- parameters and variables ---    
	int M;							// number of documents
	int V;							// number of unique words
	int K;							// number of topics	
	
	// Dirichlet priors
	double* alpha;					// document-topic prior, size K
	double* beta;					// topic-word prior, size V
	double Kalpha;					// the sum of alpha
	double Vbeta;					// the sum of beta

	// Distributions	
	double* theta;					// document-topic multinomial distribution, size K
	double** phi;					// topic-word multinomial distribution, size K * V

	// interations
	int maxIteration;							// maximum iteration number
	double epsilon;								// convergence condition

	// collapsed variational inference
	double** q;									// variational distributions for topic assignments of documents, size M x K
	double** numOfTopicWord;					// soft number of words assigned to topics, size K x V
	double* totalNumOfTopic;;					// total soft number of words assigned to topics, size K
	double* numOfTopicDoc;						// soft number of docuemnts assigned to topics, size K

	// manifold regularization
	double gamma;								// step parameter in manifold update
	double lambda;								// tuning paramter for the manifold regularization term
	int numOfNN;								// number of the nearest neighbour for documents
	int** idxOfDocGraph;						// index matrix of neighbour doucments, size M * numOfNN
	int mrIteration;							// max iteration number of manifold optimization

	double** tmpq;
	double** tmpNumOfTopicWord;
	double* tmpTotalNumOfTopic;
	double* tmpNumOfTopicDoc;

	// flag
	int ManifoldRegularization;
	int GraphTF;

	// --------------------------------------

	LapDMM() 
	{
		SetDefaultValues();
	}

	~LapDMM();

	// set default values for variables
	void SetDefaultValues();   

	// init for inference or estimation
	int InitEst(string fileOfTrainData);

	int ReadDocGraph(string filename);

	// model training using CVB
	void Train();
	void qEstimate();					
	void qUpdate(int m);
	double LapDMM_lik(dataset* corpus);

	void ComputeTheta(double* numOfTopicDoc);
	void ComputePhi(double** numOfTopicWord, double* totalNumOfTopic);

	// manifold regularization
	void GraphConstructionTF();
	void GraphConstructionWMD();

	void UpdateTmpq();
	void CopyTmpq();

	double MR_lik(double** q);

	// topic representations of documents
	void ComputeTopicRepresentationSW();

	// common functions
	void NormalizeMatrix(double **dst, double **src ,int R, int L);

	void Evaluation();
};

#endif

