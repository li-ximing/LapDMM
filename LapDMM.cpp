/*
* Copyright (C) 2019 by
* 
* 	Ximing Li
*	liximing86@gmail.com
* 	College of Computer Science and Technology, JiLin University
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "dataset.h"
#include "LapDMM.h"

LapDMM::~LapDMM() 
{
/*
	if (trainData) 
	{
		delete trainData;
	}
*/

	if(alpha)
	{
		delete [] alpha;
	}

	if(beta)
	{
		delete [] beta;
	}

	if(q) 
	{
		for(int m=0 ; m<M ; m++)
		{
			if(q[m])
			{
				delete [] q[m];
			}
		}
		delete [] q;
	}

	if(numOfTopicWord) 
	{
		for(int k=0 ; k<K ; k++) 
		{
			if(numOfTopicWord[k]) 
			{
				delete [] numOfTopicWord[k];
			}
		}
		delete [] numOfTopicWord;
	}

	if(totalNumOfTopic) 
	{
		delete [] totalNumOfTopic;
	}  

	if(numOfTopicDoc)
	{
		delete [] numOfTopicDoc;
	}

	if(theta) 
	{
		delete [] theta;
	}

	if(phi) 
	{
		for(int k=0 ; k<K ; k++) 
		{
			if(phi[k]) 
			{
				delete [] phi[k];
			}
		}
		delete [] phi;
	}

	if(idxOfDocGraph)
	{
		for(int m=0 ; m<M ; m++)
		{
			if(idxOfDocGraph[m])
			{
				delete [] idxOfDocGraph[m];
			}
		}
		delete [] idxOfDocGraph;
	}
}

void LapDMM::SetDefaultValues() 
{   
	trainData = NULL;
	evaluation = NULL;

	M = 0;
	V = 0;
	K = 6;

	ManifoldRegularization = true;
	GraphTF = true;
	gamma = 0.1;
	lambda = 0.1;
	numOfNN = 9;
	mrIteration = 10;

	maxIteration = 100;
	epsilon = 1e-6;

	q = NULL;
	numOfTopicWord = NULL;
	totalNumOfTopic = NULL;
	numOfTopicDoc = NULL;
	theta = NULL;
	phi = NULL;
	idxOfDocGraph = NULL;
}

int LapDMM::InitEst(string fileOfTrainData) 
{
	//	randomize a seed
	srand((unsigned)time(NULL));

	// + read training data
	this->trainData = new dataset;

	if(trainData->InitDataset(fileOfTrainData) != 0)
	{
		printf("error in initialization of dataset");
		exit(1);
	}

	fileBasis = fileOfTrainData.substr(0,fileOfTrainData.length()-4);

	M = trainData->M;
	V = trainData->V;

	evaluation = new EvaluationTopic(fileBasis);

	//	initialize parameters	 
	alpha = new double[K];
	for(int k=0 ; k<K ; k++)
	{
		alpha[k] = ALPHA;
	}

	beta = new double[V];
	for(int v=0 ; v<V ; v++)
	{
		beta[v] = BETA;
	}

	numOfTopicWord = new double*[K];
	for(int k=0 ; k<K ; k++) 
	{
		numOfTopicWord[k] = new double[V];
		for(int v=0 ; v<V ; v++)
		{
			numOfTopicWord[k][v] = 0.0;
		}
	}

	totalNumOfTopic = new double[K];
	numOfTopicDoc = new double[K];
	for(int k=0 ; k<K ; k++) 
	{
		totalNumOfTopic[k] = 0.0;
		numOfTopicDoc[k] = 0.0;
	}

	q = new double*[M];
	for(int m=0 ; m<M ; m++) 
	{
		q[m] = new double[K];
		for(int k=0 ; k<K ; k++)
		{
			q[m][k] = 0.0;
		}

		int N = trainData->docs[m]->uniqueword;	
		int topic = (int)(rand() % K);

		q[m][topic] = 1.0;

		// initialize for z
		for(int n=0 ; n<N ; n++) 
		{			
			numOfTopicWord[topic][trainData->docs[m]->words[n]] += (double)trainData->docs[m]->nwords[n];
		} 

		totalNumOfTopic[topic] += (double)trainData->docs[m]->length;  
		numOfTopicDoc[topic] += 1.0;
	}

	theta = new double[K];
	for(int k=0 ; k<K ; k++)
	{
		theta[k] = 1.0/(double)K;
	}

	phi = new double*[K];
	for(int k=0 ; k<K ; k++)
	{
		phi[k] = new double[V];
	}  

	NormalizeMatrix(phi, numOfTopicWord, K, V);

	if(ManifoldRegularization)
	{
		tmpq = new double*[M];
		for(int m=0 ; m<M ; m++) 
		{
			tmpq[m] = new double[K];
			for(int k=0 ; k<K ; k++)
			{
				tmpq[m][k] = 0.0;
			}
		}

		tmpNumOfTopicWord = new double*[K];
		for(int k=0 ; k<K ; k++) 
		{
			tmpNumOfTopicWord[k] = new double[V];
			for(int v=0 ; v<V ; v++)
			{
				tmpNumOfTopicWord[k][v] = 0.0;
			}
		}

		tmpTotalNumOfTopic = new double[K];
		tmpNumOfTopicDoc = new double[K];
		for(int k=0 ; k<K ; k++) 
		{
			tmpTotalNumOfTopic[k] = 0.0;
			tmpNumOfTopicDoc[k] = 0.0;
		}

		if(GraphTF == true)
		{
			idxOfDocGraph = new int*[M];

			for(int m=0 ; m<M ; m++)
			{
				idxOfDocGraph[m] = new int[numOfNN];
			}

			string filename = fileBasis + "_dgTF.txt";

			if(ReadDocGraph(filename) == 0)
			{
				GraphConstructionTF();
			}
		}
	}

	return 0;
}

int LapDMM::ReadDocGraph(string filename)
{
	FILE *file;	

	file = fopen(filename.c_str(), "r");

	if(file == NULL)
	{
		return 0;
	}

	int numOfDoc = 0;
	int maxNumOfNN = 0;
	fscanf(file, "%d %d", &numOfDoc, &maxNumOfNN);

	if(numOfDoc != M)
	{
		return 0;
	}

	int idxOfDoc = 0;
	double simOfDoc = 0.0;
	for(int m=0 ; m<M ; m++)
	{
		for(int n=0 ; n<maxNumOfNN ; n++)
		{
			fscanf(file, "%d:%lf", &idxOfDoc, &simOfDoc);
			if(n < numOfNN)
			{
				idxOfDocGraph[m][n] = idxOfDoc;
			}
		}
	}

	fclose(file);

	return 1;
}

void LapDMM::Train()
{
	double tmpLik, lik, plik = 0.0;

	int start, elapsed;

	printf("Number of documents          = %d\n", M);
	printf("Number of words              = %d\n", V);
	printf("Number of topics             = %d\n", K);
	printf("Number of outer EM iteration = %d\n", maxIteration);

	start = clock();
	for(int t=0 ; t<maxIteration ; t++)
	{
		lik = 0.0;
		printf("iteration %d/%d..\n", t + 1, maxIteration);

		//compute the sum of the hyperparameters: alpha and beta
		Kalpha = vec_sum(K,alpha);
		Vbeta = vec_sum(V,beta);

		qEstimate();

		ComputeTheta(numOfTopicDoc);
		ComputePhi(numOfTopicWord,totalNumOfTopic);

		lik = LapDMM_lik(trainData);
		
		if(ManifoldRegularization)
		{
			lik += MR_lik(q);

			for(int tt=0 ; tt<mrIteration ; tt++)
			{
				UpdateTmpq();
				ComputeTheta(tmpNumOfTopicDoc);
				ComputePhi(tmpNumOfTopicWord,tmpTotalNumOfTopic);

				tmpLik = LapDMM_lik(trainData) + MR_lik(tmpq);

				if(tmpLik > lik)
				{
					lik = tmpLik;
					CopyTmpq();
				}

				else
				{
					break;
				}
			}
		}

		printf("likelihood = %g\t", lik); 

		if ((t>5)&&(fabs((lik - plik)/lik) < epsilon)) 
		{
			elapsed = clock() - start;
			printf("\nconverged. [%fs]\n", (double)elapsed/CLOCKS_PER_SEC);
			break;
		}

		plik = lik;		
	}

	ComputeTopicRepresentationSW();

	elapsed = clock() - start;
	printf("\nruning time is %f\n", (double)(elapsed/CLOCKS_PER_SEC));
}

void LapDMM::qEstimate() 
{	
	for(int m=0 ; m<M ; m++)
	{
		qUpdate(m);
	}
}

void LapDMM::qUpdate(int m) 
{
	int N = trainData->docs[m]->uniqueword;

	for(int k=0 ; k<K ; k++)
	{
		numOfTopicDoc[k] -= q[m][k];
		for(int n=0 ; n<N ; n++)
		{
			numOfTopicWord[k][trainData->docs[m]->words[n]] -= q[m][k] * trainData->docs[m]->nwords[n];
		}
		totalNumOfTopic[k] -= q[m][k] * trainData->docs[m]->length;
	}
	
	double sum = 0.0;
	for(int k=0 ; k<K ; k++) 
	{
		q[m][k] = 1.0;
		double count = 0.0;
		for(int n=0 ; n<N ; n++)
		{
			int v = trainData->docs[m]->words[n];
			for(int u=0 ; u<trainData->docs[m]->nwords[n] ; u++)
			{
				q[m][k] *= (numOfTopicWord[k][v] + u + beta[v]) / (totalNumOfTopic[k] + count + Vbeta);
				count++;
			}
		}

		q[m][k] *= numOfTopicDoc[k] + alpha[k];

		sum += q[m][k];
	}

	for(int k=0 ; k<K ; k++) 
	{
		q[m][k] = q[m][k] / sum;
	}

	for(int k=0 ; k<K ; k++)
	{
		numOfTopicDoc[k] += q[m][k];
		for(int n=0 ; n<N ; n++)
		{
			numOfTopicWord[k][trainData->docs[m]->words[n]] += q[m][k] * trainData->docs[m]->nwords[n];
		}
		totalNumOfTopic[k] += q[m][k] * trainData->docs[m]->length;
	}
}

void LapDMM::ComputeTheta(double* numOfTopicDoc) 
{
	for(int k=0 ; k<K ; k++) 
	{
		theta[k] = (numOfTopicDoc[k] + alpha[k]) / ((double)M + Kalpha);
	}
}

void LapDMM::ComputePhi(double** numOfTopicWord, double* totalNumOfTopic) 
{
	for(int k=0 ; k<K ; k++) 
	{
		for(int v=0 ; v<V ; v++) 
		{
			phi[k][v] = (numOfTopicWord[k][v] + beta[v]) / (totalNumOfTopic[k] + Vbeta);
		}
	}
}

double LapDMM::LapDMM_lik(dataset* corpus)
{
	double zz, lik;
	lik = 0;
	
	for(int m=0 ; m<corpus->M ; m++)
	{
		int N = corpus->docs[m]->uniqueword;

		for(int k=0 ; k<K ; k++)
		{
			for(int n=0 ; n<N ; n++) 
			{
				zz = 0.0;
				zz += phi[k][corpus->docs[m]->words[n]];

				if(zz != 0.0)
				{
					lik += q[m][k] * corpus->docs[m]->nwords[n] * log(zz);
				}
			}

			lik += q[m][k] * log(theta[k]);
		}
	}

	return lik;

}

void LapDMM::GraphConstructionTF()
{
	double sim = 0.0;
	double norm1 = 0.0;
	double norm2 = 0.0;
	idxOfDocGraph = new int*[M];
	double* tmpSimValue = new double[numOfNN];

	for(int m=0 ; m<M ; m++)
	{
		idxOfDocGraph[m] = new int[numOfNN];
		for(int k=0 ; k<numOfNN ; k++)
		{
			tmpSimValue[k] = -1000000.0;
		}

		int N1 = trainData->docs[m]->uniqueword;

		for(int mm=0 ; mm<M ; mm++)
		{
			if(m != mm)
			{
				sim = 0.0;
				norm1 = 0.0;
				norm2 = 0.0;

				int N2 = trainData->docs[mm]->uniqueword;
				for(int n=0 ; n<N1 ; n++)
				{
					int word = trainData->docs[m]->words[n];

					double value1 = (double)trainData->docs[m]->nwords[n] / (double)trainData->docs[m]->length;
					norm1 += value1 * value1;

					for(int nn=0 ; nn<N2 ; nn++)
					{
						if(word == trainData->docs[mm]->words[nn])
						{
							double value2 = (double)trainData->docs[mm]->nwords[nn] / (double)trainData->docs[mm]->length;
							sim += value1 * value2;
						}
					}
				}

				for(int n=0 ; n<N2 ; n++)
				{
					double value2 = (double)trainData->docs[mm]->nwords[n] / (double)trainData->docs[mm]->length;
					norm2 += value2 * value2;
				}

				norm1 = sqrt(norm1);
				norm2 = sqrt(norm2);
				sim = sim / (norm1 * norm2);

				if(sim > tmpSimValue[numOfNN-1])
				{
					for(int k=0 ; k<numOfNN ; k++)
					{
						if(sim > tmpSimValue[k])
						{
							for(int kk=numOfNN-1 ; kk>k ; kk--)
							{
								tmpSimValue[kk] = tmpSimValue[kk-1];
								idxOfDocGraph[m][kk] = idxOfDocGraph[m][kk-1];
							}
							tmpSimValue[k] = sim;
							idxOfDocGraph[m][k] = mm;
							break;
						}
					}
				}
			}
		}
/*
for(int k=0 ; k<numOfNN ; k++)
{
	printf("%d:%f\t",idxOfDocGraph[m][k], tmpSimValue[k]);
}
printf("\n");
*/
	}

	delete [] tmpSimValue;
}

void LapDMM::UpdateTmpq()
{
	for(int k=0 ; k<K ; k++) 
	{
		for(int v=0 ; v<V ; v++)
		{
			tmpNumOfTopicWord[k][v] = 0.0;
		}
	}

	for(int k=0 ; k<K ; k++) 
	{
		tmpTotalNumOfTopic[k] = 0.0;
		tmpNumOfTopicDoc[k] = 0.0;
	}

	for(int m=0 ; m<M ; m++)
	{
		for(int k=0 ; k<K ; k++)
		{
			tmpq[m][k] = (1.0 - gamma) * q[m][k];

			for(int kk=0 ; kk<numOfNN ; kk++)
			{
				tmpq[m][k] += gamma / (double)numOfNN * q[idxOfDocGraph[m][kk]][k];
			}
		}

		int N = trainData->docs[m]->uniqueword;

		for(int k=0 ; k<K ; k++)
		{
			tmpNumOfTopicDoc[k] += tmpq[m][k] * trainData->docs[m]->length;
			for(int n=0 ; n<N ; n++)
			{
				tmpNumOfTopicWord[k][trainData->docs[m]->words[n]] += tmpq[m][k] * trainData->docs[m]->nwords[n];
			}
			tmpTotalNumOfTopic[k] += tmpq[m][k] * trainData->docs[m]->length;
		}
	}
}

void LapDMM::CopyTmpq()
{
	for(int k=0 ; k<K ; k++) 
	{
		for(int v=0 ; v<V ; v++)
		{
			numOfTopicWord[k][v] = 0.0;
		}
	}

	for(int k=0 ; k<K ; k++) 
	{
		totalNumOfTopic[k] = 0.0;
		numOfTopicDoc[k] = 0.0;
	}

	for(int m=0 ; m<M ; m++)
	{
		for(int k=0 ; k<K ; k++)
		{
			q[m][k] = tmpq[m][k];
		}

		int N = trainData->docs[m]->uniqueword;

		for(int k=0 ; k<K ; k++)
		{
			numOfTopicDoc[k] += q[m][k];
			for(int n=0 ; n<N ; n++)
			{
				numOfTopicWord[k][trainData->docs[m]->words[n]] += q[m][k] * trainData->docs[m]->nwords[n];
			}
			totalNumOfTopic[k] += q[m][k] * trainData->docs[m]->length;
		}
	}
}

double LapDMM::MR_lik(double** q)
{
	double lik = 0.0;

	for(int m=0 ; m<M ; m++)
	{
		for(int k=0 ; k<K ; k++)
		{
			for(int kk=0 ; kk<numOfNN ; kk++)
			{
				lik += pow((q[m][k] - q[idxOfDocGraph[m][kk]][k]),2.0);
			}
		}
	}

	return -lambda * (double)M * lik;
}

void LapDMM::ComputeTopicRepresentationSW()
{
	double sum = 0.0;

	for(int v=0 ; v<V ; v++)
	{
		sum = 0.0;
		for(int k=0 ; k<K ; k++)
		{
			sum += theta[k] * phi[k][v];
		}
		for(int k=0 ; k<K ; k++)
		{
			phi[k][v] = theta[k] * phi[k][v] / sum;
		}
	}

	for(int m=0 ; m<M ; m++)
	{
		sum = 0.0;
		int N = trainData->docs[m]->uniqueword;
		for(int k=0 ; k<K ; k++)
		{
			q[m][k] = 0.0;
			for(int n=0 ; n<N ; n++)
			{
				q[m][k] += phi[k][trainData->docs[m]->words[n]] * trainData->docs[m]->nwords[n];
			}
			sum += q[m][k];
		}
		for(int k=0 ; k<K ; k++)
		{
//			q[m][k] = q[m][k] / sum;
		}
	}
}

void LapDMM::NormalizeMatrix(double **dst, double **src ,int R, int L)
{
	/* column-wise normalize from src -> dst */
	double z = 0.0;
	int i, j;

	for(int k=0 ; k<R ; k++) 
	{
		z = 0.0;
		for(int v=0 ; v<L ; v++)
		{
			z += src[k][v];
		}

		for(int v=0 ; v<L ; v++)
		{
			dst[k][v] = src[k][v] / z;
		}
	}
}

void LapDMM::Evaluation()
{
	int* z = new int[M];
	for(int m=0 ; m<M ; m++)
	{
		double max = -10000000.0;
		for(int k=0 ; k<K ; k++)
		{
			if(q[m][k] > max)
			{
				max = q[m][k];
				z[m] = k;
			}
		}
	}

	double clusterNMI = evaluation->ClusteringNMI(trainData,z,K);
	printf("NMI of clustering = %g\n",clusterNMI);

	double clusterACC = evaluation->ClusteringACC(trainData,z,K);
	printf("ACC of clustering = %g\n",clusterACC);

	delete [] z;

	return;
}