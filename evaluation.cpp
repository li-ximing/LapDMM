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
#include "evaluation.h"

EvaluationTopic::EvaluationTopic(string fileBasis)
{
	topword = 10;

	filename = fileBasis + "_voc.txt";
}

double EvaluationTopic::ClusteringNMI(dataset* trainData, int* z, int K)
{
	double NMI = 0.0;
	double entropyOfCluster = 0.0;
	double entropyOfLabel = 0.0;

	int numOfCluster = K;
	double numOfDoc = (double)trainData->M;
	int numOfLabel = trainData->ulabel;
	
	int* sizeOfCluster = new int[numOfCluster];
	int* sizeOfLabel = new int[numOfLabel];
	int** mat = new int*[numOfCluster];
	for(int c=0 ; c<numOfCluster ; c++)
	{
		sizeOfCluster[c] = 0;
		mat[c] = new int[numOfLabel];
		for(int l=0 ; l<numOfLabel ; l++)
		{
			mat[c][l] = 0;
		}
	}
	for(int l=0 ; l<numOfLabel ; l++)
	{
		sizeOfLabel[l] = 0;
	}

	for(int m=0 ; m<numOfDoc ; m++)
	{
		sizeOfCluster[z[m]]++;
		sizeOfLabel[trainData->docs[m]->label]++;
		mat[z[m]][trainData->docs[m]->label]++;
	}

	for(int c=0 ; c<numOfCluster ; c++)
	{
		entropyOfCluster -= (double)sizeOfCluster[c]/numOfDoc * log10((double)sizeOfCluster[c]/numOfDoc+1e-6)
							/ log10(2.0);
		
		for(int l=0 ; l<numOfLabel ; l++)
		{
			if(mat[c][l] != 0)
			{
				NMI += (double)mat[c][l] / numOfDoc *
					log10(numOfDoc*(double)mat[c][l]/((double)sizeOfCluster[c]*(double)sizeOfLabel[l]) + 1e-6)
					/ log10(2.0);
			}
		}
	}

	for(int l=0 ; l<numOfLabel ; l++)
	{
		entropyOfLabel -= (double)sizeOfLabel[l]/numOfDoc * log10((double)sizeOfLabel[l]/numOfDoc + 1e-6)
						/ log10(2.0);
	}

	delete [] sizeOfCluster;
	delete [] sizeOfLabel;
	for(int c=0 ; c<numOfCluster ; c++)
	{
		delete [] mat[c];
	}
	delete [] mat;

	return NMI / sqrt(entropyOfCluster*entropyOfLabel);
}

double EvaluationTopic::ClusteringACC(dataset* trainData, int* z, int K)
{
	double ACC = 0.0;

	int numOfCluster = K;
	int numOfDoc = trainData->M;
	int numOfLabel = trainData->ulabel;

	if(numOfCluster == numOfLabel)
	{
		int* map = new int[numOfCluster];
		int* mat = new int[numOfCluster*numOfLabel];
		for(int c=0 ; c<numOfCluster ; c++)
		{
			for(int l=0 ; l<numOfLabel ; l++)
			{
				mat[c*numOfCluster+l] = 0;
			}
		}

		for(int m=0 ; m<numOfDoc ; m++)
		{
			mat[z[m]*numOfCluster+trainData->docs[m]->label]--;
		}

		hungarian_t prob;

		hungarian_init(&prob,mat,numOfCluster,numOfLabel,HUNGARIAN_MIN);
		hungarian_solve(&prob);

		for(int c=0 ; c<numOfCluster ; c++)
		{
			for(int l=0 ; l<numOfLabel ; l++)
			{
				if(prob.a[c] == l)
				{
					map[c] = l;
					break;
				}
			}
		}

		for(int m=0 ; m<numOfDoc ; m++)
		{
			if(map[z[m]] == trainData->docs[m]->label)
			{
				ACC++;
			}
		}

		hungarian_fini(&prob);
		delete [] map;
		delete [] mat;
	}
	
	return ACC / (double)numOfDoc;
}