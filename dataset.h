/*
* Copyright (C) 2015 by
* 
* 	Ximing Li
*	liximing86@gmail.com
* 	College of Computer Science and Technology, JiLin University
*
*/

#ifndef	_DATASET_H
#define	_DATASET_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

using namespace std;
class synonym
{
public:
	int maxNumOfSynonym;				// the max number of synonyms saved 
	int* numOfSynonym;					// the number of synonyms saved for each word type
	int** idxOfSynonym;					// the indexes of synonyms
	double** simOfSynonym;				// the similarity scores of synonyms

	synonym(string filename, int V, double threshold)
	{
		maxNumOfSynonym = 20;
		FILE *file;	
		file = fopen(filename.c_str(), "r");

		numOfSynonym = new int[V];
		idxOfSynonym = new int*[V];
		simOfSynonym = new double*[V];

		int numOfSynonymFile = 0;
		int word = 0;
		double sim = 0.0;
		int* tmpIdxOfSynonym = new int[maxNumOfSynonym];
		double* tmpSimOfSynonym = new double[maxNumOfSynonym];

		for(int v=0 ; v<V ; v++)
		{
			fscanf(file, "%10d", &numOfSynonymFile);
		
			if(numOfSynonymFile == 0)
			{
				numOfSynonym[v] = 0;
				continue;
			}
		
			int numOfSynonymSaved = 0;
			for(int n=0 ; n<numOfSynonymFile ; n++) 
			{
				fscanf(file, "%10d:%20lf", &word, &sim);
				if(sim >= threshold)
				{
					tmpIdxOfSynonym[numOfSynonymSaved] = word;
					tmpSimOfSynonym[numOfSynonymSaved] = sim;
					numOfSynonymSaved++;
				}
			}
			
			numOfSynonym[v] = numOfSynonymSaved;
			if(numOfSynonymSaved != 0)
			{
				idxOfSynonym[v] = new int[numOfSynonymSaved];
				simOfSynonym[v] = new double[numOfSynonymSaved];
			}

			for(int n=0 ; n<numOfSynonymSaved ; n++)
			{
				idxOfSynonym[v][n] = tmpIdxOfSynonym[n];
				simOfSynonym[v][n] = tmpSimOfSynonym[n];
			}
		}

		fclose(file);

		delete [] tmpIdxOfSynonym;
		delete [] tmpSimOfSynonym;
	}
};

class document 
{
public:
	int length;
	int uniqueword;						// the number of unique word in the document 

	int uniquelabel;
	int label;

	int* words;							// the index of each word occured in the document: size uniqueword 
	int* nwords;						// the number of each word occured in the document: size uniqueword 

	document() 
	{
		this->words = NULL;
		this->nwords = NULL;
		this->uniqueword = 0;	
		this->length = 0;
	}

	document(int uniqueword) 
	{
		this->length = 0;
		this->uniqueword = uniqueword;
		this->words = new int[uniqueword];	
		this->nwords = new int[uniqueword];
	}

	document(int uniqueword, int* words) 
	{
		this->uniqueword = uniqueword;
		this->words = new int[uniqueword];
		for(int i=0 ; i<uniqueword ; i++) 
		{
			this->words[i] = words[i];
		}
	}

	document(int uniqueword, int* words, int* nwords) 
	{
		this->uniqueword = uniqueword;
		this->words = new int[uniqueword];
		for(int i=0 ; i<uniqueword ; i++) 
		{
			this->words[i] = words[i];
		}
		this->nwords = new int[uniqueword];
		for(int i=0 ; i<uniqueword ; i++) 
		{
			this->nwords[i] = nwords[i];
		}
	}

	~document() 
	{
		if (words) 
		{
			delete [] words;
		}

		if (nwords) 
		{
			delete [] nwords;
		}
	}
};

class dataset 
{
public:
	document** docs;				// the corpus

	int M;							// the number of documents in the corpus
	int V;							// the number of unique words in the corpus

	int ulabel;	//number of labels

	int* labels;	// number for each label
	int sumlabel;	// the total number of label

	// the variables of tf-idf
	int* length;					// the length of each documents in the corpus: size M
	int* docword;					// the number of documents contains the word: size V
	double** tfidf;					// the tf-idf matrix of the corpus: size M * doc.uniqueword

	dataset() 
	{
		docs = NULL;
		length = NULL;
		docword = NULL;
		tfidf = NULL;
		M = 0;
		V = 0;
	}

	dataset(int M) 
	{
		this->M = M;
		this->V = 0;
		docs = new document*[M];	
		length = new int [M];
	}   

	~dataset() 
	{
		if(docs) 
		{
			for(int i=0 ; i<M ; i++) 
			{
				delete [] docs[i];
			}
		}
		delete [] docs;

		if(tfidf) 
		{
			for(int i=0 ; i<M ; i++) 
			{
				delete [] tfidf[i];
			}
		}

		if (length) 
		{
			delete [] length;
		}

		if (docword) 
		{
			delete [] docword;
		}
		
	}

	void deallocate() 
	{
		if (docs) 
		{
			for(int i=0 ; i<M ; i++) 
			{
				delete [] docs[i];
			}
		}
		delete [] docs;
		docs = NULL;
	}

	void add_doc(document* doc, int idx) 
	{
//		if (0<=idx && idx<M) 
		{
			docs[idx] = doc;
		}
	}   

	// read corpus 
	int ReadCorpus(string filename);
	
	// initialization
	int InitDataset(string filename);
};

#endif

