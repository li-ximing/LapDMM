/*
* Copyright (C) 2019 by
* 
* 	Ximing Li
*	liximing86@gmail.com
* 	College of Computer Science and Technology, JiLin University
*
*/

#include "dataset.h"

/*
// read corpus
int dataset::ReadCorpus(char* filename)
{	
	int uniqueword, uniquelabel, label, count, word;

	printf("reading data from %s\n", filename);

	FILE *file;	

	file = fopen(filename, "r");

	int nd = 0; 
	int nw = 0;
	while((fscanf(file, "%10d", &uniqueword) != EOF))
	{
		fscanf(file, "%d", &label);

		for (int n = 0; n < uniqueword; n++) 
		{
			fscanf(file, "%10d:%10d", &word, &count);
//			word = word - OFFSET;
			if (word >= nw) 
			{ 
				nw = word + 1; 
			}
		}
		nd++;
	}

	M = nd;
	V = nw;
	docs = new document*[M];
	length = new int[M];

	rewind(file);

	for(int nd=0 ; nd<M ; nd++)
	{
		length[nd] = 0;
		fscanf(file, "%10d", &uniqueword);
		fscanf(file, "%10d", &label);

		for (int n = 0; n < uniqueword; n++) 
		{
			fscanf(file, "%10d:%10d", &word, &count);
			length[nd] += count; 
		}

		document* pdoc = new document(length[nd]);	

		add_doc(pdoc, nd);
	}

	rewind(file);

	for(int nd=0 ; nd<M ; nd++)
	{
		fscanf(file, "%10d", &uniqueword);
		fscanf(file, "%10d", &label);

		this->docs[nd]->label = label;

		int temp = 0;
		for (int n = 0; n < uniqueword; n++) 
		{
			fscanf(file, "%10d:%10d", &word, &count);
		
			for(int nn=0 ; nn<count ; nn++)
			{
				this->docs[nd]->words[temp] = word;
				temp++;
			}
		}
	}

	delete [] length;

	fclose(file);

//	printf("number of docs    : %d\n", M);
//	printf("number of terms   : %d\n", V);

	return 0;
}
*/

// read corpus
int dataset::ReadCorpus(string filename)
{	
	int uniqueword, uniquelabel, label, count, word;
	int ncount = 0;

	printf("reading data from %s\n", filename.c_str());

	FILE *file;	

	file = fopen(filename.c_str(), "r");

	int nd = 0; 
	int nw = 0;
	while((fscanf(file, "%10d", &uniqueword) != EOF))
	{
		fscanf(file, "%10d", &label);

		for (int n = 0; n < uniqueword; n++) 
		{
			fscanf(file, "%10d:%10d", &word, &count);
			//			word = word - OFFSET;
			if (word > nw) 
			{ 
				nw = word; 
			}
		}

		nd++;		
	}

	M = nd;
	V = nw + 1;
	docs = new document*[M];

	rewind(file);

	uniquelabel = 0;
	M = 0;

	for(int m=0 ; m<nd ; m++)
	{
		fscanf(file, "%10d", &uniqueword);
		fscanf(file, "%10d", &label);

		document* pdoc = new document(uniqueword);

		if(label+1 > ulabel)
		{
			ulabel = label+1;
		}

		pdoc->label = label;

		for (int n = 0; n < uniqueword; n++) 
		{
			fscanf(file, "%10d:%10d", &word, &count);
			pdoc->words[n] = word;
			pdoc->nwords[n] = count;
			pdoc->length += count;
		}

		ncount += pdoc->length;

		add_doc(pdoc, M);
		M++;
	}

	fclose(file);

/*
printf("number of word tokens    : %d\n", ncount);
printf("averge document length   : %f\n", (double)ncount/(double)M);
*/

/*
FILE * fout = fopen("result", "w");
for(int m=0 ; m<M ; m++)
{
	if(this->docs[m]->uniqueword != 1)
	{
		fprintf(fout, "%d %d ", this->docs[m]->uniqueword, this->docs[m]->label);
		for (int n = 0; n < this->docs[m]->uniqueword ; n++) 
		{
			fprintf(fout, "%d:%d ", this->docs[m]->words[n], this->docs[m]->nwords[n]);
		}
		fprintf(fout, "\n");
	}
}
*/
//	printf("number of docs    : %d\n", M);
//	printf("number of terms   : %d\n", V);

	return 0;
}


// initialization 
int dataset::InitDataset(string filename)
{
	if(ReadCorpus(filename) != 0)
	{
		printf("read dataset error");
		return 1;
	}

	return 0;
}
