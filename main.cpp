/*
* Copyright (C) 2019 by
* 
* 	Ximing Li
*	liximing86@gmail.com
* 	College of Computer Science and Technology, JiLin University
*
*/

#include "LapDMM.h"
#include "dataset.h"

int main()
{
	string fileOfTrainData = "./dataset/Trec.txt";

	LapDMM TM;
	TM.InitEst(fileOfTrainData);	
	TM.Train();
	TM.Evaluation();

   	return 0;
}