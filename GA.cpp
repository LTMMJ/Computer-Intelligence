/************************************************************
*      ga.cpp - GA Program for CSCI964 - Ass2
*      Written by: Koren Ward May 2010
*      Modified by: <Ting Li & Student Number: 6286987 2019.5.20>
*      Changes: 1. Define a structure CityINFO to store cities' information.
					2. Define a structure TSPPara and put all relative parameters in, which makes it more convenient to call these parameteres in functions than using arrays.
					3. Change cNumGen to be 3800(every 200th generation is displayed), and the PopSize is 300. Add MaxValue(value is 1000000000) and FIX(value is 100000) to implement some functions.
					4. Define CityDist and CityCost array to store result.
					5. Define GetData() function to read data from file. Define a roulette wheel selection function RWselection().
					6. Change InitPop(), EvaluateFitness(), Crossover(), and Mutate() a lot, they call TSPPara City instead of arrays. And some details have also changed. Make a little change in Change() function.
*************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <bits/stdc++.h>
#include <stdlib.h>
using namespace std;
#define MaxValue 1000000000 //Define the max value of cost
typedef struct CityINFO	{//to store city's coordinates (x, y), and target t.
	double x;
	double y;
	int t;
}CityINFO;
CityINFO CityInfo[1000];

//======================PARAMETERS=========================
const double cCrossoverRate = 0.75;//0.75
const double cMutationRate = 0.001;//0.001
const int    cNumGens = 3800;
const int    PopSize = 300; // must be an even number
const int    Seed = 1234;
const int cIndividualLength = 100;// CityNum
int CityNum;
double CityDist[cIndividualLength][cIndividualLength];
double CityCost[cIndividualLength][cIndividualLength];//This is the lookup table of the cost of traveling between cities
const int FIX = 100000;
double SHUF[cIndividualLength + 1];

typedef struct TSPPara {//relative parameters
	int pop[PopSize][cIndividualLength + 1];
	double Fitness[PopSize];
	double Cost[PopSize];
	int BestRoute[cIndividualLength + 1];
	double BestFitness;
	double BestValue;
	int BestNum;
}TSPPara;

//=========================FUNCTIONS======================
void GetData();//Get data from file
void InitPop(TSPPara &City);//Load the initial pop members with random cities' serial numbers
void RWselection(TSPPara &City);//Roulette wheel selection
void EvaluateFitness(TSPPara &City);
void Crossover(TSPPara &City, double CR);
void Copy(int *P1, int *P2);
double Getfitness(int a[cIndividualLength + 1]);
void Mutate(TSPPara & City, double MR);
void CalculateCost();

//=========================MAIN==================================
int main(int argc, char *argv[]) {
	TSPPara City;
	srand(Seed);
	GetData();
	CalculateCost();
	InitPop(City);
	EvaluateFitness(City);
	for (int i = 0; i < cNumGens; i++) {
		RWselection(City);
		Crossover(City, cCrossoverRate);
		Mutate(City, cMutationRate);
		EvaluateFitness(City);
		//Output(i, City);
		if ((i == 0) || ((i + 1) % 200 == 0)) {
			cout << "Epoch: " << i + 1 << "\t The best fitness: " << City.BestFitness << "\t The minimum tour cost: " << City.BestValue << endl;
		}
	}
	getchar();
	return 0;
}
//===========================================================
void GetData() {		//Get data from file.
	fstream fin("tsp100.txt");
	if (fin) {
		int rowNum = 0, c;
		double num, a, b;
		fin >> num;
		CityNum = num;
		while (fin >> a >> b >> c) {
			CityInfo[rowNum++] = { a, b, c };
		}
		fin.close();
	}
	else {
		cerr << "There is no file found£¡" << endl;
	}
}

bool check(TSPPara &City, int pop, int num, int k)
{
	int i;
	for (i = 0; i <= num; i++) {
		if (k == City.pop[pop][i])
			return true;//The newly generated node exists in the path that was already generated.
	}
	return false;//Dose not exist.
}

void InitPop(TSPPara &City) {
	int i, j, r;
	srand(Seed);
	for (i = 0; i < PopSize; i++) {	//Define original value
		City.pop[i][0] = 0;
		City.pop[i][cIndividualLength] = 0;
		City.BestValue = MaxValue;
		City.BestFitness = 0;
	}
	for (i = 0; i < PopSize; i++){
		for (j = 1; j < cIndividualLength; j++){
			r = rand() % (cIndividualLength - 1) + 1;//Produce a random from 1~²úÉú1¡«cIndividualLength-1
			while (check(City, i, j, r)){//Produce city serial number randomly.
				r = rand() % (cIndividualLength - 1) + 1;
			}
			City.pop[i][j] = r;
		}
	}
}

void RWselection(TSPPara &City) {//Roulette wheel selection
	int i, j, k;
	int tpop[PopSize][cIndividualLength + 1];
	double s, sum = 0;
	double Assi[PopSize], SelectP[PopSize + 1];
	for (i = 0; i < PopSize; i++) {
		sum += City.Fitness[i];
	}
	for (i = 0; i < PopSize; i++) {
		Assi[i] = City.Fitness[i] / sum;
	}
	SelectP[0] = 0;
	for (i = 0; i < PopSize; i++) {
		SelectP[i + 1] = SelectP[i] + Assi[i] * RAND_MAX;
	}
	memcpy(tpop[0], City.pop[City.BestNum], sizeof(tpop[0]));//To copy data from City.BestNum to tpop[0].
	for (k = 1; k < PopSize; k++) {
		double ran = rand() % RAND_MAX + 1;
		s = (double)ran / 100.0;
		for (i = 1; i < PopSize; i++) {
			if (SelectP[i] >= s) { break; }
		}
		memcpy(tpop[k], City.pop[i - 1], sizeof(tpop[k]));
	}
	for (i = 0; i < PopSize; i++) {
		memcpy(City.pop[i], tpop[i], sizeof(tpop[i]));
	}
}

void EvaluateFitness(TSPPara &City) {//Evaluates fitness
	int i, j, s, e, best = 0;
	for (i = 0; i < PopSize; i++) {
		City.Cost[i] = 0;
		for (j = 1; j <= cIndividualLength - 1; j++) {
			s = City.pop[i][j - 1];
			e = City.pop[i][j];
			City.Cost[i] = City.Cost[i] + CityCost[s][e];	//total cost
		}
		City.Fitness[i] = FIX / City.Cost[i];
		if (City.Fitness[i] > City.Fitness[best]) {	//choose the biggest fitness value
			best = i;
		}
	}
	Copy(City.BestRoute, City.pop[best]); //copy the best one to City.BestRoute
	City.BestFitness = City.Fitness[best];
	City.BestValue = City.Cost[best];
	City.BestNum = best;
}

void Crossover(TSPPara &City, double CR)//CR is crossover rate 0.75
{
	int i, j, k, l, m, n, cm, cn;
	int Temp1[cIndividualLength + 1];
	for (i = 0; i < PopSize; i++) {
		double s = ((double)(rand() % RAND_MAX)) / RAND_MAX;
		if (s < CR) {
			cn = rand() % PopSize;
			cm = cn;
			if (cm == City.BestNum || cn == City.BestNum) { continue; }//If the optimal is encountered, the next cycle is directly carried out
			l = rand() % (cIndividualLength / 2) + 1;  //1~first half part
			m = rand() % (cIndividualLength - l) + 1; //1~cIndividualLength
			memset(SHUF, 0, sizeof(SHUF));//replace SHUF's current sizeof(SHUF) bits with 0.
			Temp1[0] = Temp1[cIndividualLength] = 0;
			for (j = 1; j <= l; j++) {//Shuffling order(is random), and the selected one is marked as 1.
				Temp1[j] = City.pop[cn][m + j - 1];
				SHUF[Temp1[j]] = 1;
			}
			for (k = 1; k < cIndividualLength; k++) {
				if (SHUF[City.pop[cm][k]] == 0) {
					Temp1[j++] = City.pop[cm][k];
					SHUF[City.pop[cm][k]] = 1;
				}
			}
			memcpy(City.pop[cm], Temp1, sizeof(Temp1));
		}
	}
}

double Getfitness(int a[cIndividualLength + 1]) {
	int i, s, e;
	double COST = 0.0;
	for (i = 0; i < cIndividualLength - 1; i++) {
		s = a[i];
		e = a[i + 1];
		COST += CityCost[s][e];
	}
	return FIX / COST;
}

void Mutate(TSPPara & City, double MR) {
	int i, n, m;
	int t[cIndividualLength + 1];
	for (i = 0; i < PopSize; i++) {
		double a = ((double)(rand() % RAND_MAX)) / RAND_MAX;
		n = rand() % PopSize;
		if (a < MR && n != City.BestNum) {// To ensure the best one do not mutate.
			int b, c, d;
			b = (rand() % (cIndividualLength - 1)) + 1;
			c = (rand() % (cIndividualLength - 1)) + 1;
			Copy(t, City.pop[i]);
			if (b > c) {//To ensure c >= b.
				d = b; b = c; c = d;
			}
			for (m = b; m < (b + c) / 2; m++) {
				d = t[m]; t[m] = t[b + c - m]; t[b + c - m] = d;
			}
			if (Getfitness(t) < Getfitness(City.pop[i])) {
				b = (rand() % (cIndividualLength - 1)) + 1;
				c = (rand() % (cIndividualLength - 1)) + 1;
				memcpy(t, City.pop[i], sizeof(t));
				if (b > c) {
					d = b; b = c; c = d;
				}
				for (m = b; m < (b + c) / 2; m++) {
					d = t[m]; t[m] = t[b + c - m]; t[b + c - m] = d;
				}
				if (Getfitness(t) < Getfitness(City.pop[i])) {
					b = (rand() % (cIndividualLength - 1)) + 1;
					c = (rand() % (cIndividualLength - 1)) + 1;
					memcpy(t, City.pop[i], sizeof(t));
					if (b > c) {
						d = b; b = c; c = d;
					}
					for (m = b; m < (b + c) / 2; m++) {
						d = t[m]; t[m] = t[b + c - m]; t[b + c - m] = d;
					}
				}
			}
			memcpy(City.pop[i], t, sizeof(t));
		}
	}
}

void CalculateCost(){
	int i, j;
	double temp1, temp2;
	for (i = 0; i < cIndividualLength; i++) {
		for (j = 0; j <= cIndividualLength; j++) {//The last city should be able to return to the departure node.
			temp1 = CityInfo[j].x - CityInfo[i].x; temp2 = CityInfo[j].y - CityInfo[i].y;
			CityDist[i][j] = sqrt(temp1 * temp1 + temp2 * temp2);
			int flag = CityInfo[j].t * CityInfo[i].t;
			switch (flag) {
			case 1: CityCost[i][j] = CityDist[i][j] * 10.0; break;
			case 2: CityCost[i][j] = CityDist[i][j] * 7.5; break;
			case 3: CityCost[i][j] = CityDist[i][j] * 5.0; break;
			case 4: CityCost[i][j] = CityDist[i][j] * 5.0; break;
			case 6: CityCost[i][j] = CityDist[i][j] * 2.5; break;
			case 9: CityCost[i][j] = CityDist[i][j] * 1.0; break;
			}
		}
	}
}

void Copy(int *P1, int *P2) {
	for (int i = 0; i < cIndividualLength; i++) {
		P1[i] = P2[i];
	}
}