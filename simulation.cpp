// Term 4 tutorial 2014 - Landscape stickyness
#include <cstdlib>   // stdlib for lots of functions     
#include <cmath>     
#include <iostream>  // input-output control
//#include <ctime>     // to initialise random number generator
#include <fstream>   // fstream for file output
#include <conio.h>   // I use for getch() to pause at the end
#include <random>    // built in C++ random numbers
#include <Windows.h> // might be needed for some functions

using namespace std;
#pragma hdrstop

#define SEED GetTickCount()  // for the random number generator, time(0) for random
#define NINDS 500000      // set the max population size 
#define START_POP 40000   // how many individuals at the start
#define FILE1 "T4_2014_vary_stickyness_Time5000_K10_M20_NegExp_indsa.csv"         // output mean values for every locus
#define FILE2 "T4_2014_vary_stickyness_Time5000_K10_M20_NegExpa.csv"         // output mean values for every locus
#define SX 50            // size of the grid (width), usually 50
#define SY 50			// size of the grid (height)
#define MUT 0.01          // mutation rate
#define FREQUENCY 100   // how often to collect the data, usually 100
#define PI2 6.28318
#define PICK_IND (int)(R1(generator)*top)   // this selects a random individual

class guy {
public:
	float age;
	float x;
	float y;
	float strategy;
}; 

class cell {
public:
	int pop;
	// will have 'stickiness' values
	// might need to keep a list of possible mates
};

// global definitions
bool header2_done = false;   // this is just for the file output, have we printed a header yet?
bool header1_done = false;   // bools just have true/false or 1/0 values
int year;
int movedist = 0;
int K = 20;   // local equilibrium density
int top = START_POP;
int rep = 0;
float mortality = (float)0.2;  // this sort of casting just prevents warnings, it will run fine without the (float) cast
float stickyness = (float)10;
// float dispersal = (float)-0.15;   // put in as minus as otherwise will be multiplied by -1 to get the distance, it's the mean of the neg exp

// set up the big arrays
guy inds[NINDS];   // make an array of individuals called 'inds'
cell cells[SX + 1][SY + 1];  // array of grid cells containing population sizes

// function prototypes
void runmodel(void);
void output(void);
void initialise(void);

// functions for random number generation	
default_random_engine generator(SEED);   // reset the random number generator using the timer
uniform_int_distribution<int> pick_cellX(0, SX - 1);     // pick_cell(0,(SX-1)); to do the whole grid, used for initialise and global movement
bernoulli_distribution coinflip(0.5);                 // 0 or 1 with given probability
uniform_real_distribution<float> R1(0.0, 1.0);
uniform_real_distribution<double> RandFloat(0.0000001, 1.0);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
int main()
{
	cout << "\nAutumn term 2014 - tutorial modelling project v1\n";

	for (rep = 1; rep <= 3; ++rep) {
		for (stickyness = 1; stickyness <= 20; stickyness += 0.5) {
			initialise();
			for (year = 1; year < 5002; ++year) {
				// print to screen to give an idea of progress
				if (!(year % 500)) cout << "\n Rep:" << rep << " of 3.   Sticky:" << stickyness << " Time:" << year << "  Top:" << top;
				runmodel();
				if (!(year%FREQUENCY)) output();
			}
		}
	}	//
	cout << "\n\n\tAny key to exit.";
	_getch();   // pauses
	return 0;
}  //


void runmodel(void) {
	int clicks = 0;
	int i;
	double dist, angle, dx, dy, dnx, dny;  // these are used for the dispersal bit

	while (clicks<(top * 2)) {
		// find an individual
		i = PICK_IND;   // this is going to be the mother   
		inds[i].age += R1(generator);  // add to the age, mean of 0.5, as you're going to get 2 events on average a year
		if (R1(generator)<0.5) {   // birth event
			if (R1(generator)>(float)cells[(int)inds[i].x][(int)inds[i].y].pop / (float)(K * 2) && R1(generator)<(1 - inds[i].strategy)) {   // density dependent and affected by the 'strategy'

				// there's a new individual
				inds[top].age = 0;
				// we could just pick a random location with distance in x and y based on strategy close to the parent.... this is what was done for the first simulations
				//	if(R1(generator)<0.5) inds[top].x = inds[i].x + (inds[i].strategy * stickyness * R1(generator));   // just picking a random position on the grid
				//	else inds[top].x = inds[i].x - (inds[i].strategy * stickyness * R1(generator));
				//	if(R1(generator)<0.5) inds[top].y = inds[i].y + (inds[i].strategy * stickyness * R1(generator));
				//	else inds[top].y = inds[i].y - (inds[i].strategy * stickyness * R1(generator));

				// negative exponential dispersal distance from parent based on strategy
				dist = (inds[i].strategy * stickyness * -1)*(log(RandFloat(generator)));	// Neg exp distribution, RandFloat starts just above zero avoid log zero 
				if (dist>30) dist = 30;  // I think very long distance moves crashed the program
				angle = R1(generator)*PI2;	// a random direction in radians
				dx = dist*(cos(angle));	// polar to rectangular
				dy = dist*(sin(angle));

				dnx = dx + inds[i].x;   // add to mother's location
				dny = dy + inds[i].y;

				inds[top].x = (float)dnx;   // then locate the new guy
				inds[top].y = (float)dny;

				// Arena is a torus
				if (inds[top].x>SX) inds[top].x -= SX;
				if (inds[top].y>SY) inds[top].y -= SY;
				if (inds[top].x<0) inds[top].x += SX;
				if (inds[top].y<0) inds[top].y += SY;

				// mutation
				if (R1(generator)<MUT) {
					inds[top].strategy = inds[i].strategy + (R1(generator) / 20) - (R1(generator) / 20);	 // +/- parent's value, max +0.05/-0.05, with a humped distribution around 0	  
					if (inds[top].strategy<0) inds[top].strategy = 0;    // mutations are capped at 0
				}
				else
				{
					inds[top].strategy = inds[i].strategy;  // you just get your mum's strategy exactly
				}

					//	if(!off) {   // only actually implement all this if it's still on the grid, only need to check under some circumstances!	
				++cells[(int)inds[top].x][(int)inds[top].y].pop;    // add one to the local population size
				++top;
			}
		}

		else  // death event
		{
			if (R1(generator)<mortality) {   // this isn't affected by density there's just a constant mortality rate

				cells[(int)inds[i].x][(int)inds[i].y].pop--;  // reduce local population size	
				inds[i] = inds[top - 1];  // move the top individual in the pile to the space just vacated (easy way to avoid gaps in the pile)

				--top;  // make the total population size one smaller
			}
		}

		++clicks;  // this just counts up the number of events this timestep

	}  // end of the 'year'

}  // end of runmodel() function


void output(void) {

	double m_strat = 0.0;
	double m_age = 0.0;

	for (int i = 0; i<top; i++) {  // loop through all individuals
		m_strat += inds[i].strategy;   // add up strategies etc.
		m_age += inds[i].age;  // this is the current age, age at death might be more interesting?
	}

	m_age = m_age / (top - 1);   // divide by pop size to give mean
	m_strat = m_strat / (top - 1);

	// could do variance here too

	ofstream out2;    // output the overall mean values
	out2.open(FILE2, ios::app);  //  app appends to a file, ios::trunc is an option
	out2.precision(0);  // number of decimal places

	if (!header2_done) {   // only if header not already written
		out2 << "sticky,rep,time,mutation,mortality,K,pop,mean_strat,mean_age\n";
		header2_done = true;
	}

	out2 << stickyness << "," << rep << "," << year << ",";
	out2.precision(4);
	out2 << MUT << "," << mortality << ",";
	out2.precision(0);
	out2 << K << "," << top << ",";
	out2.precision(4);
	out2 << m_strat << "," << m_age << ",\n";

	out2.close();

	if (year == 5000 && rep == 3) {  // only generate this bit at the end
		ofstream out1;    // output values from an individual
		out1.open(FILE1, ios::app);  //  ios::trunc
		out1.precision(0);

		if (!header1_done) {   // only if header not already written
			out1 << "sticky,rep,time,mutation,mortality,K,pop,x,y,strat,age\n";
			header1_done = true;
		}

		for (int i = 0; i<top; i++) {  // loop through all individuals
			out1 << stickyness << "," << rep << "," << year << ",";
			out1.precision(4);
			out1 << MUT << "," << mortality << ",";
			out1.precision(0);
			out1 << K << "," << top << ",";
			out1.precision(4);
			out1 << inds[i].x << "," << inds[i].y << "," << inds[i].strategy << "," << inds[i].age << ",\n";
		}

		out1.close();
	}  // end of the individual output bit

}  // end output() function


void initialise(void) {

	for (int x = 0; x<SX; x++) for (int y = 0; y<SY; y++) cells[x][y].pop = 0;

	for (int i = 0; i<NINDS; i++) {
		inds[i].x = 5.5;  // just stick any old values in here at the start, this bit is probably not needed
		inds[i].y = 5.5;
		inds[i].strategy = (float)0.1;
		inds[i].age = (float)0.1;
	}

	for (int i = 0; i<START_POP; i++) {  // usually 5000 start population

		inds[i].age = R1(generator) * 10;   // start with a uniform random age 0 to 10
		inds[i].x = (pick_cellX(generator) + R1(generator));
		inds[i].y = (pick_cellX(generator) + R1(generator));
		inds[i].strategy = R1(generator);  // start with a random strategy 0 to 1

		cells[(int)inds[i].x][(int)inds[i].y].pop++;  // add one to the local population

	}

	top = START_POP;   // this is the address of the first available space
}
