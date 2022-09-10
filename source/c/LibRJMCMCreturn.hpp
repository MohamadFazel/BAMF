//Function, Structure and Class Declarations
#define pi 3.141592
#include "math.h"
#include <stdlib.h>
#include <vector>
#include <time.h>       /* time */

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
//urand() returns a random number in the interval [0 1].
inline float urand(){ return (float)rand() / ((float)RAND_MAX + 1); };
//normpdf() returns the PDF-value for the input with the given center and sigma. 
inline float normpdf(float x, float u, float s){
	return 1.0f / ((float)sqrt(2.0f*pi)*s)*(float)exp(-(x - u)*(x - u) / (2 * s*s));
};
//randn() generates a random number from normal dist.    
inline float randn(){
	float tmp = 0; for (int nn = 0; nn<12; nn++)tmp = tmp + urand();
	return tmp - 6.0f;
};
//kround() rounds the input number.
inline int kround(float a){ return (int)floor(a + .5); };
//?
inline int krem(int a, int b){ return a - b*(a / b); };

//Structure containing some input parameters.
struct MCMC_Parameters {
	float PSFsigma;      //Size of the PSF.
	int   N_Trials;      //Number of the jumps in the second part of the chain, which will be the output 
	int   N_Burnin;      //Number of the jumps in the first part of the chain, which will be thrown away    
	float I_stdFg, I_stdBg;         //The standard deviation of the normal-PDF which the random jumps in I are being taken from.
	float X_stdFg, X_stdBg;         //The standard deviation of the norma-PDF which the random jumps in position are taken from.
	float *P_Burnin, *P_Trials;    //The probability of proposing a birth.
	std::vector<float> PV_Burnin;
	std::vector<float> PV_Trials;
	std::vector<float> P;
	std::vector<int> PSFsizeV;
	float Rho;           //Mean number of enitters in ROI 
	float Icutoff;
	float Split_std;     //The standard deviation of the norm-PDf which the random Us are being taken from.
	float Bnd_out;       //The size of the area outside the image, where a proposed particle can be placed.
	int Grid_Zoom;       //The magnification of the output image.
	//float Split_std;     //The std of the normal distributions where you pick Us from.
	float FirstCall = 1;
	float BG_std, ABBG_std;
	float DX;
	float DriftX, DriftY;
	int *PSFsize;
	int IsBackg;
};
//Stat-class deals with the input PDFarray for Intensity prior. 
class Stat{
public:
	Stat();                            //Constructor, declares the properties.
	Stat(float*, float, int);          //Constructor, assigns generated values to the declared properties.
	float findPDF(float I);            //Finds the PDF-value of the given input.
	float genRand();                   //Generates a random number from the given PDFarray.
	std::vector<float> PDFarray;       //Stores the given PDFarray.
	std::vector<float> CDFarray;       //Stores the generated CDFarray.
	std::vector<float> Xarray;         //Stores the generated Xarray.
	std::vector<float> UniqCDFarray;   //Stores the generated UniqCDFarray.
	std::vector<float> UniqXarray;     //Stores the generated UniqXarray.
};
//Emitters-class deals with the current and proposed states of the chain.
class Emitters{
public:
	Emitters();                                 //Constructor, declaring the properties.
	Emitters(int, float *, float *, float *, float, float, float);   //Constructor, assigning the input values to the properties.
	int   N;                                    //Number of the particles.                   
	float BG, ABG, BBG;
	std::vector<float> I;                       //Intensities (number of the photons)
	std::vector<float> X;                       //X-positions
	std::vector<float> Y;                       //Y-positions
	std::vector<int> Signal;                    //Signal = 1 or 0; 1 means part of signal and 0 means part of background
	void add(float, float, float, int);              //Used in the birth and split to add new particles features to the properties of the class.
	void remove(int);                           //Used in the death and merge to remove a particle from the properties of the class.
	int test(int, int, int);             //Check if the proposed new particles are not outside the allowed range (position and intensity).
	int closest(int, float);                           //Used in the merge-method to find the closest particle to be merge with the given particle.
	std::vector<int> closest2(int, float);            //returns 2 close emitters to the given emitter.
	std::vector<int> closestMulti(int, float);        //returns a random number of close emitters to the given emitter.
	float distance(int, int);                    //Used in the merge-method to find the distance betwenn the two-particles that we want to merge.
	float sum();                                //Used in the calculations to find the some of the elements of the input vector.
};
//Stencil-class computes a single row and a single column for the given state that can be comined to get the final model.
class Stencil{
public:
	Stencil();                         //Constructor, declaring the properties.
	Stencil(Emitters &, int, int, float);   //Constructor, assigning the calculated values to the properties of the class.
	int   N;                           //Number of the particles.
	std::vector<float> N_X;            //The calculated row. 
	std::vector<float> N_Y;            //The calculated column.
};
//Model-class calculates the model for the given state
class Model{
public:
	Model();                                                                      //Constructor, declaring the properties.
	Model(Emitters & E_Active, int Xsize, int Ysize, float sigma, std::vector<float> sCMOS);  //Constructor, calls Stencil() and combining the outputs of that to find the model.

	int   Xsize;                                                                  //Size of the image along the X-axis.
	int   Ysize;                                                                  //Size of the image along the Y-axis.
	std::vector<float> M;                                                         //The calculated model.
};
class NumModel{
public:
	NumModel();
	NumModel(Emitters & E_Active, int Xsize, int Ysize, std::vector<int> PSFsize, std::vector<float> SampledPSF, std::vector<float> sCMOS);

	int Xsize;
	int Ysize;
	std::vector<float> M;
};
//This is chain class with methods to add to chain. This is mostly for diagnosis purposes.
class RJMCMCChain{
public:
	RJMCMCChain();                                                             //Constructor, declaring the class properties.
	//RJMCMCChain(int ChainLength);                                              //Constructor?
	void addEmitters(Emitters & E, float LLR, float PR, int Type, int Accept, MCMC_Parameters Params);  //we only need this while running chain.
	//some things we might want to help export chain
	Emitters getEmitters(int TrialIndex);                                       // parse vectors and create an Emitter  
	float getAcceptance_Jump();                                                //return fraction of accepted 
	float getAcceptance_Split();                                               //return fraction of accepted 
	float getAcceptance_Merge();											   //return fraction of accepted 
	float getAcceptance_Birth();											   //return fraction of accepted 
	float getAcceptance_Death();											   //return fraction of accepted 
	int ChainLength;                                                           //this will be NTrials from 
	//These will all be NTrials long
	std::vector<float> LLR;                                                    //Likelihood ratio  Need to resize in constructor to NTrials
	std::vector<float> PR;                                                     //Posterior ratio  Need to resize in constructor  to NTrials
	//std::vector<float> PrI;                                                    //Prior ratio of intensity
	//std::vector<float> PrX;                                                    //Prior ration of Position
	//std::vector<float> Us;                                                     //related to auxiliary parameters Us
	//std::vector<float> Jacob;                                                  //Jacobian
	std::vector<int> JumpType;                                                 //Jump type 1=jump,2-split,3=merge,4-birth,5-death. resize to NTrials
	std::vector<int> N;                                                        //Number of emitters this step  Need to resize in constructor  to NTrials
	std::vector<int> Accepted;                                                 //0,1. 1 for accepted  Need to resize in constructor  to NTrials
	std::vector<int> Signal;                                                   //1 means part of the signal and 0 means part of the background
	//These will all be at least NTrials long, but probably longer.  resize to NTrials in contrructor, but they can grow
	// e.g. X =         [ 1.3, 1.5, 1.7] 1 emitter in first trial, 2 in second trial 
	//      Index=      [1, 2, 2]
	std::vector<int> Index;                                                   //Trial index
	std::vector<float> Photons;                                               //photons
	std::vector<float> X;                                                     //X
	std::vector<float> Y;                                                     //Y
	std::vector<float> BG;
	std::vector<float> ABG;
	std::vector<float> BBG;
};
//RJMCMC-class uses the following methods with the methods from other classes to generate the chain.
class RJMCMC{
	int Xsize;                                    //The size of the image along the X-axis.   
	int Ysize;                                    //The size of the image along the Y-axis.
	int OptModel;
	//float *SampledPSF;
	//float * Data;                                 //The input image
	//float * PGrid;                                //The output image
	//float * PBack;                                //The output background 
	std::vector<float> Data;
	std::vector<float> SampledPSF;
	std::vector<float> sCMOScorr;                            //The SCMOS camera varriance correction
	float MeanN;                                  //Mean number of the particles in the input image.
	int     MoveType;                             //0 jump,1 split, 2 merge, 3 birth, 4 death
	void    find_Birth(float*, float*, float*);       //Finds the most likely place for a new particle using the residual image.
	float   likelihood_ratio(Emitters &, std::vector<float>, MCMC_Parameters);    //Calculating the likelihood ratio
	void    update_grids();                       //Saves the accepted chain on the grids (image).
	void    jump();                               //Inside model jumps.
	void    split();                              //Split.
	void    merge();                              //Merge.
	void    birth();                              //Birth.
	void    death();                              //Death.
	void    gSplit();
	void    gMerge();
	void    backFore();                           //jumping between background and foreground
	Model	M;                                    //Model of the current state, which is saved from the previous calculations
	Model   M_Test;                               //Model of the proposed state.
	NumModel NumM;
	NumModel NumM_Test;
public:
	RJMCMC(float*, MCMC_Parameters, int, int, float*, float, float*, int);        //Constructor, declaring some properties.
	~RJMCMC();
	std::vector<float> PGridV;
	std::vector<float> PBackV;//Constructor, assigning the inputs to the properties of the class.
	void  start_Chain();                                                //jump(), split(), merge(), birth() and death() are being called inside this method. 
	Emitters E_Active;                                                  //The output of the Emitters class, which is the current state.
	MCMC_Parameters Params;                                             //The structure containing some input parameters.
	int     IsBurnin;                                                   //0 no, 1 yes
	int     N_jump, N_split, N_merge, N_birth, N_death;                 //Kepping track of the number of the different proposed jumps
	int     A_jump, A_split, A_merge, A_birth, A_death;                 //Keeping track of the number of different accepted jumps.
	int     N_jumpFg, N_splitFg, N_mergeFg, N_birthFg, N_deathFg;
	int     A_jumpFg, A_splitFg, A_mergeFg, A_birthFg, A_deathFg;
	int     N_jumpBg, N_splitBg, N_mergeBg, N_birthBg, N_deathBg;
	int     A_jumpBg, A_splitBg, A_mergeBg, A_birthBg, A_deathBg;
	RJMCMCChain   Chain;                          //The accepted chain.
	RJMCMCChain   Chain_Test;                     //The proposed chain.
	Stat OffgStat;
	Stat BackStat;
	Stat SignalStat;
};



