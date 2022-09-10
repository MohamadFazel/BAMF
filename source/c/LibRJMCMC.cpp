#include <stdio.h>
#include <math.h>
#ifdef __unix__
#include <mex.h>
#else
#include <mex.h>
#endif
#include "LibRJMCMCreturn.hpp"
//#define mexPrintf printf
//Constructor, declaring properties..
Stencil::Stencil(){
	N = 0;                                   //Number of the particles
};
//Constructor, 
Stencil::Stencil(Emitters & E, int Xsize, int Ysize, float sigma){
	N = E.N;
	float Constant = 1/(sigma*sqrt(2.0f));
	if (N){
		//make N x Xsize Stencil.
		N_X.resize(N*Xsize);
		for (int nn = 0; nn<N; nn++)
			for (int xx = 0; xx < Xsize; xx++)
			{
			//N_X.at(Xsize*nn + xx) = E.I.at(nn)*normpdf(xx + 0.5f, E.X.at(nn), sigma);
			N_X.at(Xsize*nn + xx) = E.I.at(nn)*0.5f*(erf((xx - E.X.at(nn) + 1.0f) * Constant) - erf((xx - E.X.at(nn)) * Constant));
			}
		//make N x Ysize Stencil 
		N_Y.resize(N*Ysize);
		for (int nn = 0; nn<N; nn++)
			for (int xx = 0; xx < Ysize; xx++)
			{
			//N_Y.at(Ysize*nn + xx) = normpdf(xx + 0.5f, E.Y.at(nn), sigma);
			N_Y.at(Ysize*nn + xx) = 0.5f*(erf((xx - E.Y.at(nn) + 1.0f) * Constant) - erf((xx - E.Y.at(nn)) * Constant));
			}
	}
};
//Constructor, declaring properties. 
Model::Model(){
};

//Constructor, calculating model
Model::Model(Emitters & E_A, int Xsize, int Ysize, float sigma, std::vector<float> sCMOScorr){

	//make Stencil.
	Stencil Sten_Active = Stencil(E_A, Xsize, Ysize, sigma);

	M.resize(Xsize*Ysize);

	//Add background
	for (int xx = 0; xx < Xsize; xx++)
		for (int yy = 0; yy < Ysize; yy++)
			M.at(yy + xx*Ysize) = E_A.BG + E_A.ABG*xx + E_A.BBG*yy + sCMOScorr.at(xx*Ysize + yy);

	//Add particles
	for (int nn = 0; nn<Sten_Active.N; nn++)
		for (int yy = 0; yy<Ysize; yy++)
			for (int xx = 0; xx<Xsize; xx++)
				M.at(yy + xx*Ysize) += Sten_Active.N_X.at(Xsize*nn + xx) * Sten_Active.N_Y.at(Ysize*nn + yy);
};


//Constrictor, declaring properties
NumModel::NumModel(){
};

//Constructor, calculating model
NumModel::NumModel(Emitters & E_A, int Xsize, int Ysize, std::vector<int> PSFsize, std::vector<float> SampledPSF, std::vector<float> sCMOScorr){
	float YRatio, XRatio;
	double Xinst, Yinst, OverSampling; //In the calculations in this function, we sometimes need double precision.
	int STX2, IMX2, STY2, IMY2, ST0, STX, STY, IMX, IMY, InnerInd, STXpxxPSF0, STX2pxxPSF0, PSF0123, PSF012, PSF01, PSF0, Ind1out, Ind2out, Ind3out, Ind4out;
	M.resize(Xsize*Ysize);
	OverSampling = (double)PSFsize.at(3);
	//adding the background
	for (int xx = 0; xx < Xsize; xx++)
	{
		InnerInd = xx*Ysize;
		for (int yy = 0; yy < Ysize; yy++)
			M.at(yy + InnerInd) = E_A.BG + E_A.ABG*xx + E_A.BBG*yy + sCMOScorr.at(yy + InnerInd);
	}

	PSF0123 = PSFsize.at(0) * PSFsize.at(1) * PSFsize.at(2) * PSFsize.at(3);
	PSF012 = PSFsize.at(0) * PSFsize.at(1) * PSFsize.at(2);
	PSF01 = PSFsize.at(0) * PSFsize.at(1);
	PSF0 = PSFsize.at(0);
	for (int nn = 0; nn < E_A.N; nn++)
	{
		Xinst = (double)E_A.X.at(nn);
		Yinst = (double)E_A.Y.at(nn);
		ST0 = (int)(PSFsize.at(0) / 2);
		STX = (int)ceil(ST0 - Xinst);
		IMX = (int)floor(OverSampling*(Xinst - floor(Xinst)));
		STY = (int)ceil(ST0 - Yinst);
		IMY = (int)floor(OverSampling*(Yinst - floor(Yinst)));

		if (IMX == OverSampling - 1)
		{
			STX2 = STX - 1;
			IMX2 = 0;
		}
		else
		{
			STX2 = STX;
			IMX2 = IMX + 1;
		}
		if (IMY == OverSampling - 1)
		{
			STY2 = STY - 1;
			IMY2 = 0;
		}
		else
		{
			STY2 = STY;
			IMY2 = IMY + 1;
		}

		Ind1out = IMY*PSF01 + IMX*PSF012;
		Ind2out = IMY*PSF01 + IMX2*PSF012;
		Ind3out = IMY2*PSF01 + IMX*PSF012;
		Ind4out = IMY2*PSF01 + IMX2*PSF012;
		YRatio = (float)(1 - ((Yinst - floor(Yinst)) - (double)IMY / OverSampling)*OverSampling);
		XRatio = (float)(1 - ((Xinst - floor(Xinst)) - (double)IMX / OverSampling)*OverSampling);
		for (int xx = 0; xx < Xsize; xx++)
		{
			InnerInd = xx*Ysize;
			STXpxxPSF0 = (STX + xx)*PSF0;
			STX2pxxPSF0 = (STX2 + xx)*PSF0;
			for (int yy = 0; yy < Ysize; yy++)
			{
				int Ind1, Ind2, Ind3, Ind4, Ind5;
				Ind1 = (STY + yy) + STXpxxPSF0 + Ind1out;
				if (Ind1 >= PSF0123 || Ind1 < 0)
				{
					mexPrintf("Ind1 is large: %d, X: %g, Y: %g, ST0: %d, STY: %d, STX: %d, IMY: %d, IMX: %d\n", Ind1, Xinst, Yinst, ST0, STY, STX, IMY, IMX);
					mexEvalString("drawnow;");
					Ind1 = PSF0123 - 1;
				}
				Ind2 = (STY + yy) + STX2pxxPSF0 + Ind2out;
				if (Ind2 >= PSF0123 || Ind2 < 0)
				{
					mexPrintf("Ind2 is large: %d, X: %g, Y: %g, ST0: %d, STY: %d, STX2: %d, IMY: %d, IMX2: %d\n", Ind2, Xinst, Yinst, ST0, STY, STX2, IMY, IMX2);
					mexEvalString("drawnow;");
					Ind2 = PSF0123 - 1;
				}
				Ind3 = (STY2 + yy) + STXpxxPSF0 + Ind3out;
				if (Ind3 >= PSF0123 || Ind3 < 0)
				{
					mexPrintf("Ind3 is large: %d, X: %g, Y: %g, ST0: %d, STY2: %d, STX: %d, IMY2: %d, IMX: %d\n", Ind3, Xinst, Yinst, ST0, STY2, STX, IMY2, IMX);
					mexEvalString("drawnow;");
					Ind3 = PSF0123 - 1;
				}
				Ind4 = (STY2 + yy) + STX2pxxPSF0 + Ind4out;
				if (Ind4 >= PSF0123 || Ind4 < 0)
				{
					mexPrintf("Ind4 is large: %d, X: %g, Y: %g, ST0: %d, STY2: %d, STX2: %d, IMY2: %d, IMX2: %d\n", Ind4, Xinst, Yinst, ST0, STY2, STX2, IMY2, IMX2);
					mexEvalString("drawnow;");
					Ind4 = PSF0123 - 1;
				}
				Ind5 = InnerInd + yy;
				if (Ind5 >= Ysize*Xsize || Ind5 < 0)
				{
					mexPrintf("Ind5 is large: %d\n", Ind5);
					exit(0);
				}
				float A = XRatio*SampledPSF.at(Ind1) + (1.0f - XRatio)*SampledPSF.at(Ind2);
				float B = XRatio*SampledPSF.at(Ind3) + (1.0f - XRatio)*SampledPSF.at(Ind4);
				M.at(Ind5) = M.at(Ind5) + E_A.I.at(nn) * (YRatio*A + (1.0f - YRatio)*B);
			}
		}
	}
};
//Constructor, declaring properties
Stat::Stat(){

};
//Construcotr
Stat::Stat(float *PDFarray_in, float DX, int LengthPDF){
	std::vector<float> XM;          //An intermediate parameters that assists us to find the repeated values of CDFarray. 
	CDFarray.push_back(0);                                  //The first element of the CDF is always zero.
	Xarray.push_back(0);                                    //The first element of the Xarray is zero.
	XM.push_back(0);

	for (int ii = 0; ii < LengthPDF; ii++)
	{
		PDFarray.push_back(PDFarray_in[ii]); //PDFarray
	}

	for (int ii = 1; ii < LengthPDF; ii++)
	{
		CDFarray.push_back(CDFarray.at(ii - 1) + (PDFarray.at(ii) + PDFarray.at(ii - 1))*DX / 2.0f); //CDF is the integral of PDF.
		if (CDFarray.at(ii) > 1.0f) CDFarray.at(ii) = 1.0f;
		XM.push_back(CDFarray.at(ii));                                                                 //The intermediate parameter.
		Xarray.push_back(Xarray.at(ii - 1) + DX);                                                 //Xarray
	}


	for (int ii = LengthPDF - 2; ii >= 0; ii--)
	{
		if (CDFarray.at(ii) == CDFarray.at(ii + 1)) XM.at(ii + 1) = -1;
	}
	//Filling Uniqarrays with those elements that are not -1.
	int SizUniq = 0;
	for (int ii = 0; ii < LengthPDF; ii++)
	{
		if (XM.at(ii) != -1)
		{
			SizUniq = SizUniq + 1;
			UniqXarray.push_back(Xarray.at(ii));
			UniqCDFarray.push_back(CDFarray.at(ii));
		}
	}

	//The last elemet of the CDFarray must be one, otherwise we might get an error.
	if (UniqCDFarray.at(SizUniq - 1) != 1)
		UniqCDFarray.at(SizUniq - 1) = 1;
	if (UniqCDFarray.at(0) != 0)
		UniqCDFarray.at(0) = 0;
};
//findPDF() calculates the PDF of the input using the given PDFarray.
float Stat::findPDF(float I){
	float IPDF;
	int Indd = -1;
	int LengthPDF = (int)Xarray.size();
	//IF the input is beyond the range of the available array return zero.
	if (I >= Xarray.back() || I <= 0)
		IPDF = 0;
	else
	{
		//Finding the closest element of the PDFarray to the input value.
		for (int ji = 0; ji < LengthPDF; ji++)
		{
			if (Xarray.at(ji) - I > 0)
			{
				Indd = ji;
				break;
			}
		}
		if (Indd == -1 || Indd == 0 || Indd >= Xarray.size() || Indd >= PDFarray.size())
		{
			mexPrintf("Wierd thing happening in findPDF. Indd: %d\n", Indd);
			return 0;
		}
		//The PDF value of the input is simply given by the following interpolation.
		IPDF = (PDFarray.at(Indd) * (I - Xarray.at(Indd - 1)) + PDFarray.at(Indd - 1) * (Xarray.at(Indd) - I)) / (Xarray.at(1) - Xarray.at(0));
	}
	return IPDF;
};
//genRand() generates a random number from the given PDFarray.
float Stat::genRand(){
	float RandVal, I_new;
	int SecInd = -1;
	RandVal = urand();
	int SizUniq = (int)UniqCDFarray.size();
	//Finding the closest element in the CDFarray to the random number RandVal.
	if (RandVal <= 0) return 0;
	else if (RandVal >= 1) return UniqXarray.back();
	else
	{
		for (int ii = 0; ii < SizUniq; ii++)
		{
			if (UniqCDFarray.at(ii) - RandVal > 0)
			{
				SecInd = ii;
				break;
			}
		}
	}
	if (SecInd == -1)
	{
		mexPrintf("Wierd thing happening in genRand. SecInd: %d\n", SecInd);
		return 0;
	}
	else if (SecInd == 0)
		I_new = 0;
	else if (SecInd >= SizUniq)
		I_new = UniqXarray.back();
	else
	{
		//The random number is simply given by the following interpolation.
		I_new = (UniqXarray.at(SecInd - 1) * (UniqCDFarray.at(SecInd) - RandVal) + UniqXarray.at(SecInd) * (RandVal - UniqCDFarray.at(SecInd - 1))) / (UniqCDFarray.at(SecInd) - UniqCDFarray.at(SecInd - 1));
	}
	return I_new;
};

//Constructor, declaring properties
RJMCMCChain::RJMCMCChain(){
	ChainLength = 0;                 //Length of the chain
};

//addEmitters() adds the new state to the end of the chain.
void RJMCMCChain::addEmitters(Emitters & E_in, float LLR_in, float PR_in, int Type_in, int Accept_in, MCMC_Parameters Params_in){
	int StateIndex;
	LLR.push_back(LLR_in);            //adding the likelihood ratio.
	PR.push_back(PR_in);              //adding the likelihood ratio. 
	JumpType.push_back(Type_in);      //adding the type of the jump.
	Accepted.push_back(Accept_in);    //adding the acceptance.
	N.push_back(E_in.N);              //adding the number of the particles.
	BG.push_back(E_in.BG);
	ABG.push_back(E_in.ABG);
	BBG.push_back(E_in.BBG);
	//If the chain is empty then this is the first state.
	if (Index.empty())
		StateIndex = 1;
	else
		StateIndex = Index.back() + 1;  //The state index is the index of the pervious state plus one.
	int NN = E_in.N;
	if (NN < 1)
	{
		Index.push_back(StateIndex);
		Photons.push_back(0);
		X.push_back(-100);
		Y.push_back(-100);
		Signal.push_back(-1);
		return;
	}
	for (int nn = 0; nn < E_in.N; nn++){
		Index.push_back(StateIndex);     //adding the StateIndex
		Photons.push_back(E_in.I.at(nn));   //adding the Photons
		X.push_back(E_in.X.at(nn) + Params_in.DriftX);         //adding the X-position
		Y.push_back(E_in.Y.at(nn) + Params_in.DriftY);         //adding the Y-position
		Signal.push_back(E_in.Signal.at(nn)); //adding signal or background
	}
}
//getEmitters() gets the Chain and the index of the desired state and returns that state. 
Emitters RJMCMCChain::getEmitters(int TrialIndex){
	//If the given index is larger than the size of the chain return.
	if (TrialIndex >= N.size() || TrialIndex >= BG.size() || TrialIndex >= ABG.size() || TrialIndex >= BBG.size())
	{
		float tmp[1];
		tmp[0] = 0;
		Emitters E_NTrials = Emitters(0, tmp, tmp, tmp, 0, 0, 0);
		E_NTrials.Signal.at(0) = 2;
		return E_NTrials;
	}

	int TrialPartNum = N.at(TrialIndex); //Number of the particles in that specific state
	if (TrialPartNum == 0) TrialPartNum = 1;
	int Num = 0;

	for (int nn = 0; nn < TrialIndex; nn++)
	{
		if (N.at(nn) < 1)
			Num = Num + 1;
		else
			Num = Num + N.at(nn);
	}

	if (Num + TrialPartNum > X.size() || Num + TrialPartNum > Y.size() || Num + TrialPartNum > Photons.size() || Num + TrialPartNum > Signal.size())
	{
		float tmp[1];
		tmp[0] = 0;
		Emitters E_NTrials = Emitters(0, tmp, tmp, tmp, 0, 0, 0);
		E_NTrials.Signal.at(0) = 2;
		return E_NTrials;
	}
	//Allocating enough memory to the following auxiliary parameters.
	std::vector<float> IPart;
	std::vector<float> XPart;
	std::vector<float> YPart;
	std::vector<int> SignalPart;
	//Filling the auxiliary parameters 
	for (int nn = 0; nn < TrialPartNum; nn++)
	{
		IPart.push_back(Photons.at(Num + nn));
		XPart.push_back(X.at(Num + nn));
		YPart.push_back(Y.at(Num + nn));
		SignalPart.push_back(Signal.at(Num + nn));
	}
	float BGp, ABGp, BBGp;
	BGp = BG.at(TrialIndex);
	ABGp = ABG.at(TrialIndex);
	BBGp = BBG.at(TrialIndex);
	//Making a structure of Emmiters-class with the found parameters of the desired state.
	Emitters E_NTrials = Emitters(TrialPartNum, &IPart[0], &XPart[0], &YPart[0], BGp, ABGp, BBGp);
	for (int nn = 0; nn < TrialPartNum; nn++)
		E_NTrials.Signal.at(nn) = SignalPart.at(nn);
	return E_NTrials;
}
//To be completed.
float RJMCMCChain::getAcceptance_Jump(){
	int N = (int)Accepted.size();
	float AJ = 0, J = 1;
	for (int ii = 0; ii < N; ii++)
	{
		if (JumpType.at(ii) == 1)
		{
			J = J + 1;
			if (Accepted.at(ii) == 1) AJ = AJ + 1;
		}
	}
	float Ratio = AJ / J;
	return Ratio;
}
//To be completed.
float RJMCMCChain::getAcceptance_Split(){
	int N = (int)Accepted.size();
	float AJ = 0, J = 1;
	for (int ii = 0; ii < N; ii++)
	{
		if (JumpType.at(ii) == 2)
		{
			J = J + 1;
			if (Accepted.at(ii) == 1) AJ = AJ + 1;
		}
	}
	float Ratio = AJ / J;
	return Ratio;
}
//To be completed.
float RJMCMCChain::getAcceptance_Merge(){
	int N = (int)Accepted.size();
	float AJ = 0, J = 1;
	for (int ii = 0; ii < N; ii++)
	{
		if (JumpType.at(ii) == 3)
		{
			J = J + 1;
			if (Accepted.at(ii) == 1) AJ = AJ + 1;
		}
	}
	float Ratio = AJ / J;
	return Ratio;
}
//To be completed.
float RJMCMCChain::getAcceptance_Birth(){
	int N = (int)Accepted.size();
	float AJ = 0, J = 1;
	for (int ii = 0; ii < N; ii++)
	{
		if (JumpType.at(ii) == 4)
		{
			J = J + 1;
			if (Accepted.at(ii) == 1) AJ = AJ + 1;
		}
	}
	float Ratio = AJ / J;
	return Ratio;
}
//To be completed.
float RJMCMCChain::getAcceptance_Death(){
	int N = (int)Accepted.size();
	float AJ = 0, J = 1;
	for (int ii = 0; ii < N; ii++)
	{
		if (JumpType.at(ii) == 5)
		{
			J = J + 1;
			if (Accepted.at(ii) == 1) AJ = AJ + 1;
		}
	}
	float Ratio = AJ / J;
	return Ratio;
}
//Constructor, declaring some properties
Emitters::Emitters(){
	int N = 0;
};

//Constructor
Emitters::Emitters(int N_in, float *I_in, float *X_in, float *Y_in, float BG_in, float ABG_in, float BBG_in){
	N = N_in;
	I.assign(I_in, I_in + N);
	X.assign(X_in, X_in + N);
	Y.assign(Y_in, Y_in + N);
	Signal.assign(N, 1);
	BG = BG_in;
	ABG = ABG_in;
	BBG = BBG_in;
};
//add() is called in the split and birth to add the new particle to the current state and producing the new state.
void Emitters::add(float I_in, float x, float y, int bf){
	I.push_back(I_in);   //adding the new intensity
	X.push_back(x);      //adding the new X-position
	Y.push_back(y);      //adding the new Y-position.
	Signal.push_back(bf); //always take the new amitter as part of the signal.
	N = N + 1;               //adding one to the number of the particles.
};
//remove() is called in the death and merge to remove one particle from the current state and making the new state.
void Emitters::remove(int id){
	switch (N){
	case 0:                       //return if there is no particle to remove.
		return;
		break;
	default:
		I.erase(I.begin() + id);   //removing the intensity
		X.erase(X.begin() + id);   //removing the X-position
		Y.erase(Y.begin() + id);   //removing the Y-position
		Signal.erase(Signal.begin() + id); //removing the signal
		N = N - 1;                   //subtracting one particle.
	}
};
//test() see if the intensity and the position of the new particles are not outside of the allowed range.  
int Emitters::test(int Bnd_out, int Xsize, int Ysize){
	//Intensity must be positive
	for (int nn = 0; nn<N; nn++)
		if (I.at(nn)<0) return 0;

	// x left boundary
	for (int nn = 0; nn<N; nn++)
		if (X.at(nn)<(-Bnd_out)) return 0;

	//y upper boundary
	for (int nn = 0; nn < N; nn++)
		if (Y.at(nn) < (-Bnd_out)) return 0;

	//x right boundary
	for (int nn = 0; nn<N; nn++)
		if (X.at(nn)>(Xsize + Bnd_out)) return 0;

	//y bottom boundary
	for (int nn = 0; nn<N; nn++)
		if (Y.at(nn)>(Ysize + Bnd_out)) return 0;

	return 1;
};
//used in the calculations to add up the elements of the given vector.
float Emitters::sum(){
	float tmp = 0;

	for (int nn = 0; nn<N; nn++)
		tmp += I.at(nn);

	return tmp;
};
//finding the closest particle to the given particle to be merged with that.
int Emitters::closest(int jj, float dmin){
	//float dmin = 1000; //bad
	int ii = -1;
	std::vector<int> Ind;
	std::vector<int> Ind1;
	int Num = 0;
	//searching over all the existing particles in the current state in order to find the closest one
	for (int nn = 0; nn<N; nn++)
		if (nn != jj){
		float d = sqrt(pow(X.at(nn) - X.at(jj), 2) + pow(Y.at(nn) - Y.at(jj), 2));
		if (d<dmin){
			Ind.push_back(nn);
			Num = Num + 1;
		}
		}
	if (Num > 0)
	{
		float Ua = urand();
		int mm = min((int)(Ua*Num), Num - 1);
		ii = Ind.at(mm);
	}
	else
	{
		for (int nn = 0; nn < N; nn++)
			if (nn != jj)
				Ind1.push_back(nn);

		float Ua = urand();
		int mm = min((int)(Ua*(N - 1)), N - 2);
		ii = Ind1.at(mm);
	}

	return ii;
};
std::vector<int> Emitters::closest2(int jj, float dmin){
	std::vector<int> Closest;
	Closest.push_back(-1);
	Closest.push_back(-1);
	std::vector<int> Ind;
	std::vector<int> Ind1;
	int Num = 0, mm, pp;
	//searching over all the existing particles in the current state in order to find the closest one
	for (int nn = 0; nn<N; nn++)
		if (nn != jj){
		float d = sqrt(pow(X.at(nn) - X.at(jj), 2) + pow(Y.at(nn) - Y.at(jj), 2));
		if (d<dmin){
			Ind.push_back(nn);
			Num = Num + 1;
		}
		}
	if (Num > 0)
	{
		float Ua = urand();
		mm = min((int)(Ua*Num), Num - 1);
		Closest.at(0) = Ind.at(mm);
		Ind.erase(Ind.begin() + mm);
	}
	else
	{
		for (int nn = 0; nn < N; nn++)
			if (nn != jj)
				Ind1.push_back(nn);

		float Ua = urand();
		mm = min((int)(Ua*(N - 1)), N - 2);
		Closest.at(0) = Ind1.at(mm);
		Ind1.erase(Ind1.begin() + mm);
	}
	if (Num > 1)
	{
		float Ua = urand();
		pp = min((int)(Ua*(Num - 1)), Num - 2);
		Closest.at(1) = Ind.at(pp);
	}
	else
	{
		if (Num == 1){
			for (int nn = 0; nn < N; nn++)
				if (nn != jj && nn != mm)
					Ind1.push_back(nn);
		}
		float Ua = urand();
		pp = min((int)(Ua*(N - 2)), N - 3);
		if (pp > Ind1.size())
			int adfg = 0;
		Closest.at(1) = Ind1.at(pp);
	}
	return Closest;
}

std::vector<int> Emitters::closestMulti(int jj, float dmin){
	std::vector<int> Closest;
	std::vector<int> Ind;
	int Num = 0, mm, pp, Ind1;
	//searching over all the existing particles in the current state in order to find the closest one
	for (int nn = 0; nn<N; nn++)
		if (nn != jj){
		float d = sqrt(pow(X.at(nn) - X.at(jj), 2) + pow(Y.at(nn) - Y.at(jj), 2));
		if (d<dmin){
			Ind.push_back(nn);
			Num = Num + 1;
		}
		}

	if (Num > 0)
	{
		int Np;
		float Ua = urand();
		if (Ua < 1.0f / float(Num))
		{
			Np = Num;
		}
		else
		{
			Np = min((int)(Ua*Num), Num - 1);
		}
		for (int nn = 0; nn < Np; nn++)
		{
			Ua = urand();
			Ind1 = min(int(Ua*(Num - nn)), (Num - nn - 1));
			Closest.push_back(Ind.at(Ind1));
			Ind.erase(Ind.begin() + Ind1);
		}
	}
	return Closest;
}

//used to find the distance between two given particles.
float Emitters::distance(int jj, int ii){
	float d = pow(X.at(ii) - X.at(jj), 2) + pow(Y.at(ii) - Y.at(jj), 2);
	return sqrt(d);
};
//Constructor, declaring some properties
RJMCMC::RJMCMC(float *Data_in, MCMC_Parameters Params_in, int Xsize_in, int Ysize_in, float *sCMOScorr_in, float MeanN_in, float *SampledPSF_in, int OptModel_in){                                      //input image
	int ArrNumelZoom = Xsize_in*Ysize_in*Params_in.Grid_Zoom*Params_in.Grid_Zoom;
	int ArrNumel = Xsize_in*Ysize_in;
	int PSFNumel = Params_in.PSFsizeV.at(0) * Params_in.PSFsizeV.at(1) * Params_in.PSFsizeV.at(2) * Params_in.PSFsizeV.at(3);
	Data.assign(Data_in, Data_in + ArrNumel);
	sCMOScorr.assign(sCMOScorr_in, sCMOScorr_in + ArrNumel);
	SampledPSF.assign(SampledPSF_in, SampledPSF_in + PSFNumel);
	PGridV.resize(ArrNumelZoom);
	PBackV.resize(ArrNumelZoom);
	Params = Params_in;                                    //structure containing some input parameters
	Xsize = Xsize_in;                                      //size of the input image along the X-axis
	Ysize = Ysize_in;                                      //size of the input image along the Y-axis                          //SCMOS camera variance correction
	MeanN = MeanN_in;                                    //Mean number of the particles in the input image
	OptModel = OptModel_in;
	N_jump = 1; N_split = 1; N_merge = 1; N_birth = 1; N_death = 1;    //number of the proposed jumps of each type
	A_jump = 0; A_split = 0; A_merge = 0; A_birth = 0; A_death = 0;    //number of the accepted jumps of each type.
	N_jumpFg = 1; N_splitFg = 1; N_mergeFg = 1; N_birthFg = 1; N_deathFg = 1;
	N_jumpBg = 1; N_splitBg = 1; N_mergeBg = 1; N_birthBg = 1; N_deathBg = 1;
	A_jumpFg = 0; A_splitFg = 0; A_mergeFg = 0; A_birthFg = 0; A_deathFg = 0;
	A_jumpBg = 0; A_splitBg = 0; A_mergeBg = 0; A_birthBg = 0; A_deathBg = 0;
};
//Constructor
RJMCMC::~RJMCMC(){

};

//jump(), merge(), split(), birth() and death() are implemented inside start_Chain().
void RJMCMC::start_Chain(){
	//Burnin, which is the first part of the chain and we do not save that.
	IsBurnin = 1;
	Params.P.resize(8);
	for (int pp = 0; pp < 8; pp++)
		Params.P.at(pp) = Params.P_Burnin[pp];
	for (int nn = 0; nn<Params.N_Burnin; nn++){
		//proposing either a jump, split/merge or death/birth
		float rtmp = urand();
		if (rtmp < Params.PV_Burnin.at(0)) {
			jump();
		}
		else if (rtmp >= Params.PV_Burnin.at(0) && rtmp < Params.PV_Burnin.at(1)) {
			split();
		}
		else if (rtmp >= Params.PV_Burnin.at(1) && rtmp < Params.PV_Burnin.at(2)) {
			merge();
		}
		else if (rtmp >= Params.PV_Burnin.at(2) && rtmp < Params.PV_Burnin.at(3)) {
			birth();
		}
		else if (rtmp >= Params.PV_Burnin.at(3) && rtmp < Params.PV_Burnin.at(4)){
			death();
		}
		else if (rtmp >= Params.PV_Burnin.at(4) && rtmp < Params.PV_Burnin.at(5)){
			gSplit();
		}
		else if (rtmp >= Params.PV_Burnin.at(5) && rtmp < Params.PV_Burnin.at(6)){
			gMerge();
		}
		else{
			backFore();
		}
	}
	//The second part of the chain where it has converged and will be returned as output.
	IsBurnin = 0;
	for (int pp = 0; pp < 8; pp++)
		Params.P.at(pp) = Params.PV_Trials[pp];
	for (int nn = 0; nn<Params.N_Trials; nn++){
		float rtmp = urand();
		if (rtmp < Params.PV_Trials.at(0)) {
			jump();
		}
		else if (rtmp >= Params.PV_Trials.at(0) && rtmp < Params.PV_Trials.at(1)) {
			split();
		}
		else if (rtmp >= Params.PV_Trials.at(1) && rtmp < Params.PV_Trials.at(2)) {
			merge();
		}
		else if (rtmp >= Params.PV_Trials.at(2) && rtmp < Params.PV_Trials.at(3)) {
			birth();
		}
		else if (rtmp >= Params.PV_Trials.at(3) && rtmp < Params.PV_Trials.at(4)){
			death();
		}
		else if (rtmp >= Params.PV_Burnin.at(4) && rtmp < Params.PV_Burnin.at(5)){
			gSplit();
		}
		else if (rtmp >= Params.PV_Burnin.at(5) && rtmp < Params.PV_Burnin.at(6)){
			gMerge();
		}
		else{
			backFore();
		}
	}
};

void RJMCMC::jump(){
	// within model jump

	N_jump++;

	int JumpType = 1;    //jump type
	float LLR = 0;       //Likelihood ratio
	float PR = 0;        //Posterior ratio
	int Accept = 0;      //Acceptance
	float PrI = 1;       //Intensity prior ration
	float PrBg = 1;      //Background prior ratio
	float D;
	float IPDF = 1, IPDFtest = 1;
	Emitters E_Active_Test = E_Active;

	try{

		if (E_Active.N < 1)return;

		//Picking a particle randomly
		float Ua = urand();
		int nn = min((int)(Ua*E_Active_Test.N), E_Active_Test.N - 1);
		float Ub = urand();
		if (Ub < 0.25)
		{
			//Jump in offset background
			JumpType = 2;
			E_Active_Test.BG = E_Active_Test.BG + Params.BG_std*randn();     //proposing a new value for the background
			E_Active_Test.ABG = E_Active_Test.ABG + Params.ABBG_std*randn(); //added for nonuniform background
			E_Active_Test.BBG = E_Active_Test.BBG + Params.ABBG_std*randn(); //added for nonuniform background
		}
		else if (Ub < 0.50 && Ub >= 0.25)
		{
			//finding a group of close particles and move a random subset of them so that the center of mass is conserved and the number of photons are conserved.
			if (E_Active.N < 2) return;
			JumpType = 10;
			std::vector<int> Closest;
			Closest = E_Active.closestMulti(nn, 2 * Params.PSFsigma);
			float I0;
			int Np = Closest.size();
			if (Np == 0) return;
			I0 = E_Active.I.at(nn);
			for (int ii = 0; ii < Np; ii++)
				I0 = I0 + E_Active.I.at(Closest.at(ii));
			//redistributing the photons
			float II0 = 0, XI0 = 0, YI0 = 0, XII0 = 0, YII0 = 0;
			for (int ii = 0; ii < Np; ii++)
			{
				if (E_Active.Signal.at(Closest.at(ii)) == 1)
				{
					E_Active_Test.X.at(Closest.at(ii)) = E_Active.X.at(Closest.at(ii)) + Params.X_stdFg*randn();
					E_Active_Test.Y.at(Closest.at(ii)) = E_Active.Y.at(Closest.at(ii)) + Params.X_stdFg*randn();
					//redistributing the photons among the emitters
					E_Active_Test.I.at(Closest.at(ii)) = (I0 / float(Np + 1))*(0.15f*randn() + 1);
					// calculation of intensity prior for the proposed state(Signal)
					IPDFtest = SignalStat.findPDF(E_Active_Test.I.at(Closest.at(ii)));
					//calculation of intensity prior for the current state (Signal)
					IPDF = SignalStat.findPDF(E_Active.I.at(Closest.at(ii)));
					PrI = PrI * IPDFtest / IPDF;
				}
				else
				{
					E_Active_Test.X.at(Closest.at(ii)) = E_Active.X.at(Closest.at(ii)) + Params.X_stdBg*randn();
					E_Active_Test.Y.at(Closest.at(ii)) = E_Active.Y.at(Closest.at(ii)) + Params.X_stdBg*randn();
					//redistributing the photons among the emitters
					E_Active_Test.I.at(Closest.at(ii)) = (I0 / float(Np + 1))*(0.15f*randn() + 1);
					//calculation of intensity prior for the proposed state (background)
					IPDFtest = BackStat.findPDF(E_Active_Test.I.at(Closest.at(ii)));
					//calculation of intensity prior for the current state (background)
					IPDF = BackStat.findPDF(E_Active.I.at(Closest.at(ii)));
					PrI = PrI * IPDFtest / IPDF;
				}
				XI0 = XI0 + E_Active.I.at(Closest.at(ii))*E_Active.X.at(Closest.at(ii));
				YI0 = YI0 + E_Active.I.at(Closest.at(ii))*E_Active.Y.at(Closest.at(ii));
				XII0 = XII0 + E_Active_Test.I.at(Closest.at(ii))*E_Active_Test.X.at(Closest.at(ii));
				YII0 = YII0 + E_Active_Test.I.at(Closest.at(ii))*E_Active_Test.Y.at(Closest.at(ii));
				II0 = II0 + E_Active_Test.I.at(Closest.at(ii));
			}
			E_Active_Test.I.at(nn) = I0 - II0;
			E_Active_Test.X.at(nn) = (XI0 + E_Active.I.at(nn)*E_Active.X.at(nn) - XII0) / E_Active_Test.I.at(nn);
			E_Active_Test.Y.at(nn) = (YI0 + E_Active.I.at(nn)*E_Active.Y.at(nn) - YII0) / E_Active_Test.I.at(nn);
			if (E_Active.Signal.at(nn) == 1)
			{
				IPDFtest = SignalStat.findPDF(E_Active_Test.I.at(nn));
				IPDF = SignalStat.findPDF(E_Active.I.at(nn));
				PrI = PrI * IPDFtest / IPDF;
			}
			else
			{
				IPDFtest = BackStat.findPDF(E_Active_Test.I.at(nn));
				IPDF = BackStat.findPDF(E_Active.I.at(nn));
				PrI = PrI * IPDFtest / IPDF;
			}
		}
		else if (Ub >= 0.5 && Ub < 0.75)
		{
			//Finding a group of close particles and make inside model jump for a random subset of them, the intensity and center of mass are not conserved here.
			JumpType = 3;
			std::vector<int> Ind;
			Ind = E_Active.closestMulti(nn, 2 * Params.PSFsigma);
			Ind.push_back(nn);
			int NN = Ind.size();
			for (int ii = 0; ii < NN; ii++)
			{
				if (E_Active.Signal.at(Ind.at(ii)) == 1)
				{
					N_jumpFg++;
					E_Active_Test.I.at(Ind.at(ii)) = E_Active_Test.I.at(Ind.at(ii)) + Params.I_stdFg*randn(); //proposing new intensities for the signal particles
					E_Active_Test.X.at(Ind.at(ii)) = E_Active_Test.X.at(Ind.at(ii)) + Params.X_stdFg*randn(); //proposing new x-positions for the signal particles
					E_Active_Test.Y.at(Ind.at(ii)) = E_Active_Test.Y.at(Ind.at(ii)) + Params.X_stdFg*randn(); // proposing new y-positions for the signal particles
					//calculation of intensity prior for the proposed state (Signal)
					IPDFtest = SignalStat.findPDF(E_Active_Test.I.at(Ind.at(ii)));
					//calculation of intensity prior for the current state (Signal)
					IPDF = SignalStat.findPDF(E_Active.I.at(Ind.at(ii)));
					PrI = PrI * IPDFtest / IPDF;
				}
				else
				{
					N_jumpBg++;
					E_Active_Test.I.at(Ind.at(ii)) = E_Active_Test.I.at(Ind.at(ii)) + Params.I_stdBg*randn();
					E_Active_Test.X.at(Ind.at(ii)) = E_Active_Test.X.at(Ind.at(ii)) + Params.X_stdBg*randn();
					E_Active_Test.Y.at(Ind.at(ii)) = E_Active_Test.Y.at(Ind.at(ii)) + Params.X_stdBg*randn();
					//calculation of intensity prior for the proposed state (background)
					IPDFtest = BackStat.findPDF(E_Active_Test.I.at(Ind.at(ii)));
					//calculation of intensity prior for the current state (background)
					IPDF = BackStat.findPDF(E_Active.I.at(Ind.at(ii)));
					PrI = PrI * IPDFtest / IPDF;
				}
			}
		}
		else
		{
			//picking a single random particle and make an inside model jump for that.
			if (E_Active.Signal.at(nn) == 1)
			{
				N_jumpFg++;
				E_Active_Test.I.at(nn) = E_Active_Test.I.at(nn) + Params.I_stdFg*randn(); //proposing new intensities for the signal particles
				E_Active_Test.X.at(nn) = E_Active_Test.X.at(nn) + Params.X_stdFg*randn(); //proposing new x-positions for the signal particles
				E_Active_Test.Y.at(nn) = E_Active_Test.Y.at(nn) + Params.X_stdFg*randn(); // proposing new y-positions for the signal particles
				//calculation of intensity prior for the proposed state (Signal)
				IPDFtest = SignalStat.findPDF(E_Active_Test.I.at(nn));
				//calculation of intensity prior for the current state (Signal)
				IPDF = SignalStat.findPDF(E_Active.I.at(nn));
				PrI = IPDFtest / IPDF;
			}
			else
			{
				N_jumpBg++;
				E_Active_Test.I.at(nn) = E_Active_Test.I.at(nn) + Params.I_stdBg*randn();
				E_Active_Test.X.at(nn) = E_Active_Test.X.at(nn) + Params.X_stdBg*randn();
				E_Active_Test.Y.at(nn) = E_Active_Test.Y.at(nn) + Params.X_stdBg*randn();
				//calculation of intensity prior for the proposed state (background)
				IPDFtest = BackStat.findPDF(E_Active_Test.I.at(nn));
				//calculation of intensity prior for the current state (background)
				IPDF = BackStat.findPDF(E_Active.I.at(nn));
				PrI = IPDFtest / IPDF;
			}
		}
		//Test if the particle is not outside of the allowed region.
		if (!E_Active_Test.test((int)Params.Bnd_out, Xsize, Ysize)){
			Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);            //if out of bounds, PR=0, LLR=0
			Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0
			update_grids();
			return;
		}
		// This is the likelihood ratio for within model update. 
		LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);
		// In the following, we calculate the prior ratio for the intensities.

		//The if-statement is added to see if we need to use background or signal prior.
		//Background Prior calculation
		float BgPDF, BgPDFtest;
		BgPDF = OffgStat.findPDF(E_Active.BG);           //background prior for the current state
		BgPDFtest = OffgStat.findPDF(E_Active_Test.BG);  //background prior for the proposed state

		//if (IPDF != 0 && BgPDF != 0 && IPDFtest != 0 && BgPDFtest != 0)
		//{
		//tmp = tmp*IPDFtest / IPDF;                 //prior ratio for intensity
		//PrI = tmp;

		// The following term is the prior ratio for the offset background parameters and prior ratio for the slopes of background along X and Y axes. The prior for the offset 
		// is provided by the user as a numeric array but the prior for the two other parameters are taken to be normal distribution.
		PrBg = (BgPDFtest / BgPDF)*exp(-(E_Active_Test.BBG*E_Active_Test.BBG - E_Active.BBG*E_Active.BBG) / 2)*exp(-(E_Active_Test.ABG*E_Active_Test.ABG - E_Active.ABG*E_Active.ABG) / 2);

		PR = LLR*PrI*PrBg;                           //posterior ratio

		//If accept jump then Accept =1.
		Accept = (urand() < PR);
		//}
		//else if (IPDFtest == 0 && BgPDFtest == 0)
		//{
		//	Accept = 0;
		//}
		//else
		//{
		//	Accept = 1;
		//}
		if (Accept){
			A_jump++;  //remove and caculate from Chain at end. 
			if (E_Active.Signal.at(nn) == 1)
				A_jumpFg++;
			else
				A_jumpBg++;
			E_Active = E_Active_Test;
			M = M_Test;
			NumM = NumM_Test;
		}
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
		update_grids();
	}
	catch (...){
		//print message to command line and then proceed as failed jump
		mexPrintf("Caught Exception in jump. Jump not accepted\n");
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	}
	//mexPrintf("j9\n");
}

void RJMCMC::split(){
	// split
	int JumpType = 4;   //jump type
	float LLR = 0;      //likelihood ratio
	float PR = 0;       //posterior ratio
	int Accept = 0;     //Acceptance
	float PrI = 1;    //Intensity prior ration
	float PrX = 1;    //Position prior ration
	float Us = 1;    //the term related to the Us
	float Jacob = 1;  //Jacobian
	float PrM = 0;       //Prior on models or number of the particles
	N_split++;
	Emitters E_Active_Test = E_Active;
	try {
		if (E_Active.N < 1) return;
		//Generate split parameters
		float Ua = urand();
		//Pick emitter to split
		int jj = min((int)floor((float)E_Active.N*Ua), E_Active.N - 1);
		float Ijp = E_Active.I.at(jj);
		float mux = E_Active.X.at(jj);
		float muy = E_Active.Y.at(jj);
		int BF = E_Active.Signal.at(jj);
		if (BF == 1)
			N_splitFg++;
		else
			N_splitBg++;
		//generate random vector for u
		float u1 = urand();
		float u2 = Params.Split_std*randn();
		float u3 = Params.Split_std*randn();

		//Change one emitter
		E_Active_Test.N = E_Active.N;
		E_Active_Test.I.at(jj) = u1*Ijp;
		E_Active_Test.X.at(jj) = mux + u2;
		E_Active_Test.Y.at(jj) = muy + u3;
		//float PDFfg, PDFbg;
		//PDFfg = SignalStat.findPDF(E_Active.I.at(jj));
		//PDFbg = BackStat.findPDF(E_Active.I.at(jj));
		//E_Active.Signal.at(jj) = 1;
		//if (PDFbg > PDFfg) E_Active.Signal.at(jj) = 0;

		//PDFfg = BackStat.findPDF((1 - u1)*Ijp);
		//PDFbg = BackStat.findPDF((1 - u1)*Ijp);
		//BF = 1;
		//if (PDFbg > PDFfg) BF = 0;

		//Add emitter
		E_Active_Test.add((1 - u1)*Ijp, mux - u1*u2 / (1 - u1), muy - u1*u3 / (1 - u1), BF);
		//Test
		if (!E_Active_Test.test((int)Params.Bnd_out, Xsize, Ysize)){
			Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0, LLR=0
			Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0
			update_grids();
			return;
		}

		//this is likelihood rtio times the prior for the number of the particles.
		LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);//Prior on the number of particles also included
		float tmp1 = 1, tmp2 = 1;
		//Calculations of the intensity PDF for the current state
		//for (int nn = 0; nn < E_Active.N; nn++)
		//{
		if (E_Active.Signal.at(jj) == 1)
			tmp1 = tmp1*SignalStat.findPDF(E_Active.I.at(jj));
		else
			tmp1 = tmp1*BackStat.findPDF(E_Active.I.at(jj));
		//}

		//Calculations of the intensity PDF for the proposed state
		//for (int nn = 0; nn < E_Active_Test.N; nn++)
		//{
		if (E_Active_Test.Signal.at(jj) == 1)
			tmp2 = tmp2*SignalStat.findPDF(E_Active_Test.I.at(jj));
		else
			tmp2 = tmp2*BackStat.findPDF(E_Active_Test.I.at(jj));

		if (E_Active_Test.Signal.back() == 1)
			tmp2 = tmp2*SignalStat.findPDF(E_Active_Test.I.back());
		else
			tmp2 = tmp2*BackStat.findPDF(E_Active_Test.I.back());
		//}
		//Posterior ratio
		float Norm1, Norm2;
		Norm1 = normpdf(u2, 0, Params.Split_std);
		Norm2 = normpdf(u3, 0, Params.Split_std);
		//if (tmp1 != 0 && tmp2 != 0 && Norm1 != 0 && Norm2 != 0)
		//{
		float NEm = 0;
		for (int nn = 0; nn < E_Active.N; nn++)
			NEm = NEm + E_Active.Signal.at(nn);
		PrM = MeanN / (NEm + 1);
		PrI = tmp2 / tmp1;
		float Area = (float)((Xsize + Params.Bnd_out)*(Ysize + Params.Bnd_out));
		PrX = 1 / Area;
		Us = 1 / (Norm1*Norm2);
		Jacob = Ijp / pow(1 - u1, 2);
		PR = LLR*PrI*PrX*Us*Jacob*(Params.P.at(1) / Params.P.at(2))*PrM;

		//If accept jump then Accept =1.
		Accept = (urand() < PR);
		//}
		//else if (tmp2 == 0)
		//{
		//	Accept = 0;
		//}
		//else
		//{
		//	Accept = 1;
		//}
		if (Accept){

			//float NP = 0;
			//for (int nn = 0; nn < E_Active_Test.N; nn++)
			//	NP = NP + E_Active_Test.Signal[nn];
			//if (NP == 3)

			A_split++;
			if (BF == 1)
				A_splitFg++;
			else
				A_splitBg++;
			E_Active = E_Active_Test;
			M = M_Test;
			NumM = NumM_Test;
		}
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
		update_grids();
	}
	catch (...){
		//print message to command line and then proceed as failed jump
		mexPrintf("Caught Exception in split. Split not accepted\n");
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	}
};

void RJMCMC::merge(){
	//merge
	int JumpType = 5;  //jump type
	float LLR = 0;     //likelihood ratio
	float PR = 0;      //posterior ratio
	int Accept = 0;    //Acceptance
	float PrI = 1;    //Intensity prior ration
	float PrX = 1;    //Position prior ration
	float Us = 1;    //the term related to the Us
	float Jacob = 1;  //Jacobian
	float PrM = 0;   //Prior models or number of the emitters
	N_merge++;

	Emitters E_Active_Test = E_Active;
	try {
		if (E_Active.N < 2){ update_grids(); return; }
		//Generate split parameters

		//Pick emitter to split
		float Ua = urand();
		int jj = min((int)floor((float)E_Active.N*Ua), E_Active.N - 1);
		int ii = E_Active.closest(jj, 3 * Params.PSFsigma);
		if (E_Active.Signal.at(jj) == 1)
			N_mergeFg++;
		else
			N_mergeBg++;
		//check distance
		if (E_Active.distance(ii, jj) > 3 * Params.PSFsigma){ update_grids(); return; }

		//make new particle
		float Ijp = E_Active.I.at(jj) + E_Active.I.at(ii);
		float muxjp = (E_Active.I.at(jj) * E_Active.X.at(jj) + E_Active.I.at(ii) * E_Active.X.at(ii)) / Ijp;
		float muyjp = (E_Active.I.at(jj) * E_Active.Y.at(jj) + E_Active.I.at(ii) * E_Active.Y.at(ii)) / Ijp;
		//calculating us
		float u1 = E_Active.I.at(jj) / Ijp;
		float u2 = E_Active.X.at(jj) - muxjp;
		float u3 = E_Active.Y.at(jj) - muyjp;
		E_Active_Test.I.at(jj) = Ijp;
		E_Active_Test.X.at(jj) = muxjp;
		E_Active_Test.Y.at(jj) = muyjp;
		float PDFfg = SignalStat.findPDF(Ijp);
		float PDFbg = BackStat.findPDF(Ijp);
		E_Active_Test.Signal.at(jj) = 1;
		if (PDFbg > PDFfg)
			E_Active_Test.Signal.at(jj) = 0;

		E_Active_Test.remove(ii);  //remove an old particle
		//Test if the new particle is out of the allowed limits for the position and intensity.
		if (!E_Active_Test.test((int)Params.Bnd_out, Xsize, Ysize)){
			Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0, LLR=0
			Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0
			update_grids();
			return;
		}

		//this is likelihood ratio times the prior on the number of the particles
		LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);
		//In the following, we calculate the prior ratio for the intensities.
		float tmp1 = 1, tmp2 = 1;
		//Calculations of the intensity PDF for the current state
		//for (int nn = 0; nn < E_Active.N; nn++)
		//{
		if (E_Active.Signal.at(jj) == 1)
			tmp1 = tmp1*SignalStat.findPDF(E_Active.I.at(jj));
		else
			tmp1 = tmp1*BackStat.findPDF(E_Active.I.at(jj));

		if (E_Active.Signal.at(ii) == 1)
			tmp1 = tmp1*SignalStat.findPDF(E_Active.I.at(ii));
		else
			tmp1 = tmp1*BackStat.findPDF(E_Active.I.at(ii));
		//}
		//Calculations of the intensity PDF for the proposed state
		//for (int nn = 0; nn < E_Active_Test.N; nn++)
		//{	
		if (ii > jj)
		{
			if (E_Active_Test.Signal.at(jj) == 1)
				tmp2 = tmp2*SignalStat.findPDF(E_Active_Test.I.at(jj));
			else
				tmp2 = tmp2*BackStat.findPDF(E_Active_Test.I.at(jj));
		}
		else
		{
			if (E_Active_Test.Signal.at(jj - 1) == 1)
				tmp2 = tmp2*SignalStat.findPDF(E_Active_Test.I.at(jj - 1));
			else
				tmp2 = tmp2*BackStat.findPDF(E_Active_Test.I.at(jj - 1));
		}
		//}
		float Norm1, Norm2;
		Norm1 = normpdf(u2, 0, Params.Split_std);
		Norm2 = normpdf(u3, 0, Params.Split_std);
		//if (Norm1 != 0 && Norm2 != 0 && tmp1 != 0 && tmp2 != 0)
		//{
		float NEm = 0;
		for (int nn = 0; nn < E_Active.N; nn++)
			NEm = NEm + E_Active.Signal.at(nn);

		PrI = tmp1 / tmp2;
		float Area = (float)((Xsize + Params.Bnd_out)*(Ysize + Params.Bnd_out));
		PrX = 1 / Area;
		Us = 1 / (Norm1*Norm2);
		Jacob = Ijp / pow(1 - u1, 2);
		PrM = MeanN / NEm;

		PR = LLR / (PrI*PrX*Jacob*Us*PrM*Params.P.at(1) / Params.P.at(2));

		//If accept jump then Accept =1.
		Accept = (urand() < PR);
		//}
		//else if (Norm1 == 0)
		//{
		//	Accept = 0;
		//}
		//else
		//{
		//	Accept = 1;
		//}
		//If accept merge
		if (Accept){
			//float NP = 0;
			//for (int nn = 0; nn < E_Active_Test.N; nn++)
			//	NP = NP + E_Active_Test.Signal[nn];
			//if (NP == 2)
			if (E_Active.Signal.at(jj) == 1)
				A_mergeFg++;
			else
				A_mergeBg++;
			A_merge++;
			E_Active = E_Active_Test;
			M = M_Test;
			NumM = NumM_Test;
		};
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
		update_grids();
	}
	catch (...){
		//print message to command line and then proceed as failed jump
		mexPrintf("Caught Exception in merge. Merge not accepted\n");
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	}
};


void RJMCMC::birth(){
	// Birth, posterior ratio of birth does not depend on the priors.
	int JumpType = 6;  //Jump Type
	float LLR = 0;     //likelihood ratio
	float PR = 0;      //posterior ratio
	float PrM;         //number of particles prior ratio;
	int Accept = 0;    //Acceptance
	N_birth++;
	Emitters E_Active_Test = E_Active;
	try {
		//Generate birth position
		int BF;
		float PixProb, x, y;
		find_Birth(&x, &y, &PixProb);
		if (Params.IsBackg == 1)
			BF = kround(urand());
		else
			BF = 1;
		//Add emitter
		float I_new;
		if (BF == 1)
		{
			N_birthFg++;
			I_new = SignalStat.genRand(); //picking a random intensity from the PDFarray.
		}
		else
		{
			N_birthBg++;
			I_new = BackStat.genRand();
		}
		E_Active_Test.add(I_new, (float)x, (float)y, BF);
		//Test if the new particle is out of the allowed limits for the position and intensity.
		if (!E_Active_Test.test((int)Params.Bnd_out, Xsize, Ysize)){
			Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0, LLR=0
			Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);  //if out of bounds, PR=0
			update_grids();
			return;
		}
		//posterior ratio times the prior on the number of particles
		LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);
		float NEm = 0;
		for (int nn = 0; nn < E_Active.N; nn++)
			NEm = NEm + E_Active.Signal.at(nn);
		PrM = MeanN / (NEm + 1);
		PR = LLR*Params.P.at(3) / (Params.P.at(4)* (Xsize + Params.Bnd_out)*(Ysize + Params.Bnd_out)*PixProb)*PrM; //posterior ratio is the same as likelihood ratio. 
		//If accept jump then Accept =1.
		Accept = (urand() < PR);
		//If accept birth
		if (Accept){
			//float NP = 0;
			//for (int nn = 0; nn < E_Active_Test.N; nn++)
			//	NP = NP + E_Active_Test.Signal[nn];
			//if (NP == 3)
			if (BF == 1)
				A_birthFg++;
			else
				A_birthBg++;
			A_birth++;
			E_Active = E_Active_Test;
			M = M_Test;
			NumM = NumM_Test;
		};
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
		update_grids();
	}
	catch (...){
		//print message to command line and then proceed as failed jump
		mexPrintf("Caught Exception in birth. Birth not accepted\n");
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	}
};

void RJMCMC::death(){
	// Death, the posterior ratio of death does not depend on the priors
	int JumpType = 7;  //Jump Type
	float LLR = 0;     //likelihood ratio
	float PR = 0;      //posterior ratio
	float PrM;         //Number of particles prior ratio
	int Accept = 0;    //Acceptance
	N_death++;
	Emitters E_Active_Test = E_Active;
	try {
		if (E_Active.N < 1) {
			return; //If there is no particle to kill return
		}
		//Pick death particle
		float Ua = urand();
		int jj = min((int)floor((float)E_Active.N*Ua), (E_Active.N - 1));

		if (E_Active.Signal.at(jj) == 1)
			N_deathFg++;
		else
			N_deathBg++;
		//Remove emitter
		E_Active_Test.remove(jj);
		LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);
		float NEm = 0;
		for (int nn = 0; nn < E_Active.N; nn++)
			NEm = NEm + E_Active.Signal.at(nn);
		PrM = NEm / MeanN;
		PR = LLR*(Params.P.at(4) / Params.P.at(3))*PrM; //posterior ratio is the same as likelihood ratio
		//If accept death
		Accept = (urand() < PR);
		if (Accept){
			//float NP = 0;
			//for (int nn = 0; nn < E_Active_Test.N; nn++)
			//	NP = NP + E_Active_Test.Signal[nn];
			//if (NP == 2)
			A_death++;
			if (E_Active.Signal.at(jj) == 1)
				A_deathFg++;
			else
				A_deathBg++;
			E_Active = E_Active_Test;
			M = M_Test;
			NumM = NumM_Test;
		};
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
		update_grids();
	}
	catch (...){
		//print message to command line and then proceed as failed jump
		mexPrintf("Caught Exception in death. Death not accepted\n");
		Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
		Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	}
};

void RJMCMC::gSplit(){
	//merge
	int JumpType = 8;  //jump type
	float LLR = 0;     //likelihood ratio
	float PR = 0;      //posterior ratio
	int Accept = 0;    //Acceptance
	float PrI = 1;    //Intensity prior ration
	float PrX = 1;    //Position prior ration
	float Us = 1;    //the term related to the Us
	float Jacob = 1;  //Jacobian
	float PrM = 0;   //Prior models or number of the emitters
	if (E_Active.N < 2) return;
	Emitters E_Active_Test = E_Active;
	//Generate split parameters
	float Ua = urand();
	//Pick an emitter
	int jj = min((int)floor((float)E_Active.N*Ua), E_Active.N - 1);
	//Picking a random set nearby emitter
	std::vector<int> Closest;
	Closest = E_Active.closestMulti(jj, 2 * Params.PSFsigma);
	Closest.push_back(jj);
	int Np = Closest.size();
	if (Np == 1) return;
	float I0 = 0, X0 = 0, Y0 = 0;
	for (int nn = 0; nn < Np; nn++)
	{
		I0 = I0 + E_Active.I.at(Closest.at(nn));
		X0 = X0 + E_Active.X.at(Closest.at(nn));
		Y0 = Y0 + E_Active.Y.at(Closest.at(nn));
	}
	//picking random numbers to compensate for the degrees of freedom.
	float u1 = (float)0.3*randn();
	float u2 = (float)0.3*randn();
	float u3 = (float)(1 / float(Np))*urand();
	//Proposed positions and intensities

	for (int nn = 0; nn < Np; nn++)
	{
		E_Active_Test.X.at(Closest.at(nn)) = (E_Active.X.at(Closest.at(nn)) - u3*(X0 / float(Np) + u1)) / (1 - u3);
		E_Active_Test.Y.at(Closest.at(nn)) = (E_Active.Y.at(Closest.at(nn)) - u3*(Y0 / float(Np) + u2)) / (1 - u3);
		E_Active_Test.I.at(Closest.at(nn)) = (1 - u3)*E_Active.I.at(Closest.at(nn));
	}
	E_Active_Test.add(u3*I0, X0 / float(Np) + u1, Y0 / float(Np) + u2, 1);
	//Likelihood ratio
	LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);
	//Intensity prior ratio
	float tmpCurrent, tmpProposed;
	for (int nn = 0; nn < Np; nn++)
	{
		if (E_Active.Signal.at(Closest.at(nn)) == 1)
		{
			tmpCurrent = SignalStat.findPDF(E_Active.I.at(Closest.at(nn)));
		}
		else
			tmpCurrent = BackStat.findPDF(E_Active.I.at(Closest.at(nn)));

		if (E_Active_Test.Signal.at(nn) == 1)
			tmpProposed = SignalStat.findPDF(E_Active_Test.I.at(Closest.at(nn)));
		else
			tmpProposed = BackStat.findPDF(E_Active_Test.I.at(Closest.at(nn)));

		PrI = PrI * tmpProposed / tmpCurrent;
	}
	if (E_Active_Test.Signal.back() == 1)
		tmpProposed = SignalStat.findPDF(E_Active_Test.I.back());
	else
		tmpProposed = BackStat.findPDF(E_Active_Test.I.back());
	PrI = PrI * (tmpProposed);
	//PDF of the random Us
	float PDFu1, PDFu2, PDFu3;
	PDFu1 = normpdf(u1, 0.0f, 0.4f);
	PDFu2 = normpdf(u2, 0.0f, 0.4f);
	PDFu3 = float(Np);
	Us = 1 / (PDFu1*PDFu2*PDFu3);
	//Number of emitters prior ratio
	float NEm = 0;
	for (int nn = 0; nn < E_Active.N; nn++)
		NEm = NEm + E_Active.Signal.at(nn);
	PrM = MeanN / (NEm + 1);
	//Location prior ratio
	float Area = (float)((Xsize + Params.Bnd_out)*(Ysize + Params.Bnd_out));
	PrX = 1 / Area;
	//Jacobian
	Jacob = I0 / pow(1 - u3, float(Np + 1));
	//Posterior ratio
	PR = LLR*PrI*PrX*Us*Jacob*PrM*(Params.P.at(5) / Params.P.at(6));
	Accept = (urand() < PR);
	if (Accept)
	{
		E_Active = E_Active_Test;
		M = M_Test;
		NumM = NumM_Test;
	}
	Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
	Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	update_grids();
}

void RJMCMC::gMerge(){
	int JumpType = 9;  //jump type
	float LLR = 0;     //likelihood ratio
	float PR = 0;      //posterior ratio
	int Accept = 0;    //Acceptance
	float PrI = 1;    //Intensity prior ration
	float PrX = 1;    //Position prior ration
	float Us = 1;     //the term related to the Us
	float Jacob = 1;  //Jacobian
	float PrM = 0;   //Prior models or number of the emitters
	std::vector<int> Closest;
	if (E_Active.N < 3) return;
	Emitters E_Active_Test = E_Active;
	//Generate split parameters
	float Ua = urand();
	//Pick an emitter
	int jj = min((int)floor((float)E_Active.N*Ua), E_Active.N - 1);
	//Picking a group of nearby emitter
	Closest = E_Active.closestMulti(jj, 3 * Params.PSFsigma);
	int Np = Closest.size();
	if (Np == 0) return;
	float I0 = E_Active.I.at(jj), X0 = E_Active.X.at(jj), Y0 = E_Active.Y.at(jj);
	for (int nn = 0; nn < Np; nn++)
	{
		I0 = I0 + E_Active.I.at(Closest.at(nn));
		X0 = X0 + E_Active.X.at(Closest.at(nn));
		Y0 = Y0 + E_Active.Y.at(Closest.at(nn));
	}
	float U1, U2, U3, XX0 = 0, YY0 = 0;
	U3 = E_Active.I.at(jj) / I0;
	for (int nn = 0; nn < Np; nn++)
	{
		E_Active_Test.I.at(Closest.at(nn)) = E_Active.I.at(Closest.at(nn)) / (1 - U3);
		E_Active_Test.X.at(Closest.at(nn)) = U3*E_Active.X.at(jj) + (1 - U3)*E_Active.X.at(Closest.at(nn));
		XX0 = XX0 + U3*E_Active.X.at(jj) + (1 - U3)*E_Active.X.at(Closest.at(nn));
		E_Active_Test.Y.at(Closest.at(nn)) = U3*E_Active.Y.at(jj) + (1 - U3)*E_Active.Y.at(Closest.at(nn));
		YY0 = YY0 + U3*E_Active.Y.at(jj) + (1 - U3)*E_Active.Y.at(Closest.at(nn));
	}
	U1 = E_Active.X.at(jj) - XX0 / Np;
	U2 = E_Active.Y.at(jj) - YY0 / Np;
	//intensity ratio
	float tmpProposed, tmpCurrent;
	for (int nn = 0; nn < Np; nn++)
	{
		if (E_Active_Test.Signal.at(Closest.at(nn)) == 1.0f)
			tmpProposed = SignalStat.findPDF(E_Active_Test.I.at(Closest.at(nn)));
		else
			tmpProposed = BackStat.findPDF(E_Active_Test.I.at(Closest.at(nn)));

		if (E_Active.Signal.at(Closest.at(nn)) == 1.0f)
			tmpCurrent = SignalStat.findPDF(E_Active.I.at(Closest.at(nn)));
		else
			tmpCurrent = BackStat.findPDF(E_Active.I.at(Closest.at(nn)));

		PrI = PrI * (tmpCurrent / tmpProposed); //Note that this should be the opposite because it appears in the bottom.
	}

	if (E_Active.Signal.at(jj) == 1.0f)
		tmpCurrent = SignalStat.findPDF(E_Active.I.at(jj));
	else
		tmpCurrent = BackStat.findPDF(E_Active.I.at(jj));

	PrI = PrI * tmpCurrent; //Note that this should be the opposite because it appears in the bottom.
	E_Active_Test.remove(jj);
	LLR = likelihood_ratio(E_Active_Test, SampledPSF, Params);

	float PDFu1, PDFu2, PDFu3;
	PDFu1 = normpdf(U1, 0.0f, 0.4f);
	PDFu2 = normpdf(U2, 0.0f, 0.4f);
	PDFu3 = (float)Np;
	Us = 1 / (PDFu1*PDFu2*PDFu3);
	//Number of emitters prior ratio
	float NEm = 0;
	for (int nn = 0; nn < E_Active.N; nn++)
		NEm = NEm + E_Active.Signal.at(nn);
	PrM = NEm / MeanN;
	//Location prior ratio
	float Area = (float)((Xsize + Params.Bnd_out)*(Ysize + Params.Bnd_out));
	PrX = 1 / Area;
	//Jacobian
	Jacob = I0 / pow(1 - U3, float(Np + 1));
	//Posterior ratio
	PR = LLR / (PrI*PrX*Us*Jacob*PrM*(Params.P.at(5) / Params.P.at(6)));
	if (PDFu1 == 0 || PDFu2 == 0)
		PR = -1;

	Accept = (urand() < PR);
	if (Accept)
	{
		E_Active = E_Active_Test;
		M = M_Test;
		NumM = NumM_Test;
	}

	Chain.addEmitters(E_Active, LLR, PR, JumpType, Accept, Params);
	Chain_Test.addEmitters(E_Active_Test, LLR, PR, JumpType, Accept, Params);
	update_grids();
}

void RJMCMC::backFore()
{
	try {
		if (E_Active.N == 0) return;
		if (Params.IsBackg == 0) return;
		float Ua = urand();
		int nn = min((int)floor((float)E_Active.N*Ua), E_Active.N - 1);
		float BackP = BackStat.findPDF(E_Active.I.at(nn));
		float SignalP = SignalStat.findPDF(E_Active.I.at(nn));
		float PrSig;
		if (E_Active.Signal.at(nn) == 1)
		{
			if (SignalP != 0)
			{
				PrSig = BackP / SignalP;
				if (PrSig > urand())
					E_Active.Signal.at(nn) = 0;
			}
			else
			{
				E_Active.Signal.at(nn) = 0;
			}
		}
		else
		{
			if (BackP != 0)
			{
				PrSig = SignalP / BackP;
				if (PrSig > urand())
					E_Active.Signal.at(nn) = 1;
			}
			else
			{
				E_Active.Signal.at(nn) = 1;
			}
		}
	}
	catch (...){
		//print message to command line and then proceed as failed jump
		mexPrintf("Caught Exception in backFore. BackFore not accepted\n");
	}
};
//add the accepted state of the second part of the chain to the output image (PGrid)
void RJMCMC::update_grids(){
	if (IsBurnin) return;                                     //If the first part of the chain return
	int N = E_Active.N;                                         //number of the particles in this state and also the number of iterations for the following loop.
	int x, y;
	float i;
	for (int nn = 0; nn<E_Active.N; nn++){
		x = kround((float)Params.Grid_Zoom*(E_Active.X.at(nn) + Params.DriftX));    //zooming the X-position by a factor of Grid_Zoom
		y = kround((float)Params.Grid_Zoom*(E_Active.Y.at(nn) + Params.DriftY));    //zooming the Y-position by a factor of Grid_Zoom
		i = E_Active.I.at(nn);
		if ((x >= 0) && (x < Xsize*Params.Grid_Zoom) &&               //if the particle position is inside the range of the image which is [0 Zoom*Xsize] and [0 Zoom*Ysize] then add the particle.
			(y >= 0) && (y<Ysize*Params.Grid_Zoom) && (i>Params.Icutoff))
		{
			if (E_Active.Signal.at(nn) == 1)
				PGridV.at(x*Ysize*Params.Grid_Zoom + y) = PGridV.at(x*Ysize*Params.Grid_Zoom + y) + 1;
			else
				PBackV.at(x*Ysize*Params.Grid_Zoom + y) = PBackV.at(x*Ysize*Params.Grid_Zoom + y) + 1;
		}
		;
	}
};
// This function computes the likelihood ratio, which is the same for every mechanism
float RJMCMC::likelihood_ratio(Emitters & E_Test, std::vector<float> SampledPSF, MCMC_Parameters Params){
	if (Params.FirstCall == 1) // make model the first time. 
	{
		switch (OptModel){
		case 0:
			M = Model(E_Active, Xsize, Ysize, Params.PSFsigma, sCMOScorr);
			break;
		case 1:
			NumM = NumModel(E_Active, Xsize, Ysize, Params.PSFsizeV, SampledPSF, sCMOScorr);
			break;
		}
		Params.FirstCall = 0;
	}
	//coputing the model of the proposed state 
	switch (OptModel){
	case 0:
		M_Test = Model(E_Test, Xsize, Ysize, Params.PSFsigma, sCMOScorr);
		break;
	case 1:
		NumM_Test = NumModel(E_Test, Xsize, Ysize, Params.PSFsizeV, SampledPSF, sCMOScorr);
		break;
	}

	int N = Xsize*Ysize; //the whole number of the pixel inside the image.
	float tmp = 0;
	switch (OptModel){
	case 0:{
		for (int nn = 0; nn < N; nn++)
			tmp += M.M.at(nn) - M_Test.M.at(nn) + Data.at(nn) * log(M_Test.M.at(nn) / M.M.at(nn));   //log-likelihood.
		break; }
	case 1:{
		for (int nn = 0; nn < N; nn++)
			tmp += NumM.M.at(nn) - NumM_Test.M.at(nn) + Data.at(nn) * log(NumM_Test.M.at(nn) / NumM.M.at(nn));
		break; }
	}//log-likelihood.
	return exp(tmp);
};
//finding the most likely place for the birth of the new particle.
void RJMCMC::find_Birth(float* x, float*y, float *PixProb){
	int N = Xsize*Ysize;
	std::vector<float> tmp;
	tmp.resize(N);
	int ii;
	if (Params.FirstCall == 1) // make model the first time if necessary 
	{
		switch (OptModel){
		case 0:
			M = Model(E_Active, Xsize, Ysize, Params.PSFsigma, sCMOScorr);
			break;
		case 1:
			NumM = NumModel(E_Active, Xsize, Ysize, Params.PSFsizeV, SampledPSF, sCMOScorr);
			break;
		}
		Params.FirstCall = 0;
	}
	// This samples from the Residum like a 2D Prob. Dist. 
	switch (OptModel){
	case 0:{
		for (int nn = 0; nn < N; nn++) //Residum
			tmp.at(nn) = max(0, Data.at(nn) - M.M.at(nn));
		break; }
	case 1:{
		for (int nn = 0; nn < N; nn++) //Residum
			tmp.at(nn) = max(0, Data.at(nn) - NumM.M.at(nn));
		break; }
	}
	for (int nn = 1; nn<N; nn++) //cumulative dist
		tmp.at(nn) += tmp.at(nn - 1);
	int Ind = -1;
	float r = (tmp.at(N - 1) * urand());
	for (ii = 0; ii < N; ii++)
	{
		Ind = Ind + 1;
		if (r <= tmp.at(ii)) break;
	}
	if (tmp.at(N - 1) != 0 && Ind != 0)
		PixProb[0] = (tmp.at(Ind) - tmp.at(Ind - 1)) / tmp.at(N - 1);
	else if (Ind == 0)
		PixProb[0] = tmp.at(Ind) / tmp.at(N - 1);
	else
		PixProb[0] = 1;

	PixProb[0] = 1;
	y[0] = krem(ii, Ysize)+urand()-0.5f;
	x[0] = (int)floor(((float)ii) / Ysize)+urand()-0.5f;
	return;
};
