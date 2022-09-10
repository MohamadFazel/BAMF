#pragma comment(lib, "kernel32.lib")
#ifdef __unix__
#include <mex.h>
#else
#include <mex.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LibRJMCMCreturn.hpp"

//functions for creating structures
void make_Params(mxArray *InArray, MCMC_Parameters &Params){
	//mexPrintf("P0: %g, P1: %g, P2: %g, P3: %g, P4: %g\n", InArray.P[0], InArray.P[1], InArray.P[2], InArray.P[3], InArray.P[4]);
	Params.PSFsigma = (float)mxGetScalar(mxGetField(InArray, 0, "PSFsigma"));
	Params.N_Trials = (int)mxGetScalar(mxGetField(InArray, 0, "N_Trials"));
	Params.N_Burnin = (int)mxGetScalar(mxGetField(InArray, 0, "N_Burnin"));
	Params.I_stdFg = (float)mxGetScalar(mxGetField(InArray, 0, "I_stdFg"));
	Params.I_stdBg = (float)mxGetScalar(mxGetField(InArray, 0, "I_stdBg"));
	Params.X_stdFg = (float)mxGetScalar(mxGetField(InArray, 0, "X_stdFg"));
	Params.X_stdBg = (float)mxGetScalar(mxGetField(InArray, 0, "X_stdBg"));
	Params.P_Burnin = (float *)mxGetData(mxGetField(InArray, 0, "P_Burnin"));
	Params.P_Trials = (float *)mxGetData(mxGetField(InArray, 0, "P_Trials"));
	Params.Rho = (float)mxGetScalar(mxGetField(InArray, 0, "Rho"));
	Params.Bnd_out = (float)mxGetScalar(mxGetField(InArray, 0, "Bnd_out"));
	Params.Grid_Zoom = (int)mxGetScalar(mxGetField(InArray, 0, "Grid_Zoom"));
	Params.Split_std = (float)mxGetScalar(mxGetField(InArray, 0, "Split_std"));
	Params.Icutoff = (float)mxGetScalar(mxGetField(InArray, 0, "Icutoff"));
	Params.BG_std = (float)mxGetScalar(mxGetField(InArray, 0, "BG_std"));
	Params.ABBG_std = (float)mxGetScalar(mxGetField(InArray, 0, "ABBG_std"));
	Params.DX = (float)mxGetScalar(mxGetField(InArray, 0, "DX"));
	Params.DriftX = (float)mxGetScalar(mxGetField(InArray, 0, "DriftX"));
	Params.DriftY = (float)mxGetScalar(mxGetField(InArray, 0, "DriftY"));
	Params.PSFsize = (int *)mxGetData(mxGetField(InArray, 0, "PSFsize"));
	Params.IsBackg = (int)mxGetScalar(mxGetField(InArray, 0, "IsBackg"));
	for (int ii = 0; ii < 8; ii++)
	{
		float PBtemp = 0, PTtemp = 0;
		for (int mm = 0; mm <= ii; mm++)
		{
			PBtemp = PBtemp + Params.P_Burnin[mm];
			PTtemp = PTtemp + Params.P_Trials[mm];
		}
		Params.PV_Burnin.push_back(PBtemp);
		Params.PV_Trials.push_back(PTtemp);
	}
	for (int ii = 0; ii < 4; ii++)
		Params.PSFsizeV.push_back(Params.PSFsize[ii]);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//mexPrintf("0.05\n");
	//mexEvalString("drawnow;");
	mxArray *mxMCMC, *mxE_active;
	int Xsize, Ysize, OptModel;
	float * Data, *PGrid, *PBack, *OffSetPDF, *sCMOScorr, MeanN, *SignalPDF, *BackPDF, *SampledPSF;
	MCMC_Parameters Params;
	mwSize ndims = 2;
	mwSize dims[2];
	//mexPrintf("0.1\n");
	/* Assign pointers to each input. */
	Data = (float*)mxGetData(prhs[0]);
	mxE_active = (mxArray *)(prhs[1]);
	mxMCMC = (mxArray *)(prhs[2]);
	Xsize = (int)mxGetScalar(prhs[3]);
	Ysize = (int)mxGetScalar(prhs[4]);
	OffSetPDF = (float*)mxGetData(prhs[5]);
	SignalPDF = (float*)mxGetData(prhs[6]);
	BackPDF = (float*)mxGetData(prhs[7]);
	sCMOScorr = (float*)mxGetData(prhs[8]);
	SampledPSF = (float*)mxGetData(prhs[9]);
	OptModel = (int)mxGetScalar(prhs[10]); //0 means using the analytical normal function as model and 1 means use the input (4D) PSF 
	//mexPrintf("0.2\n");
	int DimSize0, DimSize1, LengthPDF, SignalPDFlength, BackPDFlength, SignalPDFSize0, SignalPDFSize1, BackPDFSize0, BackPDFSize1;
	DimSize0 = (int)mxGetM(prhs[5]);
	DimSize1 = (int)mxGetN(prhs[5]);
	SignalPDFSize0 = (int)mxGetM(prhs[6]);
	SignalPDFSize1 = (int)mxGetN(prhs[6]);
	BackPDFSize0 = (int)mxGetM(prhs[7]);
	BackPDFSize1 = (int)mxGetN(prhs[7]);
	//mexPrintf("DimSize0: %d, DimSize1: %d, DimSize2: %d \n",DimSize[0],DimSize[1],DimSize[2]);
	LengthPDF = DimSize1;
	if (LengthPDF<DimSize0)
	{
		LengthPDF = DimSize0;
	}
	SignalPDFlength = SignalPDFSize0;
	if (SignalPDFlength < SignalPDFSize1)
	{
		SignalPDFlength = SignalPDFSize1;
	}
	BackPDFlength = BackPDFSize0;
	if (BackPDFlength < BackPDFSize1)
	{
		BackPDFlength = BackPDFSize1;
	}
	//mexPrintf("0.3\n");
	//mexPrintf("LengthPDF:%d, SignalPDFlength:%d, BackPDFlength:%d",LengthPDF,SignalPDFlength,BackPDFlength);
	//mexPrintf("0\n");
	//mexEvalString("drawnow;");
	make_Params(mxMCMC, Params);
	//mexPrintf("P0: %g, P1: %g, P2: %g, P3: %g, P4: %g\n", Params.PJump, Params.PSplit, Params.PMerge, Params.PBirth, Params.PDeath);
	//mexPrintf("P0: %g, P1: %g, P2: %g, P3: %g, P4: %g\n", Params.P[0], Params.P[1], Params.P[2], Params.P[3], Params.P[4]);
	//mexPrintf("Sigma: %g \n",Params.PSFsigma);
	//mexPrintf("0.4\n");
	//create PGRid output variable-----------------------------------------
	// This is the PGrid. It is written to directly by RJMCMC.
	dims[0] = (int)(Ysize)*Params.Grid_Zoom;
	dims[1] = (int)(Xsize)*Params.Grid_Zoom;
	plhs[0] = mxCreateNumericArray(ndims, dims, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(ndims, dims, mxSINGLE_CLASS, mxREAL);
	PGrid = (float *)mxGetData(plhs[0]);
	PBack = (float *)mxGetData(plhs[1]);
	MeanN = Params.Rho*Xsize*Ysize;
	//--------------------------------------------------------------
	//mexPrintf("0.5\n");
	//make Active Emitters. This gets setup and copied to RJMCMC
	int N = (int)mxGetScalar(mxGetField(mxE_active, 0, "N"));
	//int *Signal;
	//for (int nn = 0; nn < N; nn++)   //???????????????????
	//	Signal[nn] = 1;
	Emitters E_Active = Emitters(N,
		(float *)mxGetData(mxGetField(mxE_active, 0, "I")),
		(float *)mxGetData(mxGetField(mxE_active, 0, "X")),
		(float *)mxGetData(mxGetField(mxE_active, 0, "Y")),
		(float)mxGetScalar(mxGetField(mxE_active, 0, "BG")), 0, 0
		);
	//mexPrintf("0.6\n");
	// Note that if the offset background is zero then we will have a division by zero in the likelihood ratio.
	if (E_Active.BG == 0)
		E_Active.BG = 1;
	//make Stat object. This gets setup and copied to RJMCMC
	Stat OffgStat = Stat(OffSetPDF, Params.DX, LengthPDF);
	Stat SignalStat = Stat(SignalPDF, Params.DX, SignalPDFlength);
	Stat BackStat = Stat(BackPDF, Params.DX, BackPDFlength);
	//mexPrintf("0.7\n");
	RJMCMC RJ = RJMCMC(Data, Params, Xsize, Ysize, sCMOScorr, MeanN, SampledPSF, OptModel);
	//mexPrintf("0.8\n");
	RJ.IsBurnin = 0;
	RJ.E_Active = E_Active;
	RJ.OffgStat = OffgStat;
	RJ.BackStat = BackStat;
	RJ.SignalStat = SignalStat;
	//mexPrintf("1\n");
	//mexEvalString("drawnow;");
	RJ.start_Chain();
	//mexPrintf("2\n");
	//mexEvalString("drawnow;");
	int Iter = Xsize*Ysize*Params.Grid_Zoom*Params.Grid_Zoom;
	for (int ii = 0; ii < Iter; ii++)
	{
		PGrid[ii] = RJ.PGridV.at(ii);
		PBack[ii] = RJ.PBackV.at(ii);
	}
	// These are other outputs ---------------------------------
	//make a Chain and Chain_Test structure and export the vectors .
	mwSize NFields = (mwSize)RJ.Chain.N.size();
	int NLoop = (int)RJ.Chain.N.size();
	int nfields = 12; //N,X,Y,Photons,LLR,PR,JumpType,Accepted,Signal
	const char *fieldnames[12] =       /* pointers to field names */
	{ "N", "X", "Y", "Photons", "Signal", "LLR", "PR", "JumpType", "Accepted", "BG", "ABG", "BBG" };

	//mexPrintf("3\n");

	//These are arrays of structures with the above fields
	plhs[2] = mxCreateStructMatrix(NFields, 1, nfields, fieldnames);
	plhs[3] = mxCreateStructMatrix(NFields, 1, nfields, fieldnames);

	//loop over NTrials. This is for Chain. ----------------------------
	//mexPrintf("4");
	//mexEvalString("drawnow;");
	for (int nn = 0; nn < NLoop - 1; nn++)
	{
		double * pData, *sData;
		//Get emitters from the chain at trial nn

		Emitters E = RJ.Chain.getEmitters(nn);

		const mwSize NN = E.N;

		//make and write fields to the structure
		mxArray *outS = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		sData = (double*)mxGetData(outS);
		sData[0] = (double)E.N;
		if (E.I[0] == 0) sData[0] = 0;
		mxSetFieldByNumber(plhs[2], nn, 0, outS);  // 0 is field number ( for N)

		mxArray *out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++)
			pData[ii] = E.X[ii];
		mxSetFieldByNumber(plhs[2], nn, 1, out);

		out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.Y[ii];
		mxSetFieldByNumber(plhs[2], nn, 2, out);

		out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.I[ii];
		mxSetFieldByNumber(plhs[2], nn, 3, out);

		out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.Signal[ii];
		mxSetFieldByNumber(plhs[2], nn, 4, out);

		mxSetFieldByNumber(plhs[2], nn, 5, mxCreateDoubleScalar((double)RJ.Chain.LLR[nn]));
		mxSetFieldByNumber(plhs[2], nn, 6, mxCreateDoubleScalar((double)RJ.Chain.PR[nn]));
		//mxSetFieldByNumber(plhs[1], nn, 6, mxCreateDoubleScalar((double)RJ.Chain.PrI[nn]));
		//mxSetFieldByNumber(plhs[1], nn, 7, mxCreateDoubleScalar((double)RJ.Chain.PrX[nn]));
		//mxSetFieldByNumber(plhs[1], nn, 8, mxCreateDoubleScalar((double)RJ.Chain.Us[nn]));
		//mxSetFieldByNumber(plhs[1], nn, 9, mxCreateDoubleScalar((double)RJ.Chain.Jacob[nn]));
		mxSetFieldByNumber(plhs[2], nn, 7, mxCreateDoubleScalar((double)RJ.Chain.JumpType[nn]));
		mxSetFieldByNumber(plhs[2], nn, 8, mxCreateDoubleScalar((double)RJ.Chain.Accepted[nn]));
		mxSetFieldByNumber(plhs[2], nn, 9, mxCreateDoubleScalar((double)E.BG));
		mxSetFieldByNumber(plhs[2], nn, 10, mxCreateDoubleScalar((double)E.ABG));
		mxSetFieldByNumber(plhs[2], nn, 11, mxCreateDoubleScalar((double)E.BBG));
	}
	//mexPrintf("4\n");
	//End Chain ----------------------------------------------
	//mexPrintf("5");
	//mexEvalString("drawnow;");
	NLoop = (int)RJ.Chain_Test.N.size();
	//loop over NTrials. This is for Chain_Test. ----------------------------
	for (int nn = 0; nn<NLoop - 1; nn++)
	{
		//mexPrintf("nn:%d",nn);
		double * pData;
		//Get emitters from the chain at trial nn
		Emitters E = RJ.Chain_Test.getEmitters(nn);

		const mwSize NN = E.N;
		double Number;
		if (E.I[0] == 0)
			Number = 0;
		else
			Number = (double)E.N;
		//make and write fields to the structure
		mxSetFieldByNumber(plhs[3], nn, 0, mxCreateDoubleScalar(Number));  // 0 is field number ( for N)

		mxArray *out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.X[ii];
		mxSetFieldByNumber(plhs[3], nn, 1, out);

		out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.Y[ii];
		mxSetFieldByNumber(plhs[3], nn, 2, out);

		out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.I[ii];
		mxSetFieldByNumber(plhs[3], nn, 3, out);

		out = mxCreateNumericMatrix(1, NN, mxDOUBLE_CLASS, mxREAL);
		pData = (double*)mxGetData(out);
		for (int ii = 0; ii<E.N; ii++) pData[ii] = E.Signal[ii];
		mxSetFieldByNumber(plhs[3], nn, 4, out);

		mxSetFieldByNumber(plhs[3], nn, 5, mxCreateDoubleScalar((double)RJ.Chain.LLR[nn]));
		mxSetFieldByNumber(plhs[3], nn, 6, mxCreateDoubleScalar((double)RJ.Chain.PR[nn]));
		//mxSetFieldByNumber(plhs[2], nn, 6, mxCreateDoubleScalar((double)RJ.Chain.PrI[nn]));
		//mxSetFieldByNumber(plhs[2], nn, 7, mxCreateDoubleScalar((double)RJ.Chain.PrX[nn]));
		//mxSetFieldByNumber(plhs[2], nn, 8, mxCreateDoubleScalar((double)RJ.Chain.Us[nn]));
		//mxSetFieldByNumber(plhs[2], nn, 9, mxCreateDoubleScalar((double)RJ.Chain.Jacob[nn]));
		mxSetFieldByNumber(plhs[3], nn, 7, mxCreateDoubleScalar((double)RJ.Chain.JumpType[nn]));
		mxSetFieldByNumber(plhs[3], nn, 8, mxCreateDoubleScalar((double)RJ.Chain.Accepted[nn]));
		mxSetFieldByNumber(plhs[3], nn, 9, mxCreateDoubleScalar((double)E.BG));
		mxSetFieldByNumber(plhs[3], nn, 10, mxCreateDoubleScalar((double)E.ABG));
		mxSetFieldByNumber(plhs[3], nn, 11, mxCreateDoubleScalar((double)E.BBG));
	}
	//mexPrintf("5\n");
	//End Chain_Test ----------------------------------------------
	//mexPrintf("6");
	//mexEvalString("drawnow;");
	//make a Statistics structure and export 
	nfields = 16;
	const char *fieldnamesStat[16] =       /* pointers to field names */
	{ "NTrials", "Accept_Jump", "Accept_Split", "Accept_Merge", "Accept_Birth", "Accept_Death", "Accept_JumpFg", "Accept_SplitFg", "Accept_MergeFg",
	"Accept_BirthFg", "Accept_DeathFg", "Accept_JumpBg", "Accept_SplitBg", "Accept_MergeBg", "Accept_BirthBg", "Accept_DeathBg" };

	//These are arrays of structures with the above fields
	plhs[4] = mxCreateStructMatrix(1, 1, nfields, fieldnamesStat);
	mxSetFieldByNumber(plhs[4], 0, 0, mxCreateDoubleScalar((double)Params.N_Trials + Params.N_Burnin));
	//mxSetFieldByNumber(plhs[4], 0, 1, mxCreateDoubleScalar((double)RJ.Chain.getAcceptance_Jump()));
	mxSetFieldByNumber(plhs[4], 0, 1, mxCreateDoubleScalar((double)RJ.A_jump / RJ.N_jump));
	mxSetFieldByNumber(plhs[4], 0, 2, mxCreateDoubleScalar((double)RJ.A_split / RJ.N_split));
	mxSetFieldByNumber(plhs[4], 0, 3, mxCreateDoubleScalar((double)RJ.A_merge / RJ.N_merge));
	mxSetFieldByNumber(plhs[4], 0, 4, mxCreateDoubleScalar((double)RJ.A_birth / RJ.N_birth));
	mxSetFieldByNumber(plhs[4], 0, 5, mxCreateDoubleScalar((double)RJ.A_death / RJ.N_death));
	mxSetFieldByNumber(plhs[4], 0, 6, mxCreateDoubleScalar((double)RJ.A_jumpFg / RJ.N_jumpFg));
	mxSetFieldByNumber(plhs[4], 0, 7, mxCreateDoubleScalar((double)RJ.A_splitFg / RJ.N_splitFg));
	mxSetFieldByNumber(plhs[4], 0, 8, mxCreateDoubleScalar((double)RJ.A_mergeFg / RJ.N_mergeFg));
	mxSetFieldByNumber(plhs[4], 0, 9, mxCreateDoubleScalar((double)RJ.A_birthFg / RJ.N_birthFg));
	mxSetFieldByNumber(plhs[4], 0, 10, mxCreateDoubleScalar((double)RJ.A_deathFg / RJ.N_deathFg));
	mxSetFieldByNumber(plhs[4], 0, 11, mxCreateDoubleScalar((double)RJ.A_jumpBg / RJ.N_jumpBg));
	mxSetFieldByNumber(plhs[4], 0, 12, mxCreateDoubleScalar((double)RJ.A_splitBg / RJ.N_splitBg));
	mxSetFieldByNumber(plhs[4], 0, 13, mxCreateDoubleScalar((double)RJ.A_mergeBg / RJ.N_mergeBg));
	mxSetFieldByNumber(plhs[4], 0, 14, mxCreateDoubleScalar((double)RJ.A_birthBg / RJ.N_birthBg));
	mxSetFieldByNumber(plhs[4], 0, 15, mxCreateDoubleScalar((double)RJ.A_deathBg / RJ.N_deathBg));
	//mexPrintf("6\n");
	//mexEvalString("drawnow;");
}