function [ RJStruct ] = genRJStruct()
%genRJStruct Generates Property Structure for RJMCMC
% RJStruct is a structure containing all the details necessary to run the
% core rjmcmc algorithm on a single sub region. It is an input to the main 
% analysis function RJ.rjmcmc(). 
%   
% The structure has the following fields:
%
% RJStruct:              
%   SRZoom:         Posterior Image Zoom (int > 1) (Default=10) 
%   P_EmPhotons:    Prior PDF of emitters (Nx1)
%   P_BgEmPhotons:  Prior PDF of emitters modeling background (Nx1)
%   P_Offset:       Prior PDF of parameters for Bg offset (NxM)
%   P_dP:           Step size of Prior PDF (Photons)
%   PSF_Sigma:      Sigma of 2D Gaussian PSF 
%   N_Trials:       Number of jumps in RJMCMC chain (Default=3000)
%   N_Burnin:       Number of jumps in RJMCMC burnin (Default=3000)
%   P_Jump:         Probability of jump type (see below)
%   EmDensity:      Expected Emitter Density (Emitter/Pixel)
%   MinPhoton:      Min photons to contribute to posterior image
%
%   P_Jump describes jump types of:
%       [Move, Birth, Death,Split,Merge,Classification]
%   These elements must sum to one. 
%
% INPUTS:
%   None
%
% OUTPUTS:
%   SMD:    An SMD structure with all fields set to []
%
% REQUIRES:
%   MATLAB 2016 or higher versions.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%
RJStruct.SRZoom=20;
RJStruct.P_EmPhotons={}; 
RJStruct.P_BgEmPhotons={};
RJStruct.P_Offset={};
RJStruct.PSF_Sigma=1.3;
RJStruct.N_Trials = 2000;
RJStruct.N_Burnin = 3000;
RJStruct.N_TrialsMCMC = 1000;
RJStruct.N_BurninMCMC = 1000;
RJStruct.P_Burnin=[0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
RJStruct.P_Trials=[0.4,0,0,0.05,0.05,0.15,0.15,0.2];
RJStruct.P_dP=10; 
RJStruct.Rho=.01;
RJStruct.MinPhoton=0;
RJStruct.YDrift={};
RJStruct.XDrift={};
RJStruct.XstdFg = 0.1;
RJStruct.IstdFg = 10;
RJStruct.BGstd = 1;
RJStruct.IsBackg = 1;
RJStruct.BgPriorRange = 12;
end
