function [Chain,Stat,PImage,PBackg,ProposedChain]=rjmcmc(Data,RJstruct,Ysize,Xsize,DriftY,DriftX,sCMOSVar,SampledPSF_in,OptModel_in)
%rjmcmc() gets a single ROI and use rjmcmc to process the input ROI. 
%   The coere RJMCMC algorithm is implemented here. This function gets the
%   ROIs and returns the chain and the image of the chain for the signal 
%   and background emitters.
%
% INPUTS:
%   Data:     The raw data which is the ROI.
%   RJStruct: Structure containing info on Data, such as priors for signal,
%             background, offset, PSF-Sigma etc.
%   Ysize:    The size of the ROI along Y-axis.
%   Xsize:    The size of the ROI along X-axis.
%   sCMOSVar: The scaled variance of sCMOS camera corresponding to the
%             input ROI.
%
% OUTPUTS: 
%   Chain: Structure containing the following fields    
%     X:        X-positions of the found emitters.
%     Y:        Y-positions of the found emitters.
%     N:        Number of the found emitters.
%     Signal:   Clssification of the emitters into signal and background.
%     LLR:      Likelihood ratio.
%     PR:       Prior ratio.
%     JumpType: The type of the jump made in each state. 
%     BG:       The found constant for the offset background.
%     ABG:      The slope of the offset background along the X-axis.
%     BBG:      The slope of the offset background along the Y-axis.
%
% REQUIRES:
%   MATLAB 2014 or higher versions.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%
if size(Data,3)>1
    error('The input data to rjmcmc() must be a single ROI.');
end
E_active.N = int32(0); %the number of particles to start with
E_active.I = single([]); % the intensities to start with
E_active.X = single([]); %the X-positions to start with
E_active.Y = single([]); %the Y-positions to start with
E_active.BG = single(RJstruct.P_dP*(find(max(RJstruct.P_Offset)==RJstruct.P_Offset)-1)); %the uniform offset background

MCMC.PSFsigma = single(RJstruct.PSF_Sigma); %the standard deviation of the PSF, which is model as a Gaussisn.
MCMC.Grid_Zoom = int32(RJstruct.SRZoom); %Zoom
MCMC.N_Trials = int32(RJstruct.N_Trials); %The number of trials.
MCMC.N_Burnin = int32(RJstruct.N_Burnin); %The number of burn in trials.

%calculate some jump sized based in intensity prior
Ptmp=RJstruct.P_Offset;
POffsetwidth=sum(Ptmp>(max(Ptmp)/2))*RJstruct.P_dP;
Ptmp=RJstruct.P_BgEmPhotons;
PBgwidth=sum(Ptmp>(max(Ptmp)/2))*RJstruct.P_dP;
Ptmp=RJstruct.P_EmPhotons;
PEmwidth=sum(Ptmp>(max(Ptmp)/2))*RJstruct.P_dP;

MCMC.I_stdFg = single(RJstruct.IstdFg);%The jumps in the intensity are taken from a gaussian dist with this standard deviation.
MCMC.I_stdBg = single(RJstruct.IstdFg/2);
Ind = find(Ptmp==max(Ptmp));
MCMC.X_stdFg = single(RJstruct.XstdFg);%The jumps in the position are picked from a Normal dist. with this standard deviation.
MCMC.X_stdBg = single(RJstruct.XstdFg/2);
MCMC.P_Burnin = single(RJstruct.P_Burnin); %The propbabilities for proposing different jumps. (insideJump, split, merge, birth, death)
MCMC.P_Trials = single(RJstruct.P_Trials);
MCMC.Split_std = single(0.5); 
MCMC.Bnd_out = single(2); %The range of area where we are allowed to propose a particle outside the box.
MCMC.Icutoff = single(RJstruct.MinPhoton);
MCMC.Rho = single(RJstruct.Rho); %The average number of the particles per pixel.
MCMC.BG_std = single(RJstruct.BGstd);
MCMC.ABBG_std = single(0.02);
MCMC.DX = single(RJstruct.P_dP);
MCMC.DriftX = single(DriftX);
MCMC.DriftY = single(DriftY);
MCMC.IsBackg = int32(RJstruct.IsBackg);
if nargin < 7 || isscalar(sCMOSVar)
  sCMOSVar = zeros(size(Data),'single');
end      
if nargin < 8 || isscalar(SampledPSF_in) || isempty(SampledPSF_in)
   SampledPSF_in = zeros(10,10,2,2); 
   OptModel_in = 0;
end
MCMC.PSFsize = int32(size(SampledPSF_in));
SignalPDF = single(RJstruct.P_EmPhotons);
BackgPDF = single(RJstruct.P_BgEmPhotons);
OffSetPDF = single(RJstruct.P_Offset);
ROIySize = int32(Ysize);
ROIxSize = int32(Xsize);
PImage = zeros(RJstruct.SRZoom*size(Data),'single');
PBackg = zeros(RJstruct.SRZoom*size(Data),'single');
Data = single(Data(1:Ysize,1:Xsize));
sCMOSVar = single(sCMOSVar(1:Ysize,1:Xsize));
if isempty(OptModel_in)
   OptModel_in = 0; 
   SampledPSF_in = zeros(10,10,2,2); 
end
SampledPSF = single(SampledPSF_in);
OptModel = int32(OptModel_in);
[PImg,PBg,AcceptChain,ProposedChain,Stat]=cRJMCMC(Data,E_active,MCMC,ROIxSize,ROIySize,OffSetPDF,SignalPDF,BackgPDF,sCMOSVar,SampledPSF,OptModel);
Chain = AcceptChain;
PImage(1:RJstruct.SRZoom*Ysize,1:RJstruct.SRZoom*Xsize) = PImg;
PBackg(1:RJstruct.SRZoom*Ysize,1:RJstruct.SRZoom*Xsize)=PBg;
end