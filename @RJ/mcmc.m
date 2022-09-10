function [FoundMAP,SMDkmeans,SMDbg,AcceptChain,Stat]=mcmc(Data,RJstruct,Xsize,Ysize,DriftY,DriftX,sCMOSVar,SampledPSF_in,OptModel_in,ROI_Struct,SMD,SMDbg,ROIid,FrameID,Xstd_cutoff,Istd_cutoff)
%mcmc() gets a single ROI and use mcmc to process the input ROI. 
%   This function use the MAPN from RJMCMC to initialize an MCMC chain.
%   The MCMC chain then are used to calculated the MAPN results from MCMC. 
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
%The particles in the overlaping regions have to be removed.
ROI_Ind = ROI_Struct.ROI_Ind(ROIid); 
Xsub = ROI_Struct.Xsub;
Ysub = ROI_Struct.Ysub;
XsizeROI = ROI_Struct.Xsize(ROIid);
YsizeROI = ROI_Struct.Ysize(ROIid);
Overlap = ROI_Struct.Overlap;
Ind = ROI_Ind;
ix = floor((Ind-1)/Xsub) + 1;
iy = Ind - (ix-1)*Xsub;
BeginROIx = 0;
BeginROIy = 0;
EndROIx = XsizeROI;
EndROIy = YsizeROI;
BeginOverlapX = 0;
BeginOverlapY = 0;
if ix ~= Xsub %ceil(length(obj.Data')/obj.SubIm) %if this is not the last block in the row
   EndROIx = EndROIx-Overlap;
end
if iy ~= Ysub 
   EndROIy = EndROIy-Overlap;
end
if ix ~= 1 
   BeginROIx = BeginROIx+Overlap;
   BeginOverlapX = Overlap;
end
if  iy ~= 1 
    BeginROIy = BeginROIy+Overlap;
    BeginOverlapY = Overlap;
end

Ind1x_kmeans = find(SMD.X-RJstruct.XDrift(FrameID)<=BeginROIx);
Ind2x_kmeans = find(SMD.X-RJstruct.XDrift(FrameID)>EndROIx);
Ind1y_kmeans = find(SMD.Y-RJstruct.YDrift(FrameID)<=BeginROIy);
Ind2y_kmeans = find(SMD.Y-RJstruct.YDrift(FrameID)>EndROIy);
Index_kmeans = unique([Ind1x_kmeans,Ind2x_kmeans,Ind1y_kmeans,Ind2y_kmeans]);

Ind1x_BGkmeans = find(SMDbg.X-RJstruct.XDrift(FrameID)<=BeginROIx);
Ind2x_BGkmeans = find(SMDbg.X-RJstruct.XDrift(FrameID)>EndROIx);
Ind1y_BGkmeans = find(SMDbg.Y-RJstruct.YDrift(FrameID)<=BeginROIy);
Ind2y_BGkmeans = find(SMDbg.Y-RJstruct.YDrift(FrameID)>EndROIy);
Index_BGkmeans = unique([Ind1x_BGkmeans,Ind2x_BGkmeans,Ind1y_BGkmeans,Ind2y_BGkmeans]);

E_active.N = int32(length(SMD.X)); %the number of particles to start with
E_active.I = single(SMD.I); % the intensities to start with
E_active.X = single(SMD.X-DriftX); %the X-positions to start with
E_active.Y = single(SMD.Y-DriftY); %the Y-positions to start with
E_active.BG = single(SMD.Bg); %the uniform offset background

MCMC.PSFsigma = single(RJstruct.PSF_Sigma); %the standard deviation of the PSF, which is model as a Gaussisn.
MCMC.Grid_Zoom = single(RJstruct.SRZoom); %Zoom
MCMC.N_Trials = single(RJstruct.N_TrialsMCMC); %The number of trials.
MCMC.N_Burnin = single(RJstruct.N_BurninMCMC); %The number of burn in trials.

%calculate some jump sized based in intensity prior
Ptmp=RJstruct.P_Offset;
POffsetwidth=sum(Ptmp>(max(Ptmp)/2))*RJstruct.P_dP;
Ptmp=RJstruct.P_BgEmPhotons;
PBgwidth=sum(Ptmp>(max(Ptmp)/2))*RJstruct.P_dP;
Ptmp=RJstruct.P_EmPhotons;
PEmwidth=sum(Ptmp>(max(Ptmp)/2))*RJstruct.P_dP;
P = [1 0 0 0 0 0 0 0];
MCMC.I_stdFg = RJstruct.IstdFg;%PEmwidth/1000; %The jumps in the intensity are taken from a gaussian dist with this standard deviation.
MCMC.I_stdBg = RJstruct.IstdFg/1.5;%PBgwidth/100;
Ind = find(Ptmp==max(Ptmp));
MCMC.X_stdFg = RJstruct.XstdFg;%10*RJstruct.PSF_Sigma/sqrt(Ind(1)*RJstruct.P_dP); %The jumps in the position are picked from a Normal dist. with this standard deviation.
MCMC.X_stdBg = RJstruct.XstdFg/1.5;%100*RJstruct.PSF_Sigma/sqrt(Ind(1)*RJstruct.P_dP);
MCMC.P_Burnin = single(P); %The propbabilities for proposing different jumps. (insideJump, split, merge, birth, death)
MCMC.P_Trials = single(P);
MCMC.Split_std = 0.5; 
MCMC.Bnd_out = 2; %The range of area where we are allowed to propose a particle outside the box.
MCMC.Icutoff = single(RJstruct.MinPhoton);
MCMC.Rho = single(RJstruct.Rho); %The average number of the particles per pixel.
MCMC.BG_std = RJstruct.BGstd;%POffsetwidth/90;
MCMC.ABBG_std = 0.02;
MCMC.DX = single(RJstruct.P_dP);
MCMC.DriftX = single(DriftX);
MCMC.DriftY = single(DriftY);
MCMC.IsBackg = 0;
if nargin < 7 || isscalar(sCMOSVar)
  sCMOSVar = zeros(size(Data),'single');
end      
if nargin < 8 || isscalar(SampledPSF_in) || isempty(SampledPSF_in)
   SampledPSF_in = zeros(10,10,2,2); 
   OptModel_in = 0;
end
MCMC.PSFsize = uint32(size(SampledPSF_in));
SignalPDF = single(RJstruct.P_EmPhotons);
BackgPDF = single(RJstruct.P_BgEmPhotons);
OffSetPDF = single(RJstruct.P_Offset);
ROIySize = int32(Ysize);
ROIxSize = int32(Xsize);
Data = Data(1:Ysize,1:Xsize);
sCMOSVar = sCMOSVar(1:Ysize,1:Xsize);
SampledPSF = single(SampledPSF_in);
OptModel = uint32(OptModel_in);
%dipshow(Data)
%fprintf('MXSize:%g, MYsize: %g \n',size(Data,2),size(Data,1));
%fprintf('Data:%g, E_active:%g, MCMC:%g, ROIySize:%g, ROIxSize:%g, OffSetPDF:%g, SignalPDF:%g, BackgPDF:%g, sCMOSVar:%g\n',...
 %   numel(Data),numel(fieldnames(E_active)),numel(fieldnames(MCMC)),numel(ROIySize),numel(ROIxSize),numel(OffSetPDF),numel(SignalPDF),numel(BackgPDF),numel(sCMOSVar));
[PImg,PBg,AcceptChain,~,Stat]=cRJMCMC(Data,E_active,MCMC,ROIxSize,ROIySize,OffSetPDF,SignalPDF,BackgPDF,sCMOSVar,SampledPSF,OptModel);
X = [];
Y = [];
I = [];
ChainLength = length(AcceptChain);
for nn = round(ChainLength/2):ChainLength
    x = single(AcceptChain(nn).X);
    y = single(AcceptChain(nn).Y);
    i = single(AcceptChain(nn).Photons);
    X = cat(1,X,x);
    Y = cat(1,Y,y);
    I = cat(1,I,i);
end
FoundMAP.X = mean(X);
FoundMAP.Y = mean(Y);
FoundMAP.I = mean(I);
FoundMAP.Stat=Stat;
FoundMAP.XChain = X;
FoundMAP.YChain = Y;
FoundMAP.IChain = I;
for nn = 1:size(X,2)
    FoundMAP.X_SE(nn) = std(X(:,nn));
    FoundMAP.Y_SE(nn) = std(Y(:,nn));
    FoundMAP.I_SE(nn) = std(I(:,nn));
end
FoundMAP.ROIcorner = [(ix-1)*ROI_Struct.SizeData(2)/Xsub,(iy-1)*ROI_Struct.SizeData(1)/Ysub];

Ind1x = find(FoundMAP.X<=BeginROIx);
Ind2x = find(FoundMAP.X>EndROIx);
Ind1y = find(FoundMAP.Y<=BeginROIy);
Ind2y = find(FoundMAP.Y>EndROIy);
Index = unique([Ind1x,Ind2x,Ind1y,Ind2y]);

FoundMAP.X = FoundMAP.X+(ix-1)*(ROI_Struct.Xsize(1)-Overlap)-BeginOverlapX;
FoundMAP.Y = FoundMAP.Y+(iy-1)*(ROI_Struct.Xsize(1)-Overlap)-BeginOverlapY;

SMDkmeans.X = SMD.X+(ix-1)*(ROI_Struct.Xsize(1)-Overlap)-BeginOverlapX;
SMDkmeans.Y = SMD.Y+(iy-1)*(ROI_Struct.Xsize(1)-Overlap)-BeginOverlapY;
SMDkmeans.I = SMD.I;
SMDkmeans.X_SE = SMD.X_SE;
SMDkmeans.Y_SE = SMD.Y_SE;
SMDkmeans.I_SE = SMD.I_SE;
SMDkmeans.X(Index_kmeans)=[];
SMDkmeans.Y(Index_kmeans)=[];
SMDkmeans.I(Index_kmeans)=[];
SMDkmeans.X_SE(Index_kmeans)=[];
SMDkmeans.Y_SE(Index_kmeans)=[];
SMDkmeans.I_SE(Index_kmeans)=[];

SMDbg.X = SMDbg.X+(ix-1)*(ROI_Struct.Xsize(1)-Overlap)-BeginOverlapX;
SMDbg.Y = SMDbg.Y+(iy-1)*(ROI_Struct.Xsize(1)-Overlap)-BeginOverlapY;
SMDbg.I = SMDbg.I;
SMDbg.X_SE = SMDbg.X_SE;
SMDbg.Y_SE = SMDbg.Y_SE;
SMDbg.I_SE = SMDbg.I_SE;
SMDbg.X(Index_BGkmeans)=[];
SMDbg.Y(Index_BGkmeans)=[];
SMDbg.I(Index_BGkmeans)=[];
SMDbg.X_SE(Index_BGkmeans)=[];
SMDbg.Y_SE(Index_BGkmeans)=[];
SMDbg.I_SE(Index_BGkmeans)=[];
%removing the particles in the overlaping regions
FoundMAP.X(Index) = [];
FoundMAP.Y(Index) = [];
FoundMAP.I(Index) = [];
FoundMAP.X_SE(Index) = [];
FoundMAP.Y_SE(Index) = [];
FoundMAP.I_SE(Index) = [];
if nargin > 14
    FoundMAP.X(FoundMAP.X_SE>Xstd_cutoff) = [];
    FoundMAP.Y(FoundMAP.X_SE>Xstd_cutoff) = [];
    FoundMAP.I(FoundMAP.X_SE>Xstd_cutoff) = [];
    FoundMAP.Y_SE(FoundMAP.X_SE>Xstd_cutoff) = [];
    FoundMAP.I_SE(FoundMAP.X_SE>Xstd_cutoff) = [];
    FoundMAP.X_SE(FoundMAP.X_SE>Xstd_cutoff) = [];
    
    FoundMAP.X(FoundMAP.Y_SE>Xstd_cutoff) = [];
    FoundMAP.Y(FoundMAP.Y_SE>Xstd_cutoff) = [];
    FoundMAP.I(FoundMAP.Y_SE>Xstd_cutoff) = [];
    FoundMAP.I_SE(FoundMAP.Y_SE>Xstd_cutoff) = [];
    FoundMAP.X_SE(FoundMAP.Y_SE>Xstd_cutoff) = [];
    FoundMAP.Y_SE(FoundMAP.Y_SE>Xstd_cutoff) = [];

    FoundMAP.X(FoundMAP.I_SE>Istd_cutoff) = [];
    FoundMAP.Y(FoundMAP.I_SE>Istd_cutoff) = [];
    FoundMAP.I(FoundMAP.I_SE>Istd_cutoff) = [];
    FoundMAP.X_SE(FoundMAP.I_SE>Istd_cutoff) = [];
    FoundMAP.Y_SE(FoundMAP.I_SE>Istd_cutoff) = [];
    FoundMAP.I_SE(FoundMAP.I_SE>Istd_cutoff) = [];
end
end
