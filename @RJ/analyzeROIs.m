function [PIm,PBg,ClustInfo,ChainCell,StatCell,PChainCell,TimeSeriesIm,TimeSeriesBg]=...
    analyzeROIs(ROI_Struct,RJstrc,ROI_Data,ROI_Var,SampledPSF_Cell,OptModel,Keep,TimeSeries,LoopNum)
%analyseROIs() is called inside analyzeData() and sends ROIs to rjmcmc().
%   There are 2 for-loops that iterates over frames and ROIs in each frame
%   and send them to rjmcmc() to be processed. The results are then saved
%   inside a structure called ClustInfo.
%
% INPUTS:
%   ROI_Struct:    Structure containing information about ROIs (sub-regions)
%   RJstrc:        Structure containing input parameters to the core c++
%                  code of RJMCMC.
%   ROI_Data:      Sub-regions or ROIs of data.
%   ROI_Var:       Correponding read-out noise to the sub-regions.
%   SampledPSF_Cell:Cell containing found numerical PSFs for different 
%                  regions of the field of view.   
%   OptModel:      This parameter could be either 0 or 1. 0 means use 
%                  Gaussian PSF. 1 means use found numerical PSF.
%   Keep:          0 does not keep any chains, 1 keeps the chain of the
%                  accepted jumps, 2 keeps the both chains of accepted and
%                  proposed jumps.
%   TimeSeries:    0 means do not keep the posterior of the ROIs and 1
%                  means save them.
%   LoopNum:       Loop number in the parfor-loop.
%                  
% OUTPUTS:
%   PIm:           Posterior image of the Signal
%   PBg:           Posterior image of the structure background.
%   ClustInfo:     The structure containing information about emitters.
%   ChainCell:     Chain of the accepted jumps for each ROI.
%   StatCell:      Statistics of the jumps fro each ROI.
%   PChainCell:    Chain of the proposed jumps for each ROI.
%   TimeSeriesIm:  If the user want to save the posterior image of each
%                  subregion, they will be saved in this parameter.
%   TimeSeriesBg:  The posterior image for structured background for the 
%                  ROIs.
% REQUIRES:
%   MATLAB 2014a or higher versions.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%


Zs = ROI_Struct.SizeData(3);
StackSize = size(ROI_Data);
StackSize(3) = size(ROI_Data,3);
Iter_pp = max(ROI_Struct.FrameNum(:));
Iter_mm = max(ROI_Struct.ROI_Ind(:));
YsizeVec = ROI_Struct.Ysize;
XsizeVec = ROI_Struct.Xsize;          
if isempty(Keep)
    Keep = 0; 
end
PIm = zeros(StackSize(1)*RJstrc.SRZoom,StackSize(2)*RJstrc.SRZoom,Iter_mm,'single');
%PBg = zeros(StackSize(1)*RJstrc.SRZoom,StackSize(2)*RJstrc.SRZoom,Iter_mm,'single');
PBg = [];
if TimeSeries == 1
    TimeSeriesIm = cell(1,Zs);
    TimeSeriesBg = cell(1,Zs);
else
    TimeSeriesIm = [];
    TimeSeriesBg = [];
end
if isscalar(ROI_Var)
   ROI_Var = zeros(size(ROI_Data),'single');
end
ClustInfo.X = single([]);
ClustInfo.Y = single([]);
ClustInfo.I = single([]);
ClustInfo.X_SE = single([]);
ClustInfo.Y_SE = single([]);
ClustInfo.I_SE = single([]);
ClustInfo.X_kmeans = single([]);
ClustInfo.Y_kmeans = single([]);
ClustInfo.I_kmeans = single([]);
ClustInfo.X_SE_kmeans = single([]);
ClustInfo.Y_SE_kmeans = single([]);
ClustInfo.I_SE_kmeans = single([]);
ClustInfo.Xbg = single([]);
ClustInfo.Ybg = single([]);
ClustInfo.Ibg = single([]);
ClustInfo.Xbg_SE = single([]);
ClustInfo.Ybg_SE = single([]);
ClustInfo.Ibg_SE = single([]);
ClustInfo.FrameNum = single([]);
ClustInfo.FrameNumBg = single([]);
ClustInfo.FrameNum_kmeans = single([]);
ClustInfo.ROInum = single([]);
ClustInfo.Bg = single([]);
ClustInfo.ABg = single([]);
ClustInfo.BBg = single([]);
ClustInfo.JumpXRJ = single([]);
ClustInfo.JumpIRJ = single([]);
ClustInfo.JumpBgRJ = single([]);
ClustInfo.JumpXFg = single([]);
ClustInfo.JumpIFg = single([]);
ClustInfo.JumpXBg = single([]);
ClustInfo.JumpIBg = single([]);
ClustInfo.SplitRJ = single([]);
ClustInfo.MergeRJ = single([]);
ClustInfo.BirthRJ = single([]);
ClustInfo.DeathRJ = single([]);
ClustInfo.JumpMC = single([]);
ClustInfo.ROIcornerX = single([]);
ClustInfo.ROIcornerY = single([]);
if Keep == 1
    ChainCell = cell(1,StackSize(3));
    StatCell = cell(1,StackSize(3));
    PChainCell = [];
elseif Keep == 2
    ChainCell = cell(1,StackSize(3));
    PChainCell = cell(1,StackSize(3));
    StatCell = cell(1,StackSize(3));
else
    ChainCell = [];
    PChainCell = [];
    StatCell = [];
end
%make loops and call 

for pp = 1:Iter_pp
  
    if round(pp/50)==pp/50
       fprintf('frame:%g\n',pp);
    end
    if isempty(RJstrc.YDrift)
        DriftY = 0;
    else
        DriftY = RJstrc.YDrift(pp+(LoopNum-1)*ROI_Struct.TempFrame);
    end
    if isempty(RJstrc.XDrift)
        DriftX = 0;
    else
        DriftX = RJstrc.XDrift(pp+(LoopNum-1)*ROI_Struct.TempFrame);
    end
    for mm = 1:Iter_mm
        if OptModel == 1
            SampledPSF = SampledPSF_Cell{mm};
        else
            SampledPSF = [];
        end
        nn = (pp-1)*Iter_mm+mm;
        ThisSubRegion = ROI_Data(:,:,nn);
        ThisSubVar = ROI_Var(:,:,nn);
        Ysize = YsizeVec(nn);
        Xsize = XsizeVec(nn);
        [ROIChain,Statistics,PImS,PBgS,ROIPChain]=RJ.rjmcmc(ThisSubRegion,RJstrc,Ysize,Xsize,DriftY,DriftX,ThisSubVar,SampledPSF,OptModel); %setup and call c-codef
        PIm(:,:,mm) = PIm(:,:,mm)+PImS; 
        %PBg(:,:,mm) = PBg(:,:,mm)+PBgS;  
        %This is new findMAPN method
        [SMD,MostFrequent] = RJ.findMap(ROIChain,1);%INputs need frame, dataset number
        [SMDbg,~] = RJ.findMap(ROIChain,0);
        ClustInfo.Bg = cat(2,ClustInfo.Bg,SMD.Bg);
        ClustInfo.ABg = cat(2,ClustInfo.ABg,SMD.ABg);
        ClustInfo.BBg = cat(2,ClustInfo.BBg,SMD.BBg);
        if MostFrequent~=0 
            [FoundMAP,SMDkmeans,SMDbg,~,Stat]=RJ.mcmc(ThisSubRegion,RJstrc,Xsize,Ysize,DriftY,DriftX,ThisSubVar,SampledPSF,OptModel,ROI_Struct,SMD,SMDbg,mm,pp);
            ClustInfo.X = cat(2,ClustInfo.X,single(FoundMAP.X));
            ClustInfo.Y = cat(2,ClustInfo.Y,single(FoundMAP.Y));
            ClustInfo.I = cat(2,ClustInfo.I,single(FoundMAP.I));
            ClustInfo.X_SE = cat(2,ClustInfo.X_SE,single(FoundMAP.X_SE));
            ClustInfo.Y_SE = cat(2,ClustInfo.Y_SE,single(FoundMAP.Y_SE));
            ClustInfo.I_SE = cat(2,ClustInfo.I_SE,single(FoundMAP.I_SE));
            ClustInfo.X_kmeans = cat(2,ClustInfo.X_kmeans,single(SMDkmeans.X));
            ClustInfo.Y_kmeans = cat(2,ClustInfo.Y_kmeans,single(SMDkmeans.Y));
            ClustInfo.I_kmeans = cat(2,ClustInfo.I_kmeans,single(SMDkmeans.I));
            ClustInfo.X_SE_kmeans = cat(2,ClustInfo.X_SE_kmeans,single(SMDkmeans.X_SE));
            ClustInfo.Y_SE_kmeans = cat(2,ClustInfo.Y_SE_kmeans,single(SMDkmeans.Y_SE));
            ClustInfo.I_SE_kmeans = cat(2,ClustInfo.I_SE_kmeans,single(SMDkmeans.I_SE));
            ClustInfo.FrameNum_kmeans = cat(2,ClustInfo.FrameNum_kmeans,(pp+(LoopNum-1)*ROI_Struct.TempFrame)*ones(size(SMDkmeans.X),'single'));
            ClustInfo.Xbg = cat(2,ClustInfo.Xbg,single(SMDbg.X));
            ClustInfo.Ybg = cat(2,ClustInfo.Ybg,single(SMDbg.Y));
            ClustInfo.Ibg = cat(2,ClustInfo.Ibg,single(SMDbg.I));
            ClustInfo.Xbg_SE = cat(2,ClustInfo.Xbg_SE,single(SMDbg.X_SE));
            ClustInfo.Ybg_SE = cat(2,ClustInfo.Ybg_SE,single(SMDbg.Y_SE));
            ClustInfo.Ibg_SE = cat(2,ClustInfo.Ibg_SE,single(SMDbg.I_SE));
            ClustInfo.JumpXRJ = cat(2,ClustInfo.JumpXRJ,single(Statistics.Accept_Jump));
            ClustInfo.JumpIRJ = cat(2,ClustInfo.JumpIRJ,single(Statistics.Accept_Jump));
            ClustInfo.JumpBgRJ = cat(2,ClustInfo.JumpBgRJ,single(Statistics.Accept_JumpBg));
            ClustInfo.JumpXFg = cat(2,ClustInfo.JumpXFg,single(Statistics.Accept_JumpFg));
            ClustInfo.JumpIFg = cat(2,ClustInfo.JumpIFg,single(Statistics.Accept_JumpFg));
            ClustInfo.JumpXBg = cat(2,ClustInfo.JumpXBg,single(Statistics.Accept_JumpBg));
            ClustInfo.JumpIBg = cat(2,ClustInfo.JumpIBg,single(Statistics.Accept_JumpBg));
            ClustInfo.SplitRJ = cat(2,ClustInfo.SplitRJ,single(Statistics.Accept_SplitFg));
            ClustInfo.MergeRJ = cat(2,ClustInfo.MergeRJ,single(Statistics.Accept_MergeFg));
            ClustInfo.BirthRJ = cat(2,ClustInfo.BirthRJ,single(Statistics.Accept_BirthFg));
            ClustInfo.DeathRJ = cat(2,ClustInfo.DeathRJ,single(Statistics.Accept_DeathFg));
            %ClustInfo.JumpBgRJ = cat(2,ClustInfo.JumpBgRJ,single(Statistics.Accept_JumpBg));
            ClustInfo.JumpMC = cat(2,ClustInfo.JumpMC,single(Stat.Accept_JumpFg));
            ClustInfo.FrameNum = cat(2,ClustInfo.FrameNum,(pp+(LoopNum-1)*ROI_Struct.TempFrame)*ones(size(FoundMAP.X),'single'));
            ClustInfo.FrameNumBg = cat(2,ClustInfo.FrameNumBg,(pp+(LoopNum-1)*ROI_Struct.TempFrame)*ones(size(SMDbg.X),'single'));
            ClustInfo.ROInum = cat(2,ClustInfo.ROInum,LoopNum*ones(size(FoundMAP.X),'single'));
            ClustInfo.ROIcornerX = cat(2,ClustInfo.ROIcornerX,FoundMAP.ROIcorner(1)*ones(size(FoundMAP.X),'single'));
            ClustInfo.ROIcornerY = cat(2,ClustInfo.ROIcornerY,FoundMAP.ROIcorner(2)*ones(size(FoundMAP.X),'single'));
        end
        %keep some statistics by calling some other function:
        if Keep == 1
            ChainCell{nn} = ROIChain;
            StatCell{nn} = Statistics; 
        elseif Keep == 2
            ChainCell{nn} = ROIChain;
            PChainCell{nn} = ROIPChain;
            StatCell{nn} = Statistics; 
        end
        if TimeSeries 
            TimeSeriesIm{nn}=sparse(double(PImS));
            TimeSeriesBg{nn}=sparse(double(PBgS));
        end
   end

end
end 
       