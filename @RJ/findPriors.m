function [obj,SRImage]=findPriors(obj,PlotFlag,Y_N_BgPrior)
%findPrior() finds intensity priors of signal and background. 
%   This function uses single molecule fitting to find distribution of the
%   signal emitter intensities, the offset background and the background 
%   emitter intensities.
%
% INPUTS:
%   obj:            The class object containing the following parameters
%                   that are used here.
%       Data:           Gain and offset corrected data.
%       RJStruct:       Structure containing input parameters to the core c++
%                       code of RJMCMC. 
%   PlotFlag:       Whether to plot the found distributions or not. (1,0)
%   Camera_Var:     Scaled variance image use in SCMOS fitting (Optional)
%
% OUTPUTS:
%   obj:            The class object containing the following outouts.
%      SignalPDF:      Distribution of Photons, spaced by RJStruct.P_dP
%      OffsetPDF:      Distribution of the offset background.
%      BgPDF:          Distribution of the offset Background.
%
% REQUIRES:
%   Statistic Toolbox
%   Parallel Processing Toolbox
%   NVIDIA GPU
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel and Keith A. Lidke (Lidke Lab 2018)
%
Boxsize = 8;
Sigma1 = obj.RJStruct.PSF_Sigma;
Sigma2 = 2*Sigma1;
Minval=obj.SMAStruct.MeanPhotons/(4*pi*Sigma1^2)/4;
[~,ROIStack]=RJ.findROI(obj.Data,Sigma1,Sigma2,Boxsize,Minval);
SMD=RJ.gaussMLE(ROIStack,'Basic','CCD',Sigma1);
SMR=RJ.threhsSMA(SMD,obj.SMAStruct,Boxsize);

Photon = SMR.Photons;
if length(Photon)>200
    ID = quantile(Photon,0.98);
    Photon(Photon>ID)=[];
end
PdEm = fitdist(double(Photon),'Kernel','BandWidth',75);
Xarray = 0:obj.RJStruct.P_dP:max(Photon)+1000*obj.RJStruct.P_dP;
SignalPDF = pdf(PdEm,Xarray);

GrSignalPDF = gradient(SignalPDF);
Ind = find([0,diff(sign(GrSignalPDF))]);
IndMax = find(max(SignalPDF(Ind))==SignalPDF);
IndRemove = find(Ind > IndMax);
if length(IndRemove) >= 2
    Photon(Photon>Ind(IndRemove(1))*obj.RJStruct.P_dP)=[];
    PdEm = fitdist(double(Photon),'Kernel','BandWidth',75);
    Xarray = 0:obj.RJStruct.P_dP:max(Photon)+1000*obj.RJStruct.P_dP;
    SignalPDF = pdf(PdEm,Xarray);
end

%SignalPDF = gampdf(Xarray,PdEm.a,PdEm.b);
%SignalPDF = normpdf(Xarray,2000,500);
%SignalPDF = (normpdf(Xarray,2000,100)+unifpdf(Xarray,800,2000))/2;
PdBg = fitdist(double(SMR.Bg),'Gamma');
OffsetPDF = gampdf(Xarray,PdBg.a,PdBg.b);
if strcmp('Yes',Y_N_BgPrior)
    %MeanBg = PdEm.a*PdEm.b/obj.RJStruct.BgPriorRange; 
    MeanBg = mean(double(Photon))/obj.RJStruct.BgPriorRange;
    BgPDF = exppdf(Xarray,MeanBg);
else
    BgPDF = zeros(size(Xarray),'single');
end

%PlotFlag = 0;
if isunix
  PlotFlag = 0;
end
%BgPDF = (normpdf(Xarray,500,50)+unifpdf(Xarray,50,900)+cat(2,(0:0.0002:0.0018),zeros(1,length(Xarray)-10)))/2.055;
if PlotFlag
   %dipshow(SRImage);colormap('hot')
   figure;
   subplot(2,2,1)
   plot(Xarray(1:end-500),SignalPDF(1:end-500),'linewidth',2)
   hold(gca,'on'); histogram(Photon,'NumBins',100,'Normalization','pdf')
   ylabel('PDF')
   xlabel('Photons')
   title('Signal Prior')
   subplot(2,2,2)
   plot(Xarray(1:end-500),BgPDF(1:end-500),'linewidth',2)
   ylabel('PDF')
   xlabel('Photons')
   title('Background Prior')
   subplot(2,2,3)
   plot(Xarray(1:100),OffsetPDF(1:100),'linewidth',2)
   hold(gca,'on'); histogram(SMR.Bg,'NumBins',20,'Normalization','pdf')
   ylabel('PDF')
   xlabel('Photons')
   title('Offset Prior')
   subplot(2,2,4)
   plot(Xarray(1:end-500),SignalPDF(1:end-500),'linewidth',1.2)
   hold(gca,'on');plot(Xarray(1:end-500),BgPDF(1:end-500),'linewidth',1.2)
   %plot(Xarray(1:end-500),OffsetPDF(1:end-500),'linewidth',1.2)
   ylabel('PDFs')
   xlabel('Photons')
   legend('Signal','Background')
end
obj.RJStruct.P_EmPhotons = SignalPDF;
obj.RJStruct.P_BgEmPhotons = BgPDF;
obj.RJStruct.P_Offset = OffsetPDF;
% PEmHandle=histogram(SR.SMR.Photons,'BinEdges',(0:RJStruct.P_dP:max(SR.SMR.Photons)+300*RJStruct.P_dP),'Normalization','pdf');
% PEm = PEmHandle.Values;
% PoffsetHandle=histogram(SR.SMR.Bg,'BinEdges',(0:RJStruct.P_dP:max(SR.SMR.Photons)+200*RJStruct.P_dP),'Normalization','pdf');
% PBg = PoffsetHandle.Values;
end

