function obj=findPSF_SMA(obj,ROI_Struct)
%findPSF_SMA() finds the isolated emitters and sends them to findPSF().
%   This function fits the data and uses the results to identify the
%   isolated emitters and then sends them to another function to find the
%   PSF.
%
% INPUTS:
%   obj:        The class structure.
%   ROI_Struct: Structure containing information about the ROIs.
%
% OUTPUTS:
%   obj:        The class structure containing the found PSF.
% 
% REQUIRES: 
%   MATLAB 2014a or higher versions.
%   NVIDIA GPU.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%
MeanPhotons = obj.SMAStruct.MeanPhotons;
BoxSize = obj.SMAStruct.BoxSize;
PSFSigma = obj.RJStruct.PSF_Sigma;
DataROI = obj.ROI_Data;
MinInt=MeanPhotons/(4*pi*PSFSigma^2)/4; %Minimum of the photons
SMDin = [];
PSF = cell(1,ROI_Struct.Xsub*ROI_Struct.Ysub);
for nn = 1:ROI_Struct.Xsub*ROI_Struct.Ysub
  
    ROIStack = [];
    Data = DataROI(1:ROI_Struct.Ysize(nn),1:ROI_Struct.Xsize(nn),nn:ROI_Struct.Xsub*ROI_Struct.Ysub:end);
    try
    [SR,ROIStack]=SMA_Core.findROI(SMDin,single(Data),PSFSigma,2*PSFSigma,BoxSize,MinInt);
    catch
       fprintf('cuda error: Size of Data: %g, PSFSigma: %g, BoxSize: %g, MinInt: %g\n',size(Data,3),PSFSigma,BoxSize,MinInt); 
    end
    if size(ROIStack,3)>30
    [Res,~]=SMA_Core.gaussMLE(single(ROIStack),'Sigma','CCD',PSFSigma);
    Ind1 = Res.Sigma < (BoxSize-4)/2 & Res.Photons < 5*MeanPhotons & ...
        Res.Photons > MeanPhotons/2 & abs(Res.Photons_SE) < MeanPhotons/5 ...
        & Res.LogLikelihood > -(12-10/exp(abs(mean(Res.Photons(Res.Photons_SE<MeanPhotons/5))-1000)/8000))*BoxSize^2 ...
        & Res.Sigma_SE < 1 & Res.X > BoxSize/4 ...
        & Res.X < BoxSize - BoxSize/4 & Res.Y > BoxSize/4 & Res.Y < BoxSize-BoxSize/4; 
    Ind = Ind1 & Res.LogLikelihood>=mean(Res.LogLikelihood(Ind1));
    else
       fprintf('There is not enough input for gaussMLE.\n'); 
       Ind = 1;
    end
    if sum(Ind) < 30
        fprintf('There is not enough single emitters to calcuate PSF. Using just a normal PSF.\n'); 
        [Xg,Yg]=meshgrid((0.5:BoxSize-0.5),(0.5:BoxSize-0.5));
        PSF1 = MeanPhotons*normpdf(Xg,BoxSize/2,PSFSigma).*normpdf(Yg,BoxSize/2,PSFSigma);
        ROIs = repmat(PSF1,[1,1,50]);
        PSF_Sigma = PSFSigma;
    else
        Ind = Ind & Res.Photons > mean(Res.Photons(Ind))/2;
        PSF_Sigma=mean(Res.Sigma(Ind));
        if sum(Ind)>250
            ROIs = ROIStack(:,:,Ind);
            ROIs = ROIs(:,:,1:250);
        else
            ROIs = ROIStack(:,:,Ind);
        end
    end
    PSF{nn}=RJ.findPSF(ROIs,1,0,PSF_Sigma,BoxSize);
end
obj.SampledPSF = PSF;
end