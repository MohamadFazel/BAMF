function [Model,Fg,Bg]=makeModel(obj,ROINum,PlotFlag,PSFSigma)
%makeModel() makes an overlay of the data with the Model of different states. 
%   This is a method for diagnosis porpuses. It takes a ROI and
%   its corresponding chain and make a 3D stack, where each slice of the
%   stack is model corresponding to a state of the chain. If PlotFlag=1, it
%   also returns a plot of the overlay of the data with the model.
%
% Inputs:
%   ROI:        The raw ROI.
%   Chain:      The chain corresponding to the given ROI. 
%   PSFsigma:   The PSF size (PSF-sigma) of the emitters. This is an
%               optional input. (default = 1.3 pixels).
%   PlotFlag:   0 no plot and 1 shows the overlay of data and model
%
% Outputs:
%   None.
%
% RQUIRES:
%   MATLAB 2014 or higher versions.
%   DIPImage (www.diplib.org)
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%

Chain = obj.Chain{ROINum};
ROI = obj.ROI_Data(:,:,ROINum);
if nargin > 3
    PSF = PSFSigma;
else
   PSF = obj.RJStruct.PSF_Sigma; 
end
OptModel = obj.FindPSF;
if isempty(OptModel)
    OptModel = 0; 
end
Iter = length(Chain)-1;
[Ys,Xs]=size(ROI);
Model = zeros(Ys,Xs,Iter,'single');
Bg = Model;
Fg = Model;
fprintf('Number of the frames: %g\n',Iter);
for nn = 1:Iter
   SubChain = Chain(nn);
   X = SubChain.X;
   Y = SubChain.Y;
   I = SubChain.Photons;
   Signal = SubChain.Signal;
   [Xg,Yg]=meshgrid((0.5:Xs-0.5),(0.5:Ys-0.5));
   ImBg = zeros(size(Model(:,:,1)));
   ImFg = zeros(size(Model(:,:,1)));
   SubIter = length(X);
   Backg = SubChain.BG;
   XSlope = SubChain.ABG;
   YSlope = SubChain.BBG;
   for mm = 1:SubIter
       switch OptModel
           case 0
               if ~isscalar(PSF)
                   error('PSF_Signma must be a scalar.');
               end
               if Signal(mm)==1
                   ImFg = ImFg + I(mm)*normpdf(Yg,Y(mm),PSF).*normpdf(Xg,X(mm),PSF);
               elseif Signal(mm)==0
                   ImBg = ImBg + I(mm)*normpdf(Yg,Y(mm),PSF).*normpdf(Xg,X(mm),PSF);
               end
           case 1
               if isscalar(PSF)
                   error('PSF_Signma must be a 4D array.');
               end
               if Signal(mm)==1
                   ImFg = ImFg + I(mm)*genModel(single(PSF),single(size(PSF)),single(size(PSF,4)),single(size(ROI,1)),single(Y(mm)),single(X(mm)));
               elseif Signal(mm)==0
                   ImBg = ImBg + I(mm)*genModel(single(PSF),single(size(PSF)),single(size(PSF,4)),single(size(ROI,1)),single(Y(mm)),single(X(mm)));
               end
       end
   end
   Fg(:,:,nn) = ImFg;
   Bg(:,:,nn) = ImBg + XSlope*(Xg-0.5)+YSlope*(Yg-0.5)+Backg;
   Model(:,:,nn)=ImBg+ImFg+XSlope*(Xg-0.5)+YSlope*(Yg-0.5)+Backg;
end
if PlotFlag
    Data = single(repmat(ROI,[1,1,Iter]));
    ImRGB=joinchannels('RGB',Data,Model);
    dipshow(ImRGB)
end
end
