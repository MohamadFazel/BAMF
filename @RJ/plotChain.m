function [ImOf]=plotChain(obj,ROINum,N,PlotFlag,SigBg)
%plotChain displays the MAPN image from RJMCMC chain.
% plotChain() is a method for diagnosis purposes. It takes the Chain and
% size of the ROI plus the Model, which is determined by the number of the
% emitters, and the zoom factor and returns the plot of the found positions
% for the specified model.
%
% INPUTS:
%   Chain:     The Chain corresponding to the ROI.
%   N:         Number of the emitters that specifies model. For N=0, it
%              plots the whole chain.
%   ImSize:    The size of the ROI. The 1st and 2nd elements are
%              respectively the size of ROI along Y and X axes.
%   Zoom:      The magnification factor.
%
% OUTPUTS:
%   ImOf:      The image of the given model.
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
if nargin < 4
   PlotFlag = 1; 
end
if nargin < 5
   SigBg = 1; 
end
Chain = obj.Chain{ROINum};
Zoom = obj.RJStruct.SRZoom;
ImSize(1)=size(obj.ROI_Data,1);
ImSize(2)=size(obj.ROI_Data,1);

Iter = length(Chain);
ImOf = zeros(Zoom*ImSize,'single');
for ii = round(Iter/3):Iter
    SubChain = Chain(ii);
    Signal = SubChain.Signal;
    if N==0
        Y = round(Zoom*(SubChain.Y));
        X = round(Zoom*(SubChain.X));
        XYLength = length(X);
        for nn = 1:XYLength
            if X(nn)>0 && X(nn)<=Zoom*ImSize(2) && Y(nn)>0 && Y(nn)<=Zoom*ImSize(1) && Signal(nn)==SigBg
                   ImOf(X(nn),Y(nn)) = ImOf(X(nn),Y(nn))+1;
            end
        end
    else    
        if N==sum(Signal)
            Y = round(Zoom*(SubChain.Y));
            X = round(Zoom*(SubChain.X));
            XYLength = length(X);
            for nn = 1:XYLength
                if X(nn)>0 && X(nn)<=Zoom*ImSize(2) && Y(nn)>0 && Y(nn)<=Zoom*ImSize(1) && Signal(nn)==SigBg
                   ImOf(X(nn),Y(nn)) = ImOf(X(nn),Y(nn))+1;
                end
            end
        end
    end
end
if PlotFlag
    dipshow(ImOf)
    colormap('hot')
end
end