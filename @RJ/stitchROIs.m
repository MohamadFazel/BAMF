function [FinalIm]=stitchROIs(ROIs,ROI_Struct,Overlap,ImSize,BoxSize,Zoom)
%stitchROIs() stiches the chains from all the ROIs to make the final reconstruction.
% 
%
% INPUTS:
%   ROIs:       The input images of the chains from different ROIs.
%   ROI_Struct: The info on the ROIs like their sizes or their indexes.
%   Overlap:    The size of the overlapping regions.
%   ImSize:     A vector containing the size of the frame.
%   BoxSize:    The size of the ROIs without the overlapping regions.
%   Zoom:       The magnification factor.
%
% OUTPUTS:
%   FinalIm:    The final reconstruction.
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
FinalIm = zeros(Zoom*ImSize(1),Zoom*ImSize(2),ImSize(3),'single');
if isempty(ROIs)
    return;
end
ROI_Ind = ROI_Struct.ROI_Ind; 
Xsub = ROI_Struct.Xsub;
Ysub = ROI_Struct.Ysub;
Xsize = ROI_Struct.Xsize*Zoom;
Ysize = ROI_Struct.Ysize*Zoom;
Iters = size(ROIs,3);
for nn = 1:Iters
    frame=ROI_Struct.FrameNum(nn);
    Ind = ROI_Ind(nn);
    ix = floor((Ind-1)/Ysub) + 1;
    iy = Ind - (ix-1)*Ysub;
    
    sy = BoxSize*(iy-1)-(Overlap-1);
    ey = BoxSize*iy+Overlap;
    sx = BoxSize*(ix-1)-(Overlap-1);
    ex = BoxSize*ix+Overlap;
    if iy == 1 
        sy = 1;
        ey = BoxSize + Overlap;
    end

    if ix == 1
        sx = 1;
        ex = BoxSize + Overlap;
    end
    if iy == Ysub
        ey = ImSize(1);
    end

    if ix == Xsub
        ex = ImSize(2);
    end
    FinalIm(Zoom*sy-(Zoom-1):Zoom*ey,Zoom*sx-(Zoom-1):Zoom*ex,frame)=...
        FinalIm(Zoom*sy-(Zoom-1):Zoom*ey,Zoom*sx-(Zoom-1):Zoom*ex,frame)+ROIs(1:Ysize(nn),1:Xsize(nn),nn);
        
end
end         