function ModGrid=removeOverlapping(ROI_Struct,Overlap,Zoom,ROI_Stack)
% removeOverlapping() removes the chain inside the overlapping regions.
%
% INPUTS:
%   ROI_Struct:     A structure containing info on the ROIs like theri
%                   sizes etc.
%   Overlap:        The width of the overlapping region in pixels.
%   Zoom:           The magnification of the image of the chain.
%   ROI_Stack:      The stack of the images of the chain corresponding to
%                   the processed ROIs.
%   
% OUTPUTS:
%   ModGrid:        The stack of the input images where the overlapping 
%                   parts are removed.
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

ROInum = size(ROI_Stack,3);
ModGrid = ROI_Stack;
ROI_Ind = ROI_Struct.ROI_Ind; 
Xsub = ROI_Struct.Xsub;
Ysub = ROI_Struct.Ysub;
Xsize = ROI_Struct.Xsize*Zoom;
Ysize = ROI_Struct.Ysize*Zoom;
for nn = 1:ROInum
    Ind = ROI_Ind(nn);
    ix = floor((Ind-1)/Ysub) + 1;
    iy = Ind - (ix-1)*Ysub;
    if ix ~= Xsub %ceil(length(obj.Data')/obj.SubIm) %if this is not the last block in the row
       ModGrid(:,Xsize(nn)-(Zoom*Overlap-1):Xsize(nn),nn) = 0;
    end

    if iy ~= Ysub 
       ModGrid(Ysize(nn)-(Zoom*Overlap-1):Ysize(nn),:,nn)=0;
    end

    if ix ~= 1 
       ModGrid(:,1:Zoom*Overlap,nn)=0;
    end

    if  iy ~= 1 
        ModGrid(1:Zoom*Overlap,:,nn)=0;
    end
end

end