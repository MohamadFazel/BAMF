function [ROI_Data,ROI_Struct,ROI_Var]=makeSubRegions(Data,BoxSize,Overlap,Camera_Var)
%makeSubregions() split the input data frames into sub-regions. 
%   The sub-regions contain the overlapping regions but for those at the 
%   edges of the frame the overlapping region might be smaller so the ROIs  
%   do not have the same size. Zero pixels are added to the smaller ones so 
%   that they  can be saved in a 3D array. ROI_struct contains information 
%   about the  ROIs, like their sizes and their frame index and their index 
%   inside the frame.
% 
% INPUTS:
%   Data:      The raw super-resolution data.
%   BoxSize:   The size of ROIs.
%   Overlap:   The size of the overlapping region with the neigboring ROIs.
%   Camera_Var:The variance scaled by the gain for sCMOS camera, which would 
%              be split into subregions too. It is only required for sCMOS.
%
% OUTPUTS:
%   Data_ROI:  The ROIs.
%   ROI_Struct: Structure containing the following fields:
%   Xsize:     Size of the ROIs along the X-axis, or number of their columns.
%   Ysize:     Size of the ROIs along the Y-aixs, or number of their rows.
%   FrameNum:  The number of the frame that a ROI was cropped from that.
%   ROI_Ind:   The index of the ROI inside the frame.
%   Xsub:      The number of the ROIs along the X-axis.
%   Ysub:      The number of the ROIs along the Y-axis.
%
%   ROI_Var:   The subregions of the sCMOS camera scaled variance that
%              correspond to the ROIs.
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

[YPix,XPix,ZPix] = size(Data); %finding the size of the data
Xsub = ceil(XPix/BoxSize); %number of the sub-regions along the X-axis
Ysub = ceil(YPix/BoxSize); %number of the sub-regions along the Y-axis
ROI_Data = zeros(BoxSize+2*Overlap,BoxSize+2*Overlap,Xsub*Ysub*ZPix,'single');
if nargin < 4 
   ROI_Var = zeros(size(ROI_Data),'single');
elseif isscalar(Camera_Var)
   ROI_Var = zeros(size(ROI_Data),'single'); 
elseif isempty(Camera_Var)
   ROI_Var = zeros(size(ROI_Data),'single'); 
end
ROI_Struct.Xsize = zeros(1,Xsub*Ysub*ZPix);
ROI_Struct.Ysize = zeros(1,Xsub*Ysub*ZPix);
ROI_Struct.FrameNum = zeros(1,Xsub*Ysub*ZPix);
ROI_Struct.ROI_Ind = zeros(1,Xsub*Ysub*ZPix);
%The following loop goes through different frames
for frame = 1:ZPix
    %The two following loops crop the suregions in each frame
    for ix = 1:Xsub
        for iy = 1:Ysub
            ROIex = BoxSize + 2*Overlap;
            ROIey = BoxSize + 2*Overlap;
            sy = BoxSize*(iy-1)-(Overlap-1);
            ey = BoxSize*iy+Overlap;
            sx = BoxSize*(ix-1)-(Overlap-1);
            ex = BoxSize*ix+Overlap;
            if iy == 1 
                ROIey = BoxSize + Overlap;
                sy = 1;
                ey = BoxSize + Overlap;
            end

            if ix == 1
                ROIex = BoxSize + Overlap;
                sx = 1;
                ex = BoxSize + Overlap;
            end
            if iy == Ysub
                ey = YPix;
                ROIey = YPix-BoxSize*(Ysub-1)+Overlap;
            end

            if ix == Xsub
                ex = XPix;
                ROIex = XPix-BoxSize*(Xsub-1)+Overlap;
            end
            if Xsub == 1
                ROIex = XPix;
            end
            if Ysub == 1
                ROIey = YPix;
            end
            Ind = Xsub*Ysub*(frame-1) + Ysub*(ix-1) + iy;
            ROI_Struct.Xsize(Ind) = ROIex; 
            ROI_Struct.Ysize(Ind) = ROIey;
            ROI_Struct.FrameNum(Ind) = frame;
            ROI_Struct.ROI_Ind(Ind) = Ysub*(ix-1) + iy;
            ROI_Data(1:ROIey,1:ROIex,Ind) = single(Data(sy:ey,sx:ex,frame));
        end
    end
end
ROI_Struct.SizeData = [YPix,XPix,ZPix];
ROI_Struct.Xsub = Xsub;
ROI_Struct.Ysub = Ysub;
ROI_Struct.Overlap = Overlap;
%fprintf('makeSubRegionsEnd\n');
end