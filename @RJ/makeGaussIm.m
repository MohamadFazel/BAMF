function SRImage=makeGaussIm(obj,ClustInfoT)
%makeGaussIm() takes the coordinates and generates the image with blobs.
%   This function takes the coordinates and calls some other functions 
%   (which may not be available) to creates super-resolved image containing
%   Gaussian blobs with the size of the returned uncertainties.
% 
% INPUTS:
%   obj:       The class object.
%   ClustInfo: The structure containing the coordinates and uncertainties.
%
% OUTPUTS:
%   SRImage:   The reconstructed image using the found localizations.
%
% REQUIRES:
%   NVIDIA GPU
%   SMA_Vis
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
SM.X = ClustInfoT.X;
SM.Y = ClustInfoT.Y;
SM.Photons = ClustInfoT.I;
SM.X_SE = ClustInfoT.X_SE;
SM.Y_SE = ClustInfoT.Y_SE;
SM.Bg = zeros(size(SM.X));
SM.FrameNum = ClustInfoT.FrameNum;
SM.XSize=size(obj.Data,2);
SM.YSize=size(obj.Data,1);
Zoom = obj.RJStruct.SRZoom;
SRImage=SMA_Vis.gaussianImage(SM,Zoom);
end

