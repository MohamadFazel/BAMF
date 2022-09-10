function obj=makeClustInfo(obj)
%makeClustInfo() creates empty cluster called ClustInfo.
%
% INPUTS:
%   obj:   The class structure.
%
% OUTPUTS: 
%   obj:   The class structure containing empty structre, ClustInfo.
%
% REQUIRES:
%   MATLSB 2014a or higher versions.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%
obj.ClustInfo.X = single([]);
obj.ClustInfo.Y = single([]);
obj.ClustInfo.I = single([]);
obj.ClustInfo.X_SE = single([]);
obj.ClustInfo.Y_SE = single([]);
obj.ClustInfo.I_SE = single([]);
obj.ClustInfo.X_kmeans = single([]);
obj.ClustInfo.Y_kmeans = single([]);
obj.ClustInfo.I_kmeans = single([]);
obj.ClustInfo.X_SE_kmeans = single([]);
obj.ClustInfo.Y_SE_kmeans = single([]);
obj.ClustInfo.I_SE_kmeans = single([]);
obj.ClustInfo.Xbg = single([]);
obj.ClustInfo.Ybg = single([]);
obj.ClustInfo.Ibg = single([]);
obj.ClustInfo.Xbg_SE = single([]);
obj.ClustInfo.Ybg_SE = single([]);
obj.ClustInfo.Ibg_SE = single([]);
obj.ClustInfo.FrameNum = single([]);
obj.ClustInfo.FrameNumBg = single([]);
obj.ClustInfo.FrameNum_kmeans = single([]);
obj.ClustInfo.ROInum = single([]);
obj.ClustInfo.Bg = single([]);
obj.ClustInfo.ABg = single([]);
obj.ClustInfo.BBg = single([]);
obj.ClustInfo.JumpXRJ = single([]);
obj.ClustInfo.JumpIRJ = single([]);
obj.ClustInfo.JumpBgRJ = single([]);
obj.ClustInfo.JumpXFg = single([]);
obj.ClustInfo.JumpIFg = single([]);
obj.ClustInfo.JumpXBg = single([]);
obj.ClustInfo.JumpIBg = single([]);
obj.ClustInfo.SplitRJ = single([]);
obj.ClustInfo.MergeRJ = single([]);
obj.ClustInfo.BirthRJ = single([]);
obj.ClustInfo.DeathRJ = single([]);
obj.ClustInfo.JumpBgRJ = single([]);
obj.ClustInfo.JumpMC = single([]);
obj.ClustInfo.ROIcornerX = single([]);
obj.ClustInfo.ROIcornerY = single([]);
end