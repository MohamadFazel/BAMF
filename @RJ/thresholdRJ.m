function obj=thresholdRJ(obj)
%thresholdRJ() threshold the MAPN results from MCMC using the given linits.
%
% INPUTS:
%   obj:   The class object.
%
% OUTPUTS:
%   obj:   The class object.
%
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
Ind = obj.TotClust.X_SE < obj.Threshold.Max_XY_SE & obj.TotClust.Y_SE < obj.Threshold.Max_XY_SE & ...
          obj.TotClust.X_SE > obj.Threshold.Min_XY_SE & obj.TotClust.Y_SE > obj.Threshold.Min_XY_SE & ...
          obj.TotClust.Photons_SE < obj.Threshold.Max_Photons_SE & obj.TotClust.Photons > obj.Threshold.MinPhotons;
obj.TotClust.X = obj.TotClust.X(Ind);
obj.TotClust.Y = obj.TotClust.Y(Ind);
obj.TotClust.Photons = obj.TotClust.Photons(Ind);
obj.TotClust.X_SE = obj.TotClust.X_SE(Ind);
obj.TotClust.Y_SE = obj.TotClust.Y_SE(Ind);
obj.TotClust.Photons_SE = obj.TotClust.Photons_SE(Ind);
obj.TotClust.FrameNum = obj.TotClust.FrameNum(Ind);
obj.TotClust.ROInum = obj.TotClust.ROInum(Ind);
obj.TotClust.Bg = obj.TotClust.Bg;
obj.TotClust.ABg = obj.TotClust.ABg;
obj.TotClust.BBg = obj.TotClust.BBg;
obj.TotClust.JumpRJ = obj.TotClust.JumpRJ;
obj.TotClust.SplitRJ = obj.TotClust.SplitRJ;
obj.TotClust.MergeRJ = obj.TotClust.MergeRJ;
obj.TotClust.BirthRJ = obj.TotClust.BirthRJ;
obj.TotClust.DeathRJ = obj.TotClust.DeathRJ;
obj.TotClust.JumpBgRJ = obj.TotClust.JumpBgRJ;
obj.TotClust.JumpMC = obj.TotClust.JumpMC;
obj.TotClust.ROIcornerX = obj.TotClust.ROIcornerX(Ind);
obj.TotClust.ROIcornerY = obj.TotClust.ROIcornerY(Ind);
obj.TotClust.DatasetNum = obj.TotClust.DatasetNum(Ind);       
end