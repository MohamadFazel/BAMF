function obj=assembleResults(obj)
%assembleResults() puts together results from different datasets.
%   The code process each dataset and save the result from
%   them independently. assempleResults() loads these results and put them
%   in a single structure.
%
% INPUTS:
%   obj:   obj is the class object containing the address of the saved
%          results.
%
% OUTPUTS: 
%   obj:   obj is the class structure where the structure containing the
%          information from all the datasets is saved.
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
%   Mohamadreza Fazel (Lidek Lab 2017)
%

obj.TotClust.X = single([]);
obj.TotClust.Y = single([]);
obj.TotClust.Photons = single([]);
obj.TotClust.X_SE = single([]);
obj.TotClust.Y_SE = single([]);
obj.TotClust.Photons_SE = single([]);
obj.TotClust.X_kmeans = single([]);
obj.TotClust.Y_kmeans = single([]);
obj.TotClust.Photons_kmeans = single([]);
obj.TotClust.X_SE_kmeans = single([]);
obj.TotClust.Y_SE_kmeans = single([]);
obj.TotClust.Photons_SE_kmeans = single([]);
obj.TotClust.Xbg = single([]);
obj.TotClust.Ybg = single([]);
obj.TotClust.Photonsbg = single([]);
obj.TotClust.Xbg_SE = single([]);
obj.TotClust.Ybg_SE = single([]);
obj.TotClust.Photonsbg_SE = single([]);
obj.TotClust.FrameNum = single([]);
obj.TotClust.ROInum = single([]);
obj.TotClust.Bg = single([]);
obj.TotClust.ABg = single([]);
obj.TotClust.BBg = single([]);
obj.TotClust.JumpRJ = single([]);
obj.TotClust.SplitRJ = single([]);
obj.TotClust.MergeRJ = single([]);
obj.TotClust.BirthRJ = single([]);
obj.TotClust.DeathRJ = single([]);
obj.TotClust.JumpBgRJ = single([]);
obj.TotClust.JumpMC = single([]);
obj.TotClust.ROIcornerX = single([]);
obj.TotClust.ROIcornerY = single([]);
obj.TotClust.DatasetNum = single([]);
PImage = 0;

FileList = dir(fullfile(obj.SaveDir,'*.mat'));
N = size(FileList,1);

for ii = 1:N
    SP = FileList(ii).name;
    H = load(fullfile(obj.SaveDir,SP));
    PImage = PImage + H.PImg;
    obj.TotClust.X = cat(1,obj.TotClust.X,H.ClustInfo.X');
    obj.TotClust.Y = cat(1,obj.TotClust.Y,H.ClustInfo.Y');
    obj.TotClust.Photons = cat(1,obj.TotClust.Photons,H.ClustInfo.I');
    obj.TotClust.X_SE = cat(1,obj.TotClust.X_SE,H.ClustInfo.X_SE');
    obj.TotClust.Y_SE = cat(1,obj.TotClust.Y_SE,H.ClustInfo.Y_SE');
    obj.TotClust.Photons_SE = cat(1,obj.TotClust.Photons_SE,H.ClustInfo.I_SE');
    obj.TotClust.X_kmeans = cat(1,obj.TotClust.X_kmeans,H.ClustInfo.X_kmeans');
    obj.TotClust.Y_kmeans = cat(1,obj.TotClust.Y_kmeans,H.ClustInfo.Y_kmeans');
    obj.TotClust.Photons_kmeans = cat(1,obj.TotClust.Photons_kmeans,H.ClustInfo.I_kmeans');
    obj.TotClust.X_SE_kmeans = cat(1,obj.TotClust.X_SE_kmeans,H.ClustInfo.X_SE_kmeans');
    obj.TotClust.Y_SE_kmeans = cat(1,obj.TotClust.Y_SE_kmeans,H.ClustInfo.Y_SE_kmeans');
    obj.TotClust.Photons_SE_kmeans = cat(1,obj.TotClust.Photons_SE_kmeans,H.ClustInfo.I_SE_kmeans');
    obj.TotClust.Xbg = cat(1,obj.TotClust.Xbg,H.ClustInfo.Xbg');
    obj.TotClust.Ybg = cat(1,obj.TotClust.Ybg,H.ClustInfo.Ybg');
    obj.TotClust.Photonsbg = cat(1,obj.TotClust.Photonsbg,H.ClustInfo.Ibg');
    obj.TotClust.Xbg_SE = cat(1,obj.TotClust.Xbg_SE,H.ClustInfo.Xbg_SE');
    obj.TotClust.Ybg_SE = cat(1,obj.TotClust.Ybg_SE,H.ClustInfo.Ybg_SE');
    obj.TotClust.Photonsbg_SE = cat(1,obj.TotClust.Photonsbg_SE,H.ClustInfo.Ibg_SE');
    obj.TotClust.Bg = cat(1,obj.TotClust.Bg,H.ClustInfo.Bg');
    obj.TotClust.ABg = cat(1,obj.TotClust.ABg,H.ClustInfo.ABg');
    obj.TotClust.BBg = cat(1,obj.TotClust.BBg,H.ClustInfo.BBg');
    obj.TotClust.FrameNum = cat(1,obj.TotClust.FrameNum,H.ClustInfo.FrameNum');
    obj.TotClust.ROInum = cat(1,obj.TotClust.ROInum,H.ClustInfo.ROInum');
    obj.TotClust.JumpRJ = cat(1,obj.TotClust.JumpRJ,H.ClustInfo.JumpXRJ');
    obj.TotClust.SplitRJ = cat(1,obj.TotClust.SplitRJ,H.ClustInfo.SplitRJ');
    obj.TotClust.MergeRJ = cat(1,obj.TotClust.MergeRJ,H.ClustInfo.MergeRJ');
    obj.TotClust.BirthRJ = cat(1,obj.TotClust.BirthRJ,H.ClustInfo.BirthRJ');
    obj.TotClust.DeathRJ = cat(1,obj.TotClust.DeathRJ,H.ClustInfo.DeathRJ');
    obj.TotClust.JumpBgRJ = cat(1,obj.TotClust.JumpBgRJ,H.ClustInfo.JumpBgRJ');
    obj.TotClust.JumpMC = cat(1,obj.TotClust.JumpMC,H.ClustInfo.JumpMC');
    obj.TotClust.ROIcornerX = cat(1,obj.TotClust.ROIcornerX,H.ClustInfo.ROIcornerX');
    obj.TotClust.ROIcornerY = cat(1,obj.TotClust.ROIcornerY,H.ClustInfo.ROIcornerY');
    obj.TotClust.DatasetNum = cat(1,obj.TotClust.DatasetNum,ii*ones(size(H.ClustInfo.X')));
end
RJ.produceIm(PImage,fullfile(obj.SaveDir,'Posterior'));
end
