function [ClustInfo,PImage]=makeResultsRJ(obj)
%makeResultsRJ() gather the final results from different datasets.
%
% INPUTS:
%   obj:       The class object.
% 
% OUTPUTS:
%   ClustInfo: Structure containing the emitters information.
%   PImage:    The final posterior image.
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
ClustInfo.X = single([]);
ClustInfo.Y = single([]);
ClustInfo.I = single([]);
ClustInfo.X_SE = single([]);
ClustInfo.Y_SE = single([]);
ClustInfo.I_SE = single([]);
ClustInfo.X_kmeans = single([]);
ClustInfo.Y_kmeans = single([]);
ClustInfo.I_kmeans = single([]);
ClustInfo.X_SE_kmeans = single([]);
ClustInfo.Y_SE_kmeans = single([]);
ClustInfo.I_SE_kmeans = single([]);
ClustInfo.Xbg = single([]);
ClustInfo.Ybg = single([]);
ClustInfo.Ibg = single([]);
ClustInfo.Xbg_SE = single([]);
ClustInfo.Ybg_SE = single([]);
ClustInfo.Ibg_SE = single([]);
ClustInfo.Framenum = single([]);
ClustInfo.ROInum = single([]);
ClustInfo.Bg = single([]);
ClustInfo.ABg = single([]);
ClustInfo.BBg = single([]);
ClustInfo.JumpRJ = single([]);
ClustInfo.SplitRJ = single([]);
ClustInfo.MergeRJ = single([]);
ClustInfo.BirthRJ = single([]);
ClustInfo.DeathRJ = single([]);
ClustInfo.JumpBgRJ = single([]);
ClustInfo.JumpMC = single([]);
ClustInfo.ROIcornerX = single([]);
ClustInfo.ROIcornerY = single([]);
ClustInfo.DatasetNum = single([]);
PImage = 0;

FileList = dir(fullfile(obj.DataDir,'*.mat'));
N = size(FileList,1);

for ii = 1:N
    PImage = PImage + obj.PImg;
    ClustInfo.X = cat(1,ClustInfo.X,obj.ClustInfo.X');
    ClustInfo.Y = cat(1,ClustInfo.Y,obj.ClustInfo.Y');
    ClustInfo.I = cat(1,ClustInfo.I,obj.ClustInfo.I');
    ClustInfo.X_SE = cat(1,ClustInfo.X_SE,obj.ClustInfo.X_SE');
    ClustInfo.Y_SE = cat(1,ClustInfo.Y_SE,obj.ClustInfo.Y_SE');
    ClustInfo.I_SE = cat(1,ClustInfo.I_SE,obj.ClustInfo.I_SE');
    ClustInfo.X_kmeans = cat(1,ClustInfo.X_kmeans,obj.ClustInfo.X_kmeans');
    ClustInfo.Y_kmeans = cat(1,ClustInfo.Y_kmeans,obj.ClustInfo.Y_kmeans');
    ClustInfo.I_kmeans = cat(1,ClustInfo.I_kmeans,obj.ClustInfo.I_kmeans');
    ClustInfo.X_SE_kmeans = cat(1,ClustInfo.X_SE_kmeans,obj.ClustInfo.X_SE_kmeans');
    ClustInfo.Y_SE_kmeans = cat(1,ClustInfo.Y_SE_kmeans,obj.ClustInfo.Y_SE_kmeans');
    ClustInfo.I_SE_kmeans = cat(1,ClustInfo.I_SE_kmeans,obj.ClustInfo.I_SE_kmeans');
    ClustInfo.Xbg = cat(1,ClustInfo.Xbg,obj.ClustInfo.Xbg');
    ClustInfo.Ybg = cat(1,ClustInfo.Ybg,obj.ClustInfo.Ybg');
    ClustInfo.Ibg = cat(1,ClustInfo.Ibg,obj.ClustInfo.Ibg');
    ClustInfo.Xbg_SE = cat(1,ClustInfo.Xbg_SE,obj.ClustInfo.Xbg_SE');
    ClustInfo.Ybg_SE = cat(1,ClustInfo.Ybg_SE,obj.ClustInfo.Ybg_SE');
    ClustInfo.Ibg_SE = cat(1,ClustInfo.Ibg_SE,obj.ClustInfo.Ibg_SE');
    ClustInfo.Bg = cat(1,ClustInfo.Bg,obj.ClustInfo.Bg');
    ClustInfo.ABg = cat(1,ClustInfo.ABg,obj.ClustInfo.ABg');
    ClustInfo.BBg = cat(1,ClustInfo.BBg,obj.ClustInfo.BBg');
    ClustInfo.Framenum = cat(1,ClustInfo.Framenum,obj.ClustInfo.Framenum');
    ClustInfo.ROInum = cat(1,ClustInfo.ROInum,obj.ClustInfo.ROInum');
    ClustInfo.JumpRJ = cat(1,ClustInfo.JumpRJ,obj.ClustInfo.JumpRJ');
    ClustInfo.SplitRJ = cat(1,ClustInfo.SplitRJ,obj.ClustInfo.SplitRJ');
    ClustInfo.MergeRJ = cat(1,ClustInfo.MergeRJ,obj.ClustInfo.MergeRJ');
    ClustInfo.BirthRJ = cat(1,ClustInfo.BirthRJ,obj.ClustInfo.BirthRJ');
    ClustInfo.DeathRJ = cat(1,ClustInfo.DeathRJ,obj.ClustInfo.DeathRJ');
    ClustInfo.JumpBgRJ = cat(1,ClustInfo.JumpBgRJ,obj.ClustInfo.JumpBgRJ');
    ClustInfo.JumpMC = cat(1,ClustInfo.JumpMC,obj.ClustInfo.JumpMC');
    ClustInfo.ROIcornerX = cat(1,ClustInfo.ROIcornerX,obj.ClustInfo.ROIcornerX');
    ClustInfo.ROIcornerY = cat(1,ClustInfo.ROIcornerY,obj.ClustInfo.ROIcornerY');
    ClustInfo.DatasetNum = cat(1,ClustInfo.DatasetNum,ii*ones(size(obj.ClustInfo.X')));
end

end