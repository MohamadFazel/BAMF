function processData(obj,TimeID)
%processData() loads the datasets and called the calss to process them.
%
% INPUTS:
%   obj:    The class object.
%   TimeID: The data set number.
%
% OUTPUTS:
%
% REQUIRES:
%   MATLAB 2014a or higher versions.
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
DataName = fullfile(obj.DataDir,obj.FileList{TimeID});
obj.Data = RJ.loadData(DataName,'sequence');
if isempty(obj.FileListDriftX)
    obj.RJStruct.XDrift = zeros(size(obj.Data,3),1);
else
    Drift = load(fullfile(obj.DriftXDir,obj.FileListDriftX{TimeID}));
    Names = fieldname(Drift);
    obj.RJStruct.XDrift = -getfield(Drift,Names{1});
end
if isempty(obj.FileListDriftY)
    obj.RJStruct.YDrift = zeros(size(obj.Data,3),1);
else
    Drift = load(fullfile(obj.DriftYDir,obj.FileListDriftY{TimeID}));
    Names = fieldname(Drift);
    obj.RJStruct.YDrift = -getfield(Drift,Names{1});
end
obj.analyzeData();
SP = cat(2,'RJ_Result_',obj.FileList{TimeID});
ClustInfo = obj.ClustInfo;
PImg = obj.PImage;
save(fullfile(obj.SaveDir,SP),'ClustInfo','PImg')
S = sprintf('Priors_DataSet#%03g.png',TimeID);
print(gcf,'-dpng',fullfile(obj.SaveDir,S)) 
end
