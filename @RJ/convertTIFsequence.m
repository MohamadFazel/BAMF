function [Out]=convertTIFsequence(Path,Name)
%converTIFsequence() converts the TIF images to mat-files.
%
% INPUTS:
%   Path:   The path to the TIF-file.
%   Name:   Name of the TIF-file.
%
% OUTPUTS:
%   Out:    The sequence of the data.
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
Out=readtimeseries(Path);
sequence=uint32(Out);
save(Name,'sequence','-v7.3')
end