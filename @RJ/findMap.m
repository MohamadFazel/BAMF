function [SMD,MostFrequent]=findMap(Chain,FB)
%findMap() finds and extracts the most probable model from the chain.
%   findMap() finds the moset repeated model and then extracts all the 
%   states of that model from the chain. It then uses kmeans algorithm to
%   find the locations and intensities of the emitters.
%
% INPUTS:
%   Chain:   The chain of the accepted jumps.
%   FB:      It can be either 0 or 1. 1 returns the signal information of
%            the most probable model. 0 returns the structured background
%            emitters.
%
% OUTPUTS:
%   SMD:     Structure containing the emitters parameters.
%   MostFrequent: The most frequent number of the emitters.
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

Iter = length(Chain);
NN = zeros(1,Iter-round(Iter/2));
NInd = 0;
for nn=round(Iter/2):Iter
    NInd = NInd + 1;
    if Chain(nn).Signal ~= -1
       NN(NInd)=sum(Chain(nn).Signal==FB);
    end
end
MaxN=max(NN);
Frequency=hist(NN,(0:MaxN));
MostFrequent = find(Frequency==max(Frequency))-1;
if length(MostFrequent)>1
   MostFrequent = MostFrequent(1); 
elseif isempty(MostFrequent)
   MostFrequent = 0; 
end
X=[];
Y=[];
I=[];
Bg=[];
ABg=[];
BBg=[];
SMD.X=[];
SMD.X_SE=[];
SMD.Y=[];
SMD.Y_SE=[];
SMD.I=[];
SMD.I_SE=[];
SMD.Bg=[];
SMD.Bg_SE=[];
SMD.ABg=[];
SMD.BBg=[];
C=0;
for nn=round(Iter/2):Iter
    Strc = Chain(nn);
    N = sum(Strc.Signal==FB);
    if N==MostFrequent
       C=C+1;
       Xx=Strc.X(Strc.Signal==FB);
       Yy=Strc.Y(Strc.Signal==FB);
       Ii=Strc.Photons(Strc.Signal==FB);
       BGg=Strc.BG;
       ABGg=Strc.ABG;
       BBGg=Strc.BBG;
       X=cat(2,X,Xx);
       Y=cat(2,Y,Yy);
       I=cat(2,I,Ii);
       Bg=cat(2,Bg,BGg);
       ABg=cat(2,ABg,ABGg);
       BBg=cat(2,BBg,BBGg);
    end
end
ID = [];
if MostFrequent > 0 
    Points(:,1)=X';
    Points(:,2)=Y';
    ID = kmeans(Points,MostFrequent);
    XMean = zeros(1,MostFrequent);
    YMean = zeros(1,MostFrequent);
    IMean = zeros(1,MostFrequent); 
    X_SE = zeros(1,MostFrequent);
    Y_SE = zeros(1,MostFrequent);
    I_SE = zeros(1,MostFrequent);
    for nn = 1:MostFrequent
        Ind=find(ID == nn);
        XMean(nn)=mean(X(Ind));
        X_SE(nn)=std(X(Ind));
        YMean(nn)=mean(Y(Ind));
        Y_SE(nn)=std(Y(Ind));
        IMean(nn)=mean(I(Ind));
        I_SE(nn)=std(I(Ind));
    end
    SMD.X = XMean;
    SMD.Y = YMean;
    SMD.I = IMean;
    SMD.X_SE = X_SE;
    SMD.Y_SE = Y_SE;
    SMD.I_SE = I_SE;
end
SMD.Bg = mean(Bg);
SMD.Bg_SE = std(Bg);
SMD.ABg = mean(ABg);
SMD.BBg = mean(BBg);
end