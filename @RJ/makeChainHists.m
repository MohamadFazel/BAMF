function [X,Y]=makeChainHists(obj,Variable,Index,PlotFlag)
% makeHists() gets the chain and according to the input Variable extract
% the chain of a specific variable and plot its histogram. The last input
% (Index) is optional. If the user provides the Index then the function
% plots the histogram of the chain corresponding to that index otherwise it
% plots the histogram of the entire chain. This method is designed for
% the diagnosis poprpuses
%
% INPUTS:
%   Chain:     The stored chain, which is a cell array.
%   Variable:  This is a string that shows what diagram the user wants to
%              plot. The options are: 'Photons' (intensity histogram), 'N' 
%              (Number of the emitters), 'LLR' (Likelihood ratio), 'PR'
%              (Prior ratio), 'JumpType' (The type of the jumps), 'BG'
%              (offset background), 'ABG' (slope of the  background along 
%              the X-axis), 'BBG' (slope of the background along the
%              Y-axis).
%   Index:     This input is optional. If this is provided then the plots
%              are the results of the chain corresponding to the ROI with
%              this index.
%
% OUTPUTS:     
%   X:         The array which stores the values corresponding to the input
%              string.
%   Y:         For 'Photons' and 'N' there are two categories of the
%              emitters, which are signal and background. For these two
%              specific case X and Y respectively corespond to signal and 
%              background.
% REQUIRES:
%   MATLAB 2016 or higher.
% 
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2017)
%
Chain = obj.Chain;
RJstruct = obj.RJStruct;
if nargin<3
   RJstruct = []; 
end
if nargin > 3 && ~isempty(Index)
    IterEnd = Index;
    IterBeg = Index;
else
    IterEnd = length(Chain);
    IterBeg = 1;
end
if nargin < 5
   PlotFlag = 1; 
end

X = [];
if strcmp(Variable,'Photons') || strcmp(Variable,'N')
    Y = [];
end
for nn = IterBeg:IterEnd
    SubChain = Chain{nn};
    Iter = length(SubChain);
    for ii = round(Iter/4):Iter
        SubSubChain = SubChain(ii);
       if SubSubChain.N ~= 0
        if strcmp(Variable,'Photons') 
            Signal = getfield(SubSubChain,'Signal');
            x = zeros(1,sum(Signal));
            y = zeros(1,length(Signal)-sum(Signal));
            xx = getfield(SubSubChain,Variable);
            a = 0;
            b = 0;
            for mm = 1:length(Signal)
                if Signal(mm)==1
                    a=a+1;
                    x(a)=xx(mm);
                else
                    b=b+1;
                    y(b)=xx(mm);
                end
            end
            X = cat(2,X,x);
            Y = cat(2,Y,y);
        elseif strcmp(Variable,'N')
            Signal = getfield(SubSubChain,'Signal');
             x = sum(Signal);%getfield(SubSubChain,Variable);
             X = cat(2,X,x);
             y = -sum(Signal-1);%getfield(SubSubChain,Variable);
             Y = cat(2,Y,y);
        else  
            x = getfield(SubSubChain,Variable);
            X = cat(2,X,x);
        end
      end
    end
end
if PlotFlag
switch Variable
    case 'Photons'
        S = sprintf('Signal Photons');
    case 'N'
        S = sprintf('Number of Emitters in the Signal');
    case 'LLR'
        S = sprintf('Likelihood Ratio');
    case 'PR'
        S = sprintf('Prior Ratio');
    case 'JumpType'
        S = sprintf('Jump Types');
    case 'BG'
        S = sprintf('Offset Background');
    case 'ABG'
        S = sprintf('Background X-slope');
    case 'BBG'
        S = sprintf('Background Y-slope');
end

figure; HistHandle=histogram(X,'NumBins',25,'Normalization','pdf');
xlabel(S,'FontSize',15)
ylabel('Frequency','FontSize',15)
if ~isempty(RJstruct)
   Xarray = (0:RJstruct.P_dP:HistHandle.BinEdges(end)+50*RJstruct.P_dP);
end
if strcmp(Variable,'JumpType')
    hold off
    hist(X);
    xlabel(S,'FontSize',15)
    ylabel('Frequency','FontSize',15)
    ax = gca;
    ax.XTickLabel={'Jump','','Split','','Merge','','Birth','','Death'};
    %title('Number of the proposed jumps, which were not rejected initially.');
end
if strcmp(Variable,'Photons')  
   if ~isempty(RJstruct)
       hold(gca,'on'); plot(Xarray,RJstruct.P_EmPhotons(1:length(Xarray)),'lineWidth',1.5) 
       legend('Found dist.','Used dist.')
   end
   figure; histogram(Y,'NumBins',25,'Normalization','pdf');
   ylabel('Frequency','FontSize',15)
   xlabel('Background Photons','FontSize',15)
   if ~isempty(RJstruct)
       hold(gca,'on'); plot(Xarray,RJstruct.P_BgEmPhotons(1:length(Xarray)),'lineWidth',1.5) 
       legend('Found dist.','Used dist.')
   end
end
if strcmp(Variable,'N')
    if ~isempty(RJstruct)
        Xarray = 0:HistHandle.BinEdges(end)+10;
        hold(gca,'on'); plot(Xarray,poisspdf(Xarray,RJstruct.Rho*16*16),'lineWidth',1.5); 
        legend('Found dist.','Used dist.');
        figure;histogram(Y,'NumBins',25,'Normalization','pdf');
        ylabel('Frequency','FontSize',15)
        xlabel('Number of the Emitters in the Background','FontSize',15)
    else
        hold off
        hist(X,25);
        ylabel('Frequency','FontSize',15)
        xlabel(S,'FontSize',15)
        figure; hist(Y,25); 
        ylabel('Frequency','FontSize',15)
        xlabel('Number of the Emitters in the Background','FontSize',15)
    end
    
end
if strcmp(Variable,'BG')&&~isempty(RJstruct)
    Xarray = (0:RJstruct.P_dP:HistHandle.BinEdges(end)+5*RJstruct.P_dP);
    hold(gca,'on'); plot(Xarray,RJstruct.P_Offset(1:length(Xarray)),'lineWidth',1.5) 
    legend('Found dist.','Used dist.')
end
if strcmp(Variable,'ABG')||strcmp(Variable,'BBG')&&~isempty(RJstruct)
    Xarray = (HistHandle.BinEdges(1)-20:HistHandle.BinEdges(end)+20);
    hold(gca,'on'); plot(Xarray,normpdf(Xarray,0,sqrt(2)),'lineWidth',1.5) 
    legend('Found dist.','Used dist.')
end

end
end