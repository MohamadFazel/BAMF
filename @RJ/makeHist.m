function makeHist(obj,Struct,Variable,PlotFlag,SaveFlag)
%makeHist() returns the histogram of the input parameters.
%   The input parameter is stracted from one of the structures that are 
%   made at different stages of the algorithm based on the input and the 
%   function plots their histogram.
%
% INPUTS: 
%   obj:      The class object.
%   Struct:   The structure that the user wants to extract the parameter
%             from. It could be 'ClustInfo', 'TotClust', 'SMD' and 
%             'SMD_combined'
%   Variablr: The variable that the user is interested in. This can be 
%             'All', 'Photons', 'X_SE', 'Y_SE', 'Photons_SE', 'BG', 'ABG',
%             'BBG', 'JumpRJ', 'SplitRJ', 'MergeRJ', 'BirthRJ', 'DeathRJ',
%             'JumpMC'.
%   PlotFlag: 1 displays the plots and 0 does not.
%   SaveFlag: 1 saves the plots in the save directory and 0 does not.
%
% OUTPUTS:
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
if strcmp('ClustInfo',Struct)
    Clust = obj.ClustInfo;
elseif strcmp('TotClust',Struct)
    Clust = obj.TotClust;
elseif strcmp('SMD',Struct)
    Clust = obj.SMD;
elseif strcmp('SMD_combined',Struct)
    Clust = obj.SMD_combined;
else
    error('This structure does not exist.');
end
if strcmp('All',Variable)
   All = 1; 
else
   All = 0; 
end
if PlotFlag
   Visible = 'on'; 
else
   Visible = 'off'; 
end
if strcmp('Photons',Variable) || All
    figure('Visible',Visible)
    if strcmp(Struct,'ClustInfo')
        hist(Clust.I,25)
    else
        hist(Clust.Photons,25)
    end
    xlabel('Photons');ylabel('Frequency');title('Photons')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'Photons.png')) 
    end
end
if strcmp('X_SE',Variable) || All
    figure('Visible',Visible)
    hist(Clust.X_SE,25)
    xlabel('X-SE');ylabel('Frequency');title('X-SE')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'X-SE.png')) 
    end
end
if strcmp('Y_SE',Variable) || All
    figure('Visible',Visible)
    hist(Clust.Y_SE,25)
    xlabel('Y-SE');ylabel('Frequency');title('Y-SE')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'Y-SE.png')) 
    end
end
if strcmp('Photons_SE',Variable) || All
    figure('Visible',Visible)
    if strcmp(Struct,'ClustInfo')
        hist(Clust.I_SE,25)
    else
        hist(Clust.Photons_SE,25)
    end
    xlabel('Photons-SE');ylabel('Frequency');title('Photons-SE')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'Photons-SE.png')) 
    end
end
if strcmp('BG',Variable) || All
    figure('Visible',Visible)
    hist(Clust.Bg,25)
    xlabel('BG');ylabel('Frequency');title('BG')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'BG.png')) 
    end
end
if strcmp('ABG',Variable) || All
    figure('Visible',Visible)
    hist(Clust.ABg,25)
    xlabel('BG slope along X-axis');ylabel('Frequency');title('ABG')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'ABG.png')) 
    end
end
if strcmp('BBG',Variable) || All
    figure('Visible',Visible)
    hist(Clust.BBg,25)
    xlabel('BG slope along Y-axis');ylabel('Frequency');title('BBG')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'BBG.png')) 
    end
end
if strcmp('JumpRJ',Variable) || All
    figure('Visible',Visible)
    if strcmp(Struct,'ClustInfo')
        hist(Clust.JumpXRJ,25)
    else
        hist(Clust.JumpRJ,25)
    end
    xlabel('ROIs RJ-jump acceptance ratio');ylabel('Frequency');title('JumpRJ')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'JumpRJ.png')) 
    end
end
if strcmp('SplitRJ',Variable) || All
    figure('Visible',Visible)
    hist(Clust.SplitRJ,25)
    xlabel('ROIs split acceptance ratio');ylabel('Frequency');title('SplitRJ')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'SplitRJ.png')) 
    end
end
if strcmp('MergeRJ',Variable) || All
    figure('Visible',Visible)
    hist(Clust.MergeRJ,25)
    xlabel('ROIs merge acceptance ratio');ylabel('Frequency');title('MergeRJ')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'MergeRJ.png')) 
    end
end
if strcmp('BirthRJ',Variable) || All
    figure('Visible',Visible)
    hist(Clust.BirthRJ,25)
    xlabel('ROIs birth acceptance ratio');ylabel('Frequency');title('BirthRJ')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'BirthRJ.png')) 
    end
end
if strcmp('DeathRJ',Variable) || All
    figure('Visible',Visible)
    hist(Clust.DeathRJ,25)
    xlabel('ROIs death acceptance ratio');ylabel('Frequency');title('DeathRJ')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'DeathRJ.png')) 
    end
end
if strcmp('JumpMC',Variable) || All
    figure('Visible',Visible)
    hist(Clust.JumpMC,25)
    xlabel('MCMC jumps acceptance ratio');ylabel('Frequency');title('JumpMC')
    if SaveFlag
        print(gcf,'-dpng',fullfile(obj.SaveDir,'JumpMC.png')) 
    end
end
end