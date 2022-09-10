function [OutIm,Cluster] = hist_MapN(Film, MinWeight, MaxStd, CutOff)
 %hist_MapN() uses the posterior image to find the emitters coordinates.
 %   Hierarchical clustering algorithm is used to find the emitters' 
 %   locations, the locations' uncertainties (sigma) and their
 %   weights. It then removes those found clusters with the weights smaller
 %   than the input cutoff (MinWeight) and also those with the uncertainties 
 %   larger than the input limit for uncertainty (MaxStd). It returns a 
 %   histogram image of the found centers and a matrix consisting of 5  
 %   columns are X and Ypositions of the centers, the uncertainties and the 
 %   weights.
 %
 % INPUTS:
 %   Film:        The input image of the data.
 %   MinWeight:   The cutoff of the wweights. The clusters with the weights
 %                smaller than the cutoff will be removed. 
 %   MaxStd:      The cutoff of the uncertainty. The clusters with
 %                uncertainties larger than cutoff will be removed. (pixels)
 %   Cutoff:      The furthest distance between two data points belonging
 %                to a cluster. It determins the size of the clister. (pixels)
 %
 % OUTPUTS:
 %   OutIm:       The histogram image of the center of the found clusters.
 %   Clust:       Structure containing the following fields:
 %     X:         The X-positions of the found clusters. (pixels)
 %     Y:         The Y-positions of the found clusters. (pixels)
 %     X_SE:      The uncertainty in X. (pixels)
 %     Y_SE:      The uncertainty in Y. (pixels)
 %     Weight:    The weight of each cluster. Note that it is not the
 %                number of photons for the probes.
 % 
 % REQUIRES:
 %   MATLAB 2014 or higher versions.
 %
 % CITATION:
 %   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
 %   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
 %   Fitting using Reversible Jump Markov Chain Monte Carlo".
 %
 % Created by:
 %   Mohamadreza Fazel (Lidke Lab 2017)
 %

OutIm = zeros(size(Film),'single'); 
[ii] = find(Film ~= 0); %finding the indices of the non-zero pixels 
Weight = Film(ii);
Data = zeros(length(ii),2); %saving the found indices in an array.
[i,j]=ind2sub(size(Film),ii);
Data(:,1) = j;
Data(:,2) = i;
% If there is no data points then return an empty structure.
if length(Data) < 3
  Cluster.X = [];
  Cluster.Y = [];
  Cluster.X_SE = [];
  Cluster.Y_SE = [];
  Cluster.Weight = [];
else
    % There are three functions for heirarchical clustering that we have to 
    % use them in turn. The first one is pdist which calculates the 
    % distances between every two pairs of the data points. This is the 
    % input to the next function which is linkage. The out put of this 
    % function the distance between each two data points which would be the
    % input of the next function which is linkage.
    
    Y = pdist(Data, 'euclid');
    
    TC = linkage(Y,'single');
    
    % cluster is the third function that we need to use. The input to this
    % function is the output of the previous one. This function has several
    % methods that we need to decide which one is suitable for our purpose. 
    % For example the method distance cluster the data based on the distance
    % between the data points.
    %fprintf('cluster, TC:%g, CutOff:%g\n',numel(TC),numel(CutOff));
    AC = cluster(TC,'Cutoff',CutOff,'Criterion','distance'); 
    %fprintf('clusterEnd\n');
    m = max(AC);
    MeanStd = zeros(m,5);

    for nn = 1:m
       a = find(AC == nn);
       MeanStd(nn,1) = sum(j(a).*Weight(a)) / sum(Weight(a));
       MeanStd(nn,3) = sum(i(a).*Weight(a)) / sum(Weight(a));
       MeanStd(nn,2) = sqrt(sum(Weight(a).*(j(a)-MeanStd(nn,1)).^2) / sum(Weight(a)));
       MeanStd(nn,4) = sqrt(sum(Weight(a).*(i(a)-MeanStd(nn,3)).^2) / sum(Weight(a)));
       MeanStd(nn,5) = sum(Weight(a));
       
    end
    % removing the clusters with small weights.
    MeanStd(MeanStd(:,5) < MinWeight,:) = []; 
    % removing clusters with large uncertainty.
    Std = (MeanStd(:,2).^2 + MeanStd(:,4).^2).^0.5;
    MeanStd(Std > MaxStd,:) = [];
    Cluster.X = MeanStd(:,1);
    Cluster.Y = MeanStd(:,3);
    Cluster.X_SE = MeanStd(:,1);
    Cluster.Y_SE = MeanStd(:,4);
    Cluster.Weight = MeanStd(:,5);
    % making the histogram image.
    for nn = 1:size(MeanStd,1)
        OutIm(round(MeanStd(nn,3)),round(MeanStd(nn,1)))=1;
    end
end  
end
