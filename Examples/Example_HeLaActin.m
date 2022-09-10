
%RJ_demo is a demonstartion of RJ-class for a small section of
%experimental data in Fig. 3.

SaveDir = 'Actin';
addpath(fullfile('..','ptx'))
addpath(fullfile('..','mex'))
addpath('..')
%Loading data and drifts.
load(fullfile(SaveDir,'ActinData.mat'))
load(fullfile(SaveDir,'Drift'))
%Initializing the class.
B=RJ('NoGui');
%The PSF size for this data-set is 1.15 pixel.
B.FindPSF = 1; 
%Gain and offset.
B.Camera_Offset=90;
B.Camera_Gain=13;
%setting the jumps size in X to 0.05 pixel.
B.RJStruct.XstdFg = 0.05;
B.RJStruct.XDrift = XDrift;
B.RJStruct.YDrift = YDrift;
B.Data=sequence;
%running the code
B.analyzeData()

%Saving the posterior image and the histogram of some of the parameters.
RJ.produceIm(B.PImage,fullfile(SaveDir,'Posterior'));
PlotFlag = 0;
SaveFlag = 1;
B.SaveDir = SaveDir;
B.makeHist('ClustInfo','JumpRJ',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','X_SE',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','Y_SE',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','BG',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','Photons',PlotFlag,SaveFlag);

%Thresholding and reconstruction the image from the MCMC-output.
SZ = size(sequence,1);
Zoom = 20;
ImMCMC = threshRJ(B.ClustInfo,SZ,Zoom); 
RJ.produceIm(ImMCMC,fullfile(SaveDir,'MCMC'));
