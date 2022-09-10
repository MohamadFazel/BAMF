%Circle_demo is a demonstartion of RJ-class for one of the circle datasets
%in supplementary Fig. 3.
SaveDir = 'Circle';
addpath(fullfile('..','ptx'))
addpath(fullfile('..','mex'))
addpath('..')
%Loading data
load(fullfile(SaveDir,'CircleData.mat'))
%Initializing the class.
B=RJ('NoGui');
%The PSF size for this data-set is 1.2 pixel.
B.RJStruct.PSF_Sigma = 1.2; 
%Gain and offset.
B.Camera_Offset=0;
B.Camera_Gain=1;
%setting the jumps size in X and intensity to 0.05 pixel and 5 photons,
%respectively.
B.RJStruct.XstdFg = 0.05;
B.RJStruct.IstdFg = 5;
B.RJStruct.XDrift = zeros(1,size(Data,3));
B.RJStruct.YDrift = zeros(1,size(Data,3));
B.Data=Data;
%running the code
B.analyzeData()

%Saving the posterior image and the histogram of some of the parameters.
RJ.produceIm(B.PImage(41:120,41:120),fullfile(SaveDir,'Posterior'));
PlotFlag = 0;
SaveFlag = 1;
B.SaveDir = SaveDir;
B.makeHist('ClustInfo','JumpRJ',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','X_SE',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','Y_SE',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','BG',PlotFlag,SaveFlag);
B.makeHist('ClustInfo','Photons',PlotFlag,SaveFlag);

%Thresholding and reconstruction the image from the MCMC-output.
SZ = size(Data,1);
Zoom = 20;
ImMCMC = threshRJ(B.ClustInfo,SZ,Zoom); 
RJ.produceIm(ImMCMC(41:120,41:120),fullfile(SaveDir,'MCMC'));
