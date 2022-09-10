classdef RJ < handle
%RJ is the main class of the BAMF algorithm.
%  RJ is a calss to implement RJMCMC analysis code to reconstruct
%  super-resolution images with fine detailes and high resolution. The
%  class is consist of methods and static methods. A gui method is written
%  for user friendly interface. The following is a somewhat detailed
%  explanation on how the class works. The main method is called
%  analyzeData() and the other methods is being called from inside that.
%  First the gain and offset correction (gainCorrection() method) are
%  implemented. Then the findPriors() methods call a fast single emitter 
%  code to calculate the intensity and offset background priors. Next, the
%  data is sent to all the workers available on the machine using a parfor
%  loop. Inside the parfor loop the sequence is split into a number of
%  subsequences and each worker works on one of them. Then each frame is
%  split into subregions that are processed independently from each other.
%  After that, findPSF_SMA() method uses the subregions to find the numeric
%  PSF for different regions of the data. Then, the RJMCMC algorithm is
%  implemented inside the mex-function. The returned results are a chain of
%  all the jumps and the posterior image. The posterior images are put
%  together to get the posterior image for the entire dataset. The chain is
%  also used to extract the most likely model (MAPN), which is done inside
%  the method findMap(). The MAPN results are saved in a structure called 
%  ClustInfo. At the end, thses outcomes are retrieved from the workers and
%  put to gether to get the final results. 
%  There are some methods for diagnosys after the algorithm is done, for
%  most of them the user need to save the chain. Note that the chain takes
%  a large space on the memory and for large datastes it may not be saved.
%  A short description of these diagnosys methods are given here. For a
%  detailed description see the help for each method.
%  makeHist(): make histogram plots of some of the different jumps'
%  acceptance ratios, intensities and the standard deviations of the found
%  positions, intensities and background parameters.
%  makeModel(): returns a movie of the overlay of a ROI, selected by the
%  user, and the states of the chain. To use this method you need the
%  chain.
%  plotChain(): returns an image of the chain and the user can choose to
%  see different models that exist in the chain.
%  chainOverlayData(): makes a plot of all the accepted posiyions on top of the raw
%  data for a single ROI, selected by user.
%  chainAnimation(): shows an animation of the chain, including the found
%  positions and intensities and the type of attempted jump.
%
% REQUIRES:
%  MATLAB 2016 or higher versions.
%  Some methods have additional requirements.             
%
% USAGE:
%  A=RJ;
%  DataPath=which('Cell01_VeryDence01.mat'); 
%  DataStruct = load(DataPath);
%  A.Data = DataStruct.sequence;
%  A.Camera_Gain = 11;
%  A.Camera_Offset = 100;
%  A.RJStruct.PSF_Sigma = 1.3;
%  A.analyzeData();
%  dipshow(A.PImage)
%  colormap('hot')
%  dipshow(A.PBg)
%  colormap('hot')
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%
    properties
        RJStruct      %Structure containing Priors, Number of the jumps etc.  
        Data          %The raw supre-resolution data.
        ROI_Data      %The sub-regions (output of makeSubregions()).
        PImage        %The final reconstruction of the found signal.
        PBg           %The final reconstruction of the found background.
        TimeSeriesIm; %The resulted posteriors for the ROIs are saved here.
        TimeSeriesBg; %The resulted posterios for background is saved here.
        ClustInfo=[]; %The MAP results are saved in this structure.
        Chain;        %Cell containing chains of accepted jumps.
        PChain;       %Cell containing chains of the proposed jumps.
        JumpStat;     %The jump statistics are saved in this structure.
        Camera_Gain = 11;   %The camera gain. 
        Camera_Offset = 100; %The camera offset.
        Camera_Var    %The read-out noise for sCMOS camera. 
        BoxSize = 16; %The size of the ROIs (without overlapping regions).
        Overlap = 3;  %The size of the overlapping regions (in pixel).
        TimeSeries=0; %zero means do not save the posterior image for the ROIs.
        KeepChain = 0;%Because the size of the chain can be very large, you 
                      %don't want to save it always. 1 saves the chain 
                      %and 0 throws it away. To use the diagnosis methods 
                      %you must save the chain.
        UsePPToolbox=1; %Use Parallel processing toolbox. Default=1 if present  
        SMAStruct = []; %Parameters inside this structures are used by the
                        %single emitter code.
        SampledPSF=[];  %Found numberical PSFs are saved here.
        FindPSF = 0;    %Find numeric PSF and use it for analsys or not.
        FindPriors = 1; %1 means calculate the priors and 0 means use the 
                        %provided prior.
        DataDir;        %Data directory
        Threshold = []; %Thresholding parameters are inside this structure
        FrameConnect = [];%The parameters required for Frame-connect code 
                          %are stored in this structure.
        DriftXDir;      %Direction to the drift files, used in gui.
        DriftYDir;      %Direction to the drift files, used in the gui.
        FileListDriftX; %List of the drift files, used in the gui.
        FileListDriftY; %List of the drift files, used in the gui.
        TotClust = [];    
        SMD = [];
        SMD_combined = [];
        FileList;       %List of all the datasets that need to be analyzed.
        SaveDir;        %The directory that the results are saved.
        ThreshFlag = 1; %1 means do thresholding, 0 not, used in the gui.
        FCFlag = 0;     %1 means do frame connect, 0 not, used in the gui.
        DriftFlag = 0;  %1 means do drift correction, 0 not, used in the gui.
    end
    
    methods
        function obj=RJ(GuiFlag)
            %Create RJ object
            obj.RJStruct=RJ.genRJStruct();
            obj.RJStruct.XDrift = zeros(size(obj.Data,3),1);
            obj.RJStruct.YDrift = zeros(size(obj.Data,3),1);
            obj.Threshold.Max_XY_SE = 0.25;
            obj.Threshold.Min_XY_SE = 0.001;
            obj.Threshold.Max_Photons_SE = 50;
            obj.Threshold.MinPhotons = 200;
            obj.FrameConnect.MinPValue = 0.01;
            obj.FrameConnect.MaxFrame_gap = 10;
            obj.FrameConnect.MaxDistance = 4;
            obj.FrameConnect.LOS = 0.01;
            obj.SMAStruct.MeanPhotons = 800;
            obj.SMAStruct.MinPhotons = 200;
            obj.SMAStruct.MinPValue = 0.01;
            %Use PPT if present 
            obj.UsePPToolbox=license('test', 'Distrib_Computing_Toolbox');
            if nargin == 0
                obj.gui;
            end
        end
        [A,B]=makeResultsRJ(obj)
        makeHist(obj,StructString,Variable,PlotFlag,SaveFlag)
        Im = makeGaussIm(obj,ClustInfoT);
        [Model,Fg,Bg]=makeModel(obj,ROINum,PlotFlag,PSFSigma)
        [Out1,Out2]=makeChainHists(obj,Variable,Index,PlotFlag)
        [ImOf]=plotChain(obj,ROINum,N,PlotFlag,SigBg)
        chainOverlayData(obj,FrameN,SigOrBg,Zoom)
        chainAnimation(obj,ROINum,X,Y,StartState)
        obj=makeClustInfo(obj)
        processData(obj,TimeID);
        obj=assembleResults(obj);
        function gainCorrection(obj)
            % Gain and offset correction. 
            if isempty(obj.Data)
               error('Data has to be assigned to the property Data.'); 
            end
            if isempty(obj.Camera_Offset)
               error('Camera_Offset is empty.'); 
            end
            if isempty(obj.Camera_Gain)
               error('Camera_Gain is empty.'); 
            end
            obj.Data = (single(obj.Data)-obj.Camera_Offset)./obj.Camera_Gain;
        end

        function analyzeData(obj)
           %Setup and process sub regions
           %analyzeData() implements makeSubRegions() to break the input
           %data into ROIs and then find the essential priors calling
           %findPriors(). Then it calls analyzeROIs() to analyze the ROIs. 
           %Finally, removeOverlapping() and stitchROIs() are called to
           %make the final reconstruction.
           %fprintf('analyzeData\n');
           %profile on 
           obj.gainCorrection()
           if isempty(obj.RJStruct.PSF_Sigma)
              fprintf('RJStruct.PSF_Sigma is empty so we assigned 1.3 (pixels) to that.\n');
              obj.RJStruct.PSF_Sigma = 1.3;
           end
           if isempty(obj.RJStruct.P_dP)
              fprintf('RJStruct.P_dP is empty so we assign 50 (Photons) to that.');
              obj.RJStruct.P_dP = 10;
           end
           
           if obj.FindPriors
               if isempty(obj.Camera_Var) || isscalar(obj.Camera_Var)
                   obj = obj.findPriors(1,'Yes');
               else
                   obj = obj.findPriors(1,'Yes');
               end
           elseif isempty(obj.RJStruct.P_EmPhotons) || isempty(obj.RJStruct.P_BgEmPhotons) || isempty(obj.RJStruct.P_Offset)
               error('Priors Must be given.');
           end
            [obj.ROI_Data,ROIstruct,~]=RJ.makeSubRegions(obj.Data,obj.BoxSize,obj.Overlap,obj.Camera_Var);
            if ~isfield(obj.SMAStruct,'BoxSize')
               obj.SMAStruct.BoxSize = 8; 
            end
            if obj.FindPSF == 1
                obj=obj.findPSF_SMA(ROIstruct);
            end
            DataPar = obj.Data;
            BoxSizePar = obj.BoxSize;
            OverlapPar = obj.Overlap;
            SampledPSFpar = obj.SampledPSF;
            FindPSFpar = obj.FindPSF;
            Camera_VarPar = obj.Camera_Var;
            TimeSe = obj.TimeSeries;
            RJstrc = obj.RJStruct;
            Keep = obj.KeepChain;
            PImageT=0;
            PBgT=0;
            ClustInfoTemp = [];
            ChainCell = [];
            PChainCell = [];
            StatCell = [];
            TempSeriesIm=[];
            TempSeriesBg=[];
            Time = tic();
            if obj.UsePPToolbox 
                MachineInfo=parcluster();
                NumWorkers=MachineInfo.NumWorkers;
                parpool('local',NumWorkers);
                TempFrame = floor(size(DataPar,3)/NumWorkers);
                parfor pp = 1:NumWorkers
                    EndFrame = pp*TempFrame;
                    if pp == NumWorkers
                       EndFrame = size(DataPar,3);
                    end
                    [ROI_DataT,ROI_Struct,ROI_Var]=RJ.makeSubRegions(DataPar(:,:,(pp-1)*TempFrame+1:EndFrame),BoxSizePar,OverlapPar,Camera_VarPar);
                    ROI_Struct.TempFrame = TempFrame;
                    [PIm,PBgT1,ClustInfoT,ChainCellT,StatCellT,PChainCellT,TSeriesIm,TSeriesBg]=...
                        RJ.analyzeROIs(ROI_Struct,RJstrc,ROI_DataT,ROI_Var,SampledPSFpar,FindPSFpar,Keep,TimeSe,pp);
                    [PROIs_Signal]=RJ.removeOverlapping(ROI_Struct,OverlapPar,RJstrc.SRZoom,PIm);
                    [PROIs_Bg]=RJ.removeOverlapping(ROI_Struct,OverlapPar,RJstrc.SRZoom,PBgT1);
                    [PImageTemp]=RJ.stitchROIs(PROIs_Signal,ROI_Struct,OverlapPar,[size(DataPar(:,:,1)),1],BoxSizePar,RJstrc.SRZoom);
                    [PBgTemp]=RJ.stitchROIs(PROIs_Bg,ROI_Struct,OverlapPar,[size(DataPar(:,:,1)),1],BoxSizePar,RJstrc.SRZoom);
                    PImageT = PImageT + PImageTemp;
                    PBgT = PBgT + PBgTemp;
                    ClustInfoTemp = [ClustInfoTemp,ClustInfoT];
                    StatCell = [StatCell,StatCellT];
                    ChainCell = [ChainCell,ChainCellT];
                    PChainCell = [PChainCell,PChainCellT];
                    TempSeriesBg = [TempSeriesBg;TSeriesBg];
                    TempSeriesIm = [TempSeriesIm;TSeriesIm];
                end
                delete(gcp('nocreate'));
            else
                NumWorkers = 1;
                pp=1;
                TempFrame = floor(size(DataPar,3)/NumWorkers);
                [ROI_DataT,ROI_Struct,ROI_Var]=RJ.makeSubRegions(DataPar(:,:,(pp-1)*TempFrame+1:pp*TempFrame),BoxSizePar,OverlapPar,Camera_VarPar);
                ROI_Struct.TempFrame = TempFrame;
                [PIm,PBgT1,ClustInfoT,ChainCellT,StatCellT,PChainCellT,TSeriesIm,TSeriesBg]=...
                     RJ.analyzeROIs(ROI_Struct,RJstrc,ROI_DataT,ROI_Var,SampledPSFpar,FindPSFpar,Keep,TimeSe,pp);
                [PROIs_Signal]=RJ.removeOverlapping(ROI_Struct,OverlapPar,RJstrc.SRZoom,PIm);
                [PROIs_Bg]=RJ.removeOverlapping(ROI_Struct,OverlapPar,RJstrc.SRZoom,PBgT1);
                [PImageTemp]=RJ.stitchROIs(PROIs_Signal,ROI_Struct,OverlapPar,[size(DataPar(:,:,1)),1],BoxSizePar,RJstrc.SRZoom);
                [PBgTemp]=RJ.stitchROIs(PROIs_Bg,ROI_Struct,OverlapPar,[size(DataPar(:,:,1)),1],BoxSizePar,RJstrc.SRZoom);
                PImageT = PImageT + PImageTemp;
                PBgT = PBgT + PBgTemp;
                ClustInfoTemp = [ClustInfoTemp,ClustInfoT];
                StatCell = [StatCell,StatCellT];
                ChainCell = [ChainCell,ChainCellT];
                PChainCell = [PChainCell,PChainCellT];
                TempSeriesBg = [TempSeriesBg;TSeriesBg];
                TempSeriesIm = [TempSeriesIm;TSeriesIm];
            end
            obj.Chain = ChainCell;
            obj.PChain = PChainCell;
            obj.JumpStat = StatCell;
            obj.TimeSeriesIm = TempSeriesIm;
            obj.TimeSeriesBg = TempSeriesBg;
            obj = obj.makeClustInfo;

            for ii = 1:NumWorkers
                obj.ClustInfo.X = [obj.ClustInfo.X,ClustInfoTemp(ii).X];
                obj.ClustInfo.Y = [obj.ClustInfo.Y,ClustInfoTemp(ii).Y];
                obj.ClustInfo.I = [obj.ClustInfo.I,ClustInfoTemp(ii).I];
                obj.ClustInfo.X_SE = [obj.ClustInfo.X_SE,ClustInfoTemp(ii).X_SE];
                obj.ClustInfo.Y_SE = [obj.ClustInfo.Y_SE,ClustInfoTemp(ii).Y_SE];
                obj.ClustInfo.I_SE = [obj.ClustInfo.I_SE,ClustInfoTemp(ii).I_SE];
                obj.ClustInfo.X_kmeans = [obj.ClustInfo.X_kmeans,ClustInfoTemp(ii).X_kmeans];
                obj.ClustInfo.Y_kmeans = [obj.ClustInfo.Y_kmeans,ClustInfoTemp(ii).Y_kmeans];
                obj.ClustInfo.I_kmeans = [obj.ClustInfo.I_kmeans,ClustInfoTemp(ii).I_kmeans];
                obj.ClustInfo.X_SE_kmeans = [obj.ClustInfo.X_SE_kmeans,ClustInfoTemp(ii).X_SE_kmeans];
                obj.ClustInfo.Y_SE_kmeans = [obj.ClustInfo.Y_SE_kmeans,ClustInfoTemp(ii).Y_SE_kmeans];
                obj.ClustInfo.I_SE_kmeans = [obj.ClustInfo.I_SE_kmeans,ClustInfoTemp(ii).I_SE_kmeans];
                obj.ClustInfo.Xbg = [obj.ClustInfo.Xbg,ClustInfoTemp(ii).Xbg];
                obj.ClustInfo.Ybg = [obj.ClustInfo.Ybg,ClustInfoTemp(ii).Ybg];
                obj.ClustInfo.Ibg = [obj.ClustInfo.Ibg,ClustInfoTemp(ii).Ibg];
                obj.ClustInfo.Xbg_SE = [obj.ClustInfo.Xbg_SE,ClustInfoTemp(ii).Xbg_SE];
                obj.ClustInfo.Ybg_SE = [obj.ClustInfo.Ybg_SE,ClustInfoTemp(ii).Ybg_SE];
                obj.ClustInfo.Ibg_SE = [obj.ClustInfo.Ibg_SE,ClustInfoTemp(ii).Ibg_SE];
                obj.ClustInfo.FrameNum = [obj.ClustInfo.FrameNum,ClustInfoTemp(ii).FrameNum];
                obj.ClustInfo.FrameNumBg = [obj.ClustInfo.FrameNumBg,ClustInfoTemp(ii).FrameNumBg];
                obj.ClustInfo.FrameNum_kmeans = [obj.ClustInfo.FrameNum_kmeans,ClustInfoTemp(ii).FrameNum_kmeans];
                obj.ClustInfo.ROInum = [obj.ClustInfo.ROInum,ClustInfoTemp(ii).ROInum];
                obj.ClustInfo.Bg = [obj.ClustInfo.Bg,ClustInfoTemp(ii).Bg];
                obj.ClustInfo.ABg = [obj.ClustInfo.ABg,ClustInfoTemp(ii).ABg];
                obj.ClustInfo.BBg = [obj.ClustInfo.BBg,ClustInfoTemp(ii).BBg];
                obj.ClustInfo.ROIcornerX = [obj.ClustInfo.ROIcornerX,ClustInfoTemp(ii).ROIcornerX];
                obj.ClustInfo.ROIcornerY = [obj.ClustInfo.ROIcornerY,ClustInfoTemp(ii).ROIcornerY];
                obj.ClustInfo.JumpXRJ = [obj.ClustInfo.JumpXRJ,ClustInfoTemp(ii).JumpXRJ];
                obj.ClustInfo.JumpIRJ = [obj.ClustInfo.JumpIRJ,ClustInfoTemp(ii).JumpIRJ];
                obj.ClustInfo.JumpBgRJ = [obj.ClustInfo.JumpBgRJ,ClustInfoTemp(ii).JumpBgRJ];
                obj.ClustInfo.JumpXFg = [obj.ClustInfo.JumpXFg,ClustInfoTemp(ii).JumpXFg];
                obj.ClustInfo.JumpIFg = [obj.ClustInfo.JumpIFg,ClustInfoTemp(ii).JumpIFg];
                obj.ClustInfo.JumpXBg = [obj.ClustInfo.JumpXBg,ClustInfoTemp(ii).JumpXBg];
                obj.ClustInfo.JumpIBg = [obj.ClustInfo.JumpIBg,ClustInfoTemp(ii).JumpIBg];
                obj.ClustInfo.SplitRJ = [obj.ClustInfo.SplitRJ,ClustInfoTemp(ii).SplitRJ];
                obj.ClustInfo.MergeRJ = [obj.ClustInfo.MergeRJ,ClustInfoTemp(ii).MergeRJ];
                obj.ClustInfo.BirthRJ = [obj.ClustInfo.BirthRJ,ClustInfoTemp(ii).BirthRJ];
                obj.ClustInfo.DeathRJ = [obj.ClustInfo.DeathRJ,ClustInfoTemp(ii).DeathRJ];
                obj.ClustInfo.JumpBgRJ = [obj.ClustInfo.JumpBgRJ,ClustInfoTemp(ii).JumpBgRJ];
                obj.ClustInfo.JumpMC = [obj.ClustInfo.JumpMC,ClustInfoTemp(ii).JumpMC];
            end
            obj.PImage = PImageT;
            obj.PBg = PBgT;
            T = toc(Time);
            fprintf('It took %g seconds per ROI for RJMCMC to analyze this data set.\n',T/size(obj.ROI_Data,3));
        end
        
    end
    methods (Static)
        %static methods
        [PIm,PBgT1,ClustInfoT,ChainCellT,StatCellT,PChainCellT,TSeriesIm,TSeriesBg]=...
            analyzeROIs(ROI_Struct,RJstrc,ROI_DataT,ROI_Var,SampledPSF,OptModel,Keep,TimeSe,pp);
        RJStruct=genRJStruct(); %set up the RJStruct.
        [ROI_Data,ROI_Struct,ROI_Var]=makeSubRegions(Data,BoxSize,Overlap,Camera_Var); 
        SampledPSFm=findPSF(ROIStack,Gain,Offset,PSFSigma,BoxSize);
        [Chain,Stat,PImage,PBackg,ProposedChain]=rjmcmc(Data,RJStruct,Ysize,Xsize,DriftY,DriftX,CMOSVar,SampledPSF_in,OptModel_in) %setup and call c-code
        [Out]=removeOverlapping(ROI_Struct,Overlap,Zoom,PIm)
        [FinalIm]=stitchROIs(ROIs,ROI_Struct,Overlap,ImSize,BoxSize,Zoom); 
        [ImOf]=makeMapN_ROI(Chain,ROISize,Zoom);                                %not in use any more
        [OutIm,FoundMap] = hist_MapN(Film, MinWeight, MaxStd, CutOff);          %not in use any more
        [SMD, ROIStack, Res]=findROI(Image,Sigma1,Sigma2,Boxsize,Minval);
        [Results, Statistics]=gaussMLE(varargin);
        SMA=threhsSMA(SMD,SMAStruct,BoxSize);
        [OutIm_int]=convertTIFsequence(Path,Name);
        produceIm(Im,Name);
        MagnifiedIm=blockresample(Im,Zoom);
        [SMD,MostFrequent]=findMap(Chain,FB)
        ROI = makeROIData(Data,Begin,ROISize);
        [FoundMAP,SMDkmeans,SMDbg,Chain,Stat]=mcmc(ThisSubRegion,RJstrc,Xsize,Ysize,DriftY,...
            DriftX,ThisSubVar,SampledPSF,OptModel,ROI_Struct,SMD,SMDbg,ROIid,FrameID,Xstd_cutoff,Istd_cutoff);
        [Data]=loadData(DataName,DataID);
         twoPoints(Separation,Photons,PSF,SZ,Bg,NFrames,HWMHprior,WeightOfSecond)
    end
    
end
