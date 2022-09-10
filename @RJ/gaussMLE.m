% gaussMLE produces a Maximum Likelihood Estimate of Gaussian blob parameters
% defining a PSF (Point Spread Function) for a stack of N subregions (or ROIs =
% Regions Of Interest), each of extent SZ x SZ.  gaussMLE uses GPU (CUDA)
% functionality for speed.
%
% INPUTS:
%    Data:             data to be analyzed of size (SZ,SZ,N)
%    FitType:          type of fit: 'Basic', 'Sigma', 'Z', 'SigmaYX'
%    CameraType:       'CCD' or 'SCMOS'
%    Sigma:            size or initial estimate of Gaussian blob [pixels]
%    VarianceIm:       scaled Variance image for use in SCMOS fitting
%    BoxCorners:       ROI box corners used to locate ROIs within VarianceImage
%                      in SCMOS fitting [N x 2] given in the order (x, y)
%    ZStruct:          structure for Z fitting based on astigmatism
%       Ax, Ay:
%       Bx, By:
%       Gamma:
%       D:
%
% OUTPUTS:
%    Results           results structure with components:
%       Y, X:             PSF center coordinates for each subregion [pixels]
%       Y_SE, X_SE:       standard error corresponding to the square root of
%                         the CRLB (Cramer-Rao Lower Bound) [pixels]
%       Photons:          photons per subregion
%       Photons_SE:       photon standard error per subregion
%       BG, BG_SE:        background and standard error per subregion
%       LogLikelihood:    log likelihood estimates per subregion
%       Sigma, Sigma_SE:  (FitType == 'Basic' || FitType == 'Sigma')
%       Z, Z_SE:          (FitType == 'Z')
%       SigmaY:           (FitType == 'SigmaYX')
%       SigmaX:           (FitType == 'SigmaYX')
%       SigmaY_SE:        (FitType == 'SigmaYX')
%       SigmaX_SE:        (FitType == 'SigmaYX')
%    Statistics:       information on the fitting with components:
%       FitTime:          CPU seconds used for fitting
%       FitPerSecond:     number of fits per second
%       NKernelsCalls:    number of calls to the CUDA kernel
%       NFitsPerCall:     number of fits per call to the CUDA kernel
%
% REQUIRES:
%    For data generation:
%       Statistics Toolbox or DipImage Toolbox
%    For GPU processing:
%       Parallel Computing Toolbox
%       NVidia GPU
%       cuda_gaussMLEv2.cu
%       cuda_gaussMLEv2.ptx
%
% USAGE:
%    CCD Fitting of Y,X,Photons,BG,(Sigma)
%       [Results,Statistics]=gaussMLE(Data,FitType,CameraType,Sigma)
%    SCMOS Fitting of Y,X,Photons,BG,(Sigma)
%       [Results,Statistics]=gaussMLE(Data,FitType,CameraType,Sigma, ...
%                                     VarianceImage,BoxCorners)
%    CCD Fitting of Y,X,Photons,BG,Z
%       [Results,Statistics]=gaussMLE(Data,FitType,CameraType,Sigma,ZStruct)
%    SCMOS Fitting of Y,X,Photons,BG,Z
%       [Results,Statistics]=gaussMLE(Data,FitType,CameraType,Sigma, ...
%                                     VarianceImage,BoxCorners,ZStruct)
%
% CITATION:
%    Michael Wester (Lidke Lab, April 26, 2017)
function [Results, Statistics]=gaussMLE(varargin)
    
    if length(varargin)<4
        error('gaussMLE must have 4 to 7 inputs. See doc gaussMLE')
    elseif length(varargin)>7
        error('gaussMLE must have 4 to 7 inputs. See doc gaussMLE')
    end
    Data=single(varargin{1});
    FitType=varargin{2};
    CameraType=varargin{3};
    Sigma=varargin{4};
    switch CameraType
        case 'CCD'
            switch FitType
                case 'Z'
                    ZStruct=varargin{5};
            end
        case 'SCMOS'
            VarianceImage=varargin{5};
            %VarianceImage=NoiseImage.^2;
            BoxCorners=varargin{6};
            switch FitType
                case 'Z'
                    ZStruct=varargin{7};
            end
    end
    
    %Set up kernels and number of parameters
    %Number of iterations in the NR update step
    Iterations=[];
    switch FitType
        case 'Basic'
            NP=4;
            KernelID='_XYNB_';
            if isempty(Iterations);Iterations=20;end
        case 'Sigma'
            NP=5;
            KernelID='_XYNBS_';
            if isempty(Iterations);Iterations=20;end
        case 'Z'
            NP=5;
            KernelID='_XYNBZ_';
            if isempty(Iterations);Iterations=20;end
        case 'SigmaYX'
            NP=6;
            KernelID='_XYNBSXSY_';
            if isempty(Iterations);Iterations=20;end
    end
    switch CameraType
        case 'SCMOS'
            if size(VarianceImage,1)~=size(VarianceImage,2)
                error('gaussMLE: SCMOS fitting is only implemented for square Camera ROI')
            end
            KernelID=['_SCMOS' KernelID(2:end)];
    end
    k = parallel.gpu.CUDAKernel('cuda_gaussMLEv2.ptx','cuda_gaussMLEv2.cu',KernelID);
    
    N=size(Data,3);
    SZ=size(Data,1);
    G=gpuDevice();
    AM=G.AvailableMemory;
    BytesPerFloat=4;
    MemoryPerFit=BytesPerFloat*(SZ*SZ+NP+NP+1); %Data plus outputs
    MaxFits=floor(AM/MemoryPerFit/2); %Factor of 2 to be conservative
    
    BSZ=128; %Threads per block hard-coded into kernel
    NKernelsCalls = ceil(N/MaxFits);
    NFitsPerCall=min(MaxFits,N);
    
    %Setup output arrays
    Params_out=zeros(N,NP,'single');
    CRLB_out=zeros(N,NP,'single');
    LL_out=zeros(N,1,'single');
    %LL_out=zeros(N,NP,'single');
    
    tic
    for nn=1:NKernelsCalls
        StartIndex=(nn-1)*NFitsPerCall+1;
        EndIndex=min((nn+1)*NFitsPerCall,N); %Don't go past data size
        NFitsActual = EndIndex-StartIndex+1;
        SubData=Data(:,:,StartIndex:EndIndex);
        k.GridSize = [ceil(NFitsActual/BSZ) 1];
        k.ThreadBlockSize = [BSZ 1];
        d_Parameters=zeros(NFitsActual,NP,'single');
        d_CRLBs=zeros(NFitsActual,NP,'single');
        d_LogLikelihood=zeros(NFitsActual,1,'single');
        switch CameraType
            case 'CCD'
                switch FitType
                    case {'Basic','Sigma','SigmaYX'}
                        [P, CRLB,LL] = feval(k,SubData,mean(Sigma),SZ,Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                    case 'Z'
                        [P, CRLB,LL] = feval(k,SubData,Sigma(1),...
                            ZStruct.Ax,ZStruct.Ay, ZStruct.Bx, ZStruct.By,ZStruct.Gamma,ZStruct.D,Sigma(2),...
                            SZ,Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                        otherwise
                        error('gaussMLE: Unknown fit type: %s',FitType)
                end
            case 'SCMOS'
                switch FitType
                    case {'Basic','Sigma','SigmaYX'}
                        [P, CRLB,LL] = feval(k,SubData,BoxCorners,VarianceImage,...
                            mean(Sigma),SZ,size(VarianceImage,1),Iterations,...
                            d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                    case 'Z'
                        Z0=zeros(NFitsActual,1,'single'); 
                        [P, CRLB,LL] = feval(k,SubData,BoxCorners,VarianceImage,Z0,...
                            Sigma(1),ZStruct.Ax,ZStruct.Ay, ZStruct.Bx,ZStruct.By,...
                            ZStruct.Gamma,ZStruct.D,Sigma(2),...
                            SZ,size(VarianceImage,1),Iterations,d_Parameters,d_CRLBs,d_LogLikelihood,NFitsActual);
                    otherwise
                        error('gaussMLE: Unknown fit type: %s',FitType)
                end
        end
        Params_out(StartIndex:EndIndex,:)=gather(P);
        CRLB_out(StartIndex:EndIndex,:)=gather(CRLB);
        LL_out(StartIndex:EndIndex)=gather(LL);
    end
    
    Results.Y=Params_out(:,1);
    Results.X=Params_out(:,2);
    Results.Photons=Params_out(:,3);
    Results.BG=Params_out(:,4);
    Results.Y_SE=sqrt(CRLB_out(:,1));
    Results.X_SE=sqrt(CRLB_out(:,2));
    Results.Photons_SE=sqrt(CRLB_out(:,3));
    Results.BG_SE=sqrt(CRLB_out(:,4));
    switch FitType
        case 'Basic'
        case 'Sigma'
            Results.Sigma=Params_out(:,5);
            Results.Sigma_SE=sqrt(CRLB_out(:,5));
        case 'Z'
            Results.Z=Params_out(:,5);
            Results.Z_SE=sqrt(CRLB_out(:,5));
        case 'SigmaYX'
            Results.SigmaY=Params_out(:,5);
            Results.SigmaX=Params_out(:,6);
            Results.SigmaY_SE=sqrt(CRLB_out(:,5));
            Results.SigmaX_SE=sqrt(CRLB_out(:,6));
    end
    Results.LogLikelihood=LL_out;
     
    Statistics.FitTime=toc;
    Statistics.FitPerSecond=N/Statistics.FitTime;
    Statistics.NKernelsCalls = NKernelsCalls;
    Statistics.NFitsPerCall=NFitsPerCall;
end
