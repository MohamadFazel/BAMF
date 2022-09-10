function SampledPSFm=findPSF(ROIStack,Gain,Offset,PSFSigma,BoxSize)
%findPSF() uses Fourier transforms of the isolated emitters to calculate PSF
%   PSF is created by shifting and averaging over more than 100 high 
%   signal, found single emitters. A 4 sub-sampled PSF is created by 
%   padding the Fourier transform. 
%
% INPUTS:
%   ROIStack:   Stack of ROIs containing isolated emitters.
%   Gain:       Camera gain.
%   Offset:     Camera offset.
%   PSFSigma:   PSF size (pixels).
%   BoxSize:    Size of the returnd PSF stack.
%
% OUTPUT:
%   SampledPSFm: The sampled PSF.
%
% REQUIRES:
%   MATLAB 2014a or higher versions.
%   NVIDIA GPU
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel and Keith A. Lidke (Lidke Lab 2018)
%

ROIStack = single((ROIStack-Offset)/Gain);
[Results,Stats]=SMA_Core.gaussMLE(ROIStack,'Sigma','CCD',PSFSigma);

R = ROIStack;
X = Results.X;
Y = Results.Y;
% Now shift them on top of each other in fourier space

N=length(X);
Fm =0;
[XGm,YGm]=meshgrid((-ceil((BoxSize-1)/2):floor((BoxSize-1)/2)),(-ceil((BoxSize-1)/2):floor((BoxSize-1)/2)));

for nn=1:min([N,1000])
    FIMm = fftshift(fft2(fftshift(R(:,:,nn)/10)));
    %Moving the centres of all the PSFs to the same location by adding
    %phase to their Fourier transform.
    FIMm=FIMm.*exp(2*pi*1i*XGm*(X(nn)-BoxSize/2-1/2)/BoxSize).*exp(2*pi*1i*YGm*(Y(nn)-BoxSize/2-1/2)/BoxSize);
    Fm = Fm+FIMm;
end

% Oversample
OverSampling=4;
FitBoxSize=16;
%band limit
BLm = sqrt(XGm.^2+YGm.^2)<BoxSize/2;
%pad in F space
FPadm=wextend(2,'zpd',single(Fm),12);

%real space sampled and normalized
PSF_Sampledm=abs(fftshift(ifft2(FPadm)));

Am = PSF_Sampledm(1,:);
Bm = PSF_Sampledm(end,:);
Cm = PSF_Sampledm(:,end);
Dm = PSF_Sampledm(:,1);
Offsetm=mean( [mean(Am(:)), mean(Bm(:)),mean(Cm(:)),mean(Dm(:))]);

PSF_Sampledm=PSF_Sampledm-Offsetm;
PSF_Sampledm(PSF_Sampledm<0)=0;

PSF_Sampledm=PSF_Sampledm/max(PSF_Sampledm(:));

%pad to 4x fitting region size (16)x oversampling;

PSF_Sampledm=wextend(2,'zpd',PSF_Sampledm,112);
PSF_Sampledm=single(PSF_Sampledm);
%we could intgrate with convolution here...
PSFArraym=zeros(4*FitBoxSize,4*FitBoxSize,OverSampling,OverSampling);
%break that into 16 sub images

[Ys,Xs]=size(PSF_Sampledm);
[Xg,Yg]=meshgrid((-ceil(Xs/2):floor((Xs-1)/2)),(-ceil(Ys/2):floor((Ys-1)/2)));
XGm=OverSampling-Xg/OverSampling;XGm=(XGm-floor(XGm))*OverSampling;
YGm=OverSampling-Yg/OverSampling;YGm=(YGm-floor(YGm))*OverSampling;

PSFTempm = zeros(4*FitBoxSize*4*FitBoxSize,1,'single');
clear PSFArraym
for ii=0:OverSampling-1
    for jj=0:OverSampling-1
        Maskm=XGm==ii;
        Maskm=Maskm&YGm==jj;
        Maskm=single(Maskm);
        Cnt=0;
        for kk=0:size(PSF_Sampledm,1)-1
            for ll=0:size(PSF_Sampledm,2)-1
                if Maskm(ll+1,kk+1)
                    PSFTempm(Cnt+1)=PSF_Sampledm(ll+1,kk+1);
                    Cnt=Cnt+1;
                end
            end
        end
        
        PSFTempm=reshape(PSFTempm,[4*FitBoxSize,4*FitBoxSize]);
        [SZ1,SZ2] = size(PSFTempm);
        if jj>0
            PSFTempm(round(SZ1/3)+1:round(2*SZ1/3)+1,:)=single(PSFTempm(round(SZ1/3):round(2*SZ1/3),:));
        end
        if ii>0
            PSFTempm(:,round(SZ2/3)+1:round(2*SZ2/3)+1)=single(PSFTempm(:,round(SZ2/3):round(2*SZ2/3)));
        end
        
        %subtract offset
        PSFTempm=PSFTempm/sum(PSFTempm(:));
        
        %normalization
        PSFArraym(:,:,jj+1,ii+1)=PSFTempm;
        
    end
end
SampledPSFm=single(PSFArraym);
end