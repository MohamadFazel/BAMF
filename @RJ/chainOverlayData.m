function chainOverlayData(obj,FrameN,SigOrBg,Zoom)
%chainOverlayData() displays the ROI with the plot of the chain on top of it.
%   
% INPUTS:
%   obj:     The class object containing the chain of the accepted jumps.
%   FrameN:  The number of the frame that the user is intrested in.
%   SigOrBg: If the string is 'Sig', the signal chin will be plotted and
%            'Bg' will plot the structured background chain.
%   Zoom:    The zoom factor. (pixels)
%
% OUTPUTS:
%   No outputs.
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
if nargin<4
    Zoom=1;
end
if nargin<3
   SigOrBg = 'Sig'; 
end
Data = obj.Data;
BoxSize = obj.BoxSize;
Overlap = obj.Overlap;
ChainCell = obj.Chain;
ImSize=size(Data(:,:,1));
[~,ROI_Struct]=RJ.makeSubRegions(Data,BoxSize,Overlap);
MFrame = RJ.blockresample(Data(:,:,FrameN),Zoom);
Iters = find(ROI_Struct.FrameNum==FrameN);
dipshow(MFrame);
hold(gca,'on');
for mm = 1:ROI_Struct.Xsub*ROI_Struct.Ysub
    Ind = ROI_Struct.ROI_Ind(Iters(mm));
    ix = floor((Ind-1)/ROI_Struct.Xsub)+1;
    iy = Ind - (ix-1)*ROI_Struct.Xsub;
    sy = BoxSize*(iy-1)-Overlap;
    Lsy = BoxSize*(iy-1);
    Ley = BoxSize*iy;
    sx = BoxSize*(ix-1)-Overlap;
    Lsx = BoxSize*(ix-1);
    Lex = BoxSize*ix;  
    if iy == 1 
        sy = 0;
    end

    if ix == 1
        sx = 0;
    end
    if iy == ROI_Struct.Ysub
        Ley = ImSize(1);
    end

    if ix == ROI_Struct.Xsub
        Lex = ImSize(2);
    end
    Chain = ChainCell{Iters(mm)};
    LCh = length(Chain);
    switch SigOrBg
        case 'Sig'
            for nn = round(LCh/3):LCh
                SubChain=Chain(nn);
                if SubChain.N~=0
                   Signal = SubChain.Signal;
                   X = (sy+SubChain.X);
                   Y = (sx+SubChain.Y);
                   x=X;
                   y=Y;
                   X=X(x>Lsx&x<Lex&y>Lsy&y<Ley);
                   Y=Y(x>Lsx&x<Lex&y>Lsy&y<Ley);
                   Signal=Signal(x>Lsx&x<Lex&y>Lsy&y<Ley);
                   plot((X(Signal==1))*Zoom,(Y(Signal==1))*Zoom,'g.')
                end
             end
         case 'Bg'
             for nn = round(LCh/3):LCh
                 SubChain=Chain(nn);
                 if SubChain.N ~=0
                    Signal = SubChain.Signal;
                    X = (sy+SubChain.X);
                    Y = (sx+SubChain.Y);
                    x=X;
                    y=Y;
                    X=X(x>Lsx&x<Lex&y>Lsy&y<Ley);
                    Y=Y(x>Lsx&x<Lex&y>Lsy&y<Ley);
                    Signal=Signal(x>Lsx&x<Lex&y>Lsy&y<Ley);
                    plot((X(Signal==0))*Zoom,(Y(Signal==0))*Zoom,'g.')
                 end
             end
         case 'SigBg'
              for nn = round(LCh/3):LCh
                  SubChain = Chain(nn);
                  if SubChain.N~=0
                     X = (sy+SubChain.X);
                     Y = (sx+SubChain.Y);
                     x=X;
                     y=Y;
                     X=X(x>Lsx&x<Lex&y>Lsy&y<Ley);
                     Y=Y(x>Lsx&x<Lex&y>Lsy&y<Ley);
                     plot((X)*Zoom,(Y)*Zoom,'g.')  
                  end
              end
    end
    
end
end
