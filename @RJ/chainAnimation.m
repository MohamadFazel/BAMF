function chainAnimation(obj,ROINum,x,y,StartState)
%chainAnimation() is a diagnosis function to go through all the states.
%   When this function s called a gui will be popped up that lets the user
%   to look at different proposed and accepted jumps in the chain.
%   Different parameters such as the emitters parameters, likelihood ratio,
%   posterior ratio, etc., will be displayed on the gui. Since the true
%   locations of the emitters are required, this function can only be used
%   for simulated data.
%
% INPUTS:
%   obj:        The class object containing the chain of the proposed and
%               accepted jumps.
%   ROINum:     Number of the ROI that the user is interested in.
%   x:          True X-locations of the emitters (pixels).
%   y:          True Y-locations of the emitters (pixels).
%   StartState: The number of the jump to start with.
%
% OUTPUTS:
%   No output.
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

if nargin<5
   StartState=1; 
end
ChainROI = obj.Chain{ROINum};
PChainROI = obj.PChain{ROINum};
BoxSize = obj.BoxSize+2*obj.Overlap;
L = length(ChainROI)-1;
Loop = 1;
guiFig=figure();
%guiFig.Position = [530 400 1250 500];
guiFig.Position = [530 80 1300 880];
hFig1Panel = uipanel('Parent',guiFig,'Title','Proposed','Position',[0.011 0.44 0.315 0.48]);
hAxis1 = axes('Parent',hFig1Panel,'Units','pixels','XAxisLocatio','top',...
    'Ydir','reverse','Position',[45,15,350,350]);
hFig2Panel = uipanel('Parent',guiFig,'Title','Accepted','Position',[0.335 0.44 0.315 0.48]);
hAxis2 = axes('Parent',hFig2Panel,'Units','pixels','XAxisLocation','top',...
    'Ydir','reverse','Position',[45,15,350,350]);
hData1Panel = uipanel('Parent',guiFig,'Title','Proposed','Position',[0.83 0.44 0.165 0.43]);  
hData2Panel = uipanel('Parent',guiFig,'Title','Accepted','Position',[0.658 0.44 0.165 0.43]);  

hData1PanelSub = uipanel('Parent',hData1Panel,'Position',[0 -1 0.9 1.93]); 
hData2PanelSub = uipanel('Parent',hData2Panel,'Position',[0 -1 0.9 1.93]); 

hSlider1 = uicontrol('Style','Slider','Parent',hData1Panel,'Units','normalized',...
    'Position',[0.9,0,0.1,0.93],'Value',1,'Callback',{@slider1,hData1PanelSub});
hSlider2 = uicontrol('Style','Slider','Parent',hData2Panel,'Units','normalized',...
    'Position',[0.9,0,0.1,0.93],'Value',1,'Callback',{@slider2,hData2PanelSub});

htext1P = uicontrol('Parent',hData1Panel,'Style','text','String','Photons','Position',[5,340,50,20]);
htext1X=  uicontrol('Parent',hData1Panel,'Style','text','String','X','Position',[60,340,50,20]);
htext1Y=  uicontrol('Parent',hData1Panel,'Style','text','String','Y','Position',[115,340,50,20]);
htext2P = uicontrol('Parent',hData2Panel,'Style','text','String','Photons','Position',[5,340,50,20]);
htext2X=  uicontrol('Parent',hData2Panel,'Style','text','String','X','Position',[60,340,50,20]);
htext2Y=  uicontrol('Parent',hData2Panel,'Style','text','String','Y','Position',[115,340,50,20]);

htext3P = uicontrol('Parent',hData1PanelSub,'Style','text','String','','Position',[5,5,50,690]);
htext3X = uicontrol('Parent',hData1PanelSub,'Style','text','String','','Position',[60,5,50,690]);
htext3Y = uicontrol('Parent',hData1PanelSub,'Style','text','String','','Position',[115,5,50,690]);
htext4P = uicontrol('Parent',hData2PanelSub,'Style','text','String','','Position',[5,5,50,690]);
htext4X = uicontrol('Parent',hData2PanelSub,'Style','text','String','','Position',[60,5,50,690]);
htext4Y = uicontrol('Parent',hData2PanelSub,'Style','text','String','','Position',[115,5,50,690]);

htext5 = uicontrol('Style','text','String','Jump:','Position',[10,820,50,30],'FontSize',10);
htext6 = uicontrol('Style','text','String','','Position',[70,820,60,30],'FontSize',10);
htext7 = uicontrol('Style','text','String','Frame:','Position',[290,820,50,30],'FontSize',10);
htext8 = uicontrol('Style','text','String','','Position',[350,820,60,30],'FontSize',10);
htext9 = uicontrol('Style','text','String','LLR:','Position',[560,820,50,30],'FontSize',10);
htext10 = uicontrol('Style','text','String','','Position',[620,820,100,30],'FontSize',10);
htext11 = uicontrol('Style','text','String','PR:','Position',[730,820,50,30],'FontSize',10);
htext12 = uicontrol('Style','text','String','','Position',[770,820,100,30],'FontSize',10);
htext13 = uicontrol('Style','text','String','Accept:','Position',[150,820,60,30],'FontSize',10);
htext14 = uicontrol('Style','text','String','','Position',[210,820,50,30],'FontSize',10);
htext20 = uicontrol('Style','text','String','NCurrent:','Position',[420,820,60,30],'FontSize',10);
htext21 = uicontrol('Style','text','String','','Position',[490,820,50,30],'FontSize',10);

hedit1 = uicontrol('Style','edit','Position',[930,820,50,25],'CallBack',@getEdit1);
htext15 = uicontrol('Style','text','String','Frame','Position',[930,850,50,20],'CallBack',@getEdit1);
hpush1 = uicontrol('Style','pushbutton','String','go Frame','Position',[930,780,50,30],'CallBack',@getPush1);
hpush2 = uicontrol('Style','pushbutton','String','Previous','Position',[870,780,50,30],'CallBack',@getPush2);
hpush3 = uicontrol('Style','pushbutton','String','Next','Position',[990,780,50,30],'CallBack',@getPush3);

hedit2 = uicontrol('Style','edit','Position',[1060,820,50,25],'CallBack',@getEdit2);
set(hedit2,'String',1);
htext15 = uicontrol('Style','text','String','Start Frame','Position',[1050,850,70,20]);
hpush4 = uicontrol('Style','pushbutton','String','run','Position',[1060,780,50,30],'CallBack',@getPush4);

hBgPanel = uipanel('Parent',guiFig,'Title','BG','Position',[0.905 0.875 0.085 0.12]);
htext16 = uicontrol('Parent',hBgPanel,'Style','text','String','Accepted BG','Position',[5,70,70,20]);
htext17 = uicontrol('Parent',hBgPanel,'Style','text','String','','Position',[5,48,70,20]);
htext18 = uicontrol('Parent',hBgPanel,'Style','text','String','Proposed BG','Position',[5,26,70,20]);
htext19 = uicontrol('Parent',hBgPanel,'Style','text','String','','Position',[5,4,70,20]);

hPanelIm1 = uipanel('Parent',guiFig,'Title','Proposed','Position',[0.011 0.011 0.315 0.43]);
hAxesImProposed = axes('Parent',hPanelIm1,'Units','pixels','Position',[10,10,385,385],'Visible','off');
hPanelIm2 = uipanel('Parent',guiFig,'Title','Accepted','Position',[0.335 0.011 0.315 0.43]);
hAxesImAccepted = axes('Parent',hPanelIm2,'Units','pixels','Position',[10,10,385,385],'Visible','off');
hPanelIm3 = uipanel('Parent',guiFig,'Title','Residuum','Position',[0.66 0.011 0.315 0.43]);
hAxesImRes = axes('Parent',hPanelIm3,'Units','pixels','Position',[10,10,385,385],'Visible','off');

plot(hAxis2,x,y,'ok')
plot(hAxis1,x,y,'ok')
xlim(hAxis2,[0,BoxSize]);ylim(hAxis2,[0,BoxSize]);
    xlabel(hAxis2,'X(pixel)');ylabel(hAxis2,'Y(pixel)');
    set(hAxis2,'Ydir','reverse','XAxisLocation','top')
 xlim(hAxis1,[0,BoxSize]);ylim(hAxis1,[0,BoxSize]);
    xlabel(hAxis1,'X(pixel)');ylabel(hAxis1,'Y(pixel)');
    set(hAxis1,'Ydir','reverse','XAxisLocation','top')    
DataT = obj.Data(:,:,ROINum);
imshow(obj.Data(:,:,ROINum),[min(DataT(:)),max(DataT(:))],'Parent',hAxesImProposed)
imshow(obj.Data(:,:,ROINum),[min(DataT(:)),max(DataT(:))],'Parent',hAxesImAccepted)
imshow(obj.Data(:,:,ROINum),[min(DataT(:)),max(DataT(:))],'Parent',hAxesImRes)

function getEdit1(~,~)
    StartState = str2double(get(hedit1,'String'));
    L = StartState;
end
function getEdit2(~,~)
    StartState = str2double(get(hedit2,'String'));
    L = length(ChainROI)-1;
end
function getPush1(~,~)
    Loop = 1;
    getEdit1();
    runLoop()    
end
function getPush2(~,~)
    Loop = 1;
    getEdit1();
    StartState = StartState - 1;
    L = L-1;
    runLoop();
    set(hedit1,'String',StartState);
end
function getPush3(~,~)
    Loop = 1;
    getEdit1();
    StartState = StartState + 1;
    L = L+1;
    runLoop();   
    set(hedit1,'String',StartState);
end
function getPush4(~,~)
    Loop = 1;
    if strcmp(get(hpush4,'String'),'run')
        set(hpush4,'String','stop')
    else
        set(hpush4,'String','run')
        Loop = 0;
        return;
    end
    getEdit2();
    runLoop();
    
end

function slider1(src,eventdata,arg1)
    val = get(src,'Value');
    set(arg1,'Position',[0 -val 1 2])
end

function slider2(src,eventdata,arg1)
    val = get(src,'Value');
    set(arg1,'Position',[0 -val 1 2])    
end

function runLoop()
[Xg,Yg]=meshgrid((0.5:BoxSize-0.5),(0.5:BoxSize-0.5));
PSF = obj.RJStruct.PSF_Sigma;
for nn = StartState:L
    X = ChainROI(nn).X;
    Y = ChainROI(nn).Y;
    I = ChainROI(nn).Photons;
    Im = 0;
    for ii = 1:length(X)
        Im = Im + I(ii)*normpdf(Xg,X(ii),PSF).*normpdf(Yg,Y(ii),PSF); 
    end
    Im = Im + ChainROI(nn).BG;
    ImFused = imfuse(obj.Data(:,:,ROINum),Im,'Scaling','joint','ColorChannels',[1 2 0]);
    imshow(ImFused,'Parent',hAxesImAccepted);
    if X < -5
        continue;
    end
    Xp = PChainROI(nn).X;
    Yp = PChainROI(nn).Y;
    Ip = PChainROI(nn).Photons;
    Imp = 0;
    for ii = 1:length(Xp)
        Imp = Imp + Ip(ii)*normpdf(Xg,Xp(ii),PSF).*normpdf(Yg,Yp(ii),PSF); 
    end
    Imp = Imp + PChainROI(nn).BG;
    ImpFused = imfuse(obj.Data(:,:,ROINum),Imp,'Scaling','joint','ColorChannels',[1 2 0]);
    imshow(ImpFused,'Parent',hAxesImProposed);
    ResIm = obj.Data(:,:,ROINum)-Imp;
    imshow(ResIm,[min(ResIm(:)),max(ResIm(:))],'Parent',hAxesImRes);
    Type = ChainROI(nn).JumpType;
    switch Type
        case 1
            Stype =sprintf('Jump1');
        case 2
            Stype =sprintf('JumpBG');
        case 3
            Stype =sprintf('Jump2');
        case 4
            Stype =sprintf('Split');
        case 5
            Stype =sprintf('Merge');
        case 6
            Stype =sprintf('Birth');
        case 7
            Stype =sprintf('Death');
        case 8
            Stype =sprintf('Cons');
        case 9
            Stype =sprintf('Dest');
        case 10
            Stype =sprintf('Jump3');
    end
    plot(hAxis2,X,Y,'.')
    hold(hAxis2,'on');
    plot(hAxis2,x,y,'ok')
    xlim(hAxis2,[0,BoxSize]);ylim(hAxis2,[0,BoxSize]);
    xlabel(hAxis2,'X(pixel)');ylabel(hAxis2,'Y(pixel)');
    set(hAxis2,'Ydir','reverse','XAxisLocation','top')
    plot(hAxis1,Xp,Yp,'.r')
    hold(hAxis1,'on');
    plot(hAxis1,x,y,'ok')
    xlim(hAxis1,[0,BoxSize]);ylim(hAxis1,[0,BoxSize]);
    xlabel(hAxis1,'X(pixel)');ylabel(hAxis1,'Y(pixel)');
    set(hAxis1,'Ydir','reverse','XAxisLocation','top')
    
    l = length(I);
    lp = length(Ip);
    CP = cell(1,l);
    CX = cell(1,l);
    CY = cell(1,l);
    CPP = cell(1,lp);
    CXP = cell(1,lp);
    CYP = cell(1,lp);
    for mm = 1:l
        CP{1,mm} = uint32(I(mm));
        CX{1,mm} = X(mm);
        CY{1,mm} = Y(mm);
    end 
    for mm = 1:lp
        CPP{1,mm} = uint32(Ip(mm));
        CXP{1,mm} = Xp(mm);
        CYP{1,mm} = Yp(mm); 
    end
    set(htext4P,'String',sprintf('%d\n',CP{:}));
    set(htext4X,'String',sprintf('%.03f\n',CX{:}));
    set(htext4Y,'String',sprintf('%.03f\n',CY{:}));
    set(htext3P,'String',sprintf('%d\n',CPP{:}));
    set(htext3X,'String',sprintf('%.03f\n',CXP{:}));
    set(htext3Y,'String',sprintf('%.03f\n',CYP{:}));
    set(htext6,'String',Stype);
    set(htext8,'String',sprintf('%d\n',nn));
    set(htext10,'String',sprintf('%g\n',ChainROI(nn).LLR));
    set(htext12,'String',sprintf('%g\n',ChainROI(nn).PR));
    set(htext14,'String',sprintf('%g\n',ChainROI(nn).Accepted));
    set(htext17,'String',ChainROI(nn).BG);
    set(htext19,'String',PChainROI(nn).BG);
    set(htext21,'String',ChainROI(nn).N);
    pause(0.05);
    hold(hAxis2,'off');
    hold(hAxis1,'off');
    if Loop == 0
        set(hedit1,'String',nn); 
        return; 
    end
end
set(hedit1,'String',nn);
end
end