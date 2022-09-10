function gui(obj)
%Graphical User Interface for BAMF. 
%   The gui displays eight panels. The DATA-panel has buttons to load the 
%   raw data and drifts. The RJStruct-panel tabulates default values for 
%   the RJMCMC parameters, which can be altered by the user. The P 
%   Burnin-panel and P Trial-panel let the user change the probabilities of 
%   proposing different jumps in the burn-in part of the chain and the trial 
%   portion code are included in the Parameters-panel. The Threshold-panel 
%   and FrameConnect-panel contain the post-processing parameters. The 
%   post-processing can be disabled by unchecking the boxes next to their 
%   panels.
%
% INPUTS:
%   obj:   The class object.
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
guiFig = figure('Visible','on','Position',[250,320,1010,500],'MenuBar','none','Name','RJMCMC user interface','NumberTitle','off','Interruptible','off');

hDataPanel = uipanel('Parent',guiFig,'Title','DATA','Position',[0.025 0.65 0.95 0.33]);

fileButton1 = uicontrol('Parent',hDataPanel, 'Style', 'pushbutton', 'String','Choose File(s)','Enable','on','Position', [5 95 120 40],'Callback',@getFileList);
fileList1 = uicontrol('Parent',hDataPanel,'Style','listbox','String',' ','Position',[135 7 170 128],'BackgroundColor',[1 1 1]);
deleteButton1 = uicontrol('Parent',hDataPanel, 'Style', 'pushbutton', 'String','Delete Selected File','Enable','on','Position', [5 45 120 40],'Callback',@deleteFiles);

fileButton2 = uicontrol('Parent',hDataPanel, 'Style', 'pushbutton', 'String','Choose Drift-X','Enable','on','Position', [320 95 120 40],'Callback',@getDriftX);
fileList2 = uicontrol('Parent',hDataPanel,'Style','listbox','String',' ','Position',[450 7 170 128],'BackgroundColor',[1 1 1]);
deleteButton2 = uicontrol('Parent',hDataPanel, 'Style', 'pushbutton', 'String','Delete Selected File','Enable','on','Position', [320 45 120 40],'Callback',@deleteDriftX);

fileButton3 = uicontrol('Parent' ,hDataPanel, 'Style', 'pushbutton', 'String','Choose Drift-Y','Enable','on','Position', [635 95 120 40],'Callback',@getDriftY);
fileList3 = uicontrol('Parent',hDataPanel,'Style','listbox','String',' ','Position',[765 7 170 128],'BackgroundColor',[1 1 1]);
deleteButton3 = uicontrol('Parent',hDataPanel, 'Style', 'pushbutton', 'String','Delete Selected File','Enable','on','Position', [635 45 120 40],'Callback',@deleteDriftY);

hRJparametersPanel = uipanel('Parent',guiFig,'Title','RJStruct','Position',[0.025 0.285 0.3 0.35]);

NameRJ1 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','Zoom','Position',[10,130,60,20]);
ValueRJ1 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[10,110,60,20],'CallBack',@getEditRJ1);
set(ValueRJ1,'String','20');
NameRJ2 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','PSFsigma','Position',[10,80,60,20]);
ValueRJ2 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditRJ2);
set(ValueRJ2,'String','1.3');
NameRJ3 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','NTrials','Position',[10,30,60,20]);
ValueRJ3 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[10,10,60,20],'CallBack',@getEditRJ3);
set(ValueRJ3,'String','2000');
NameRJ4 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','NBurnin','Position',[85,130,60,20]);
ValueRJ4 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[85,110,60,20],'CallBack',@getEditRJ4);
set(ValueRJ4,'String','3000');
NameRJ5 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','NTrialMCMC','Position',[85,80,60,20]);
ValueRJ5 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[85,60,60,20],'CallBack',@getEditRJ5);
set(ValueRJ5,'String','1000');
NameRJ6 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','NBurninMC','Position',[85,30,60,20]);
ValueRJ6 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[85,10,60,20],'CallBack',@getEditRJ6);
set(ValueRJ6,'String','1000');
NameRJ7 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','P_dP','Position',[160,130,60,20]);
ValueRJ7 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[160,110,60,20],'CallBack',@getEditRJ7);
set(ValueRJ7,'String','10');
NameRJ8 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','Rho','Position',[160,80,60,20]);
ValueRJ8 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[160,60,60,20],'CallBack',@getEditRJ8);
set(ValueRJ8,'String','0.01');
NameRJ9 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','Xstd','Position',[160,30,60,20]);
ValueRJ9 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[160,10,60,20],'CallBack',@getEditRJ9);
set(ValueRJ9,'String','0.1');
NameRJ10 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','Istd','Position',[235,130,60,20]);
ValueRJ10 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[235,110,60,20],'CallBack',@getEditRJ10);
set(ValueRJ10,'String','10');
NameRJ11 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','BGstd','Position',[235,80,60,20]);
ValueRJ11 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[235,60,60,20],'CallBack',@getEditRJ11);
set(ValueRJ11,'String','1');
NameRJ12 = uicontrol('Parent',hRJparametersPanel,'Style','text','String','BgPriorRatio','Position',[235,30,65,20]);
ValueRJ12 = uicontrol('Parent',hRJparametersPanel,'Style','edit','Position',[235,10,60,20],'CallBack',@getEditRJ12);
set(ValueRJ12,'String','12');

hPBurninPanel = uipanel('Parent',guiFig,'Title','P_Burnin','Position',[0.445 0.185 0.16 0.45]);

NamePBin1 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_Jump','Position',[10,180,60,20]);
ValuePBin1 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[10,160,60,20],'CallBack',@getEditPBin1);
set(ValuePBin1,'String','0.3')
NamePBin2 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_Split','Position',[10,130,60,20]);
ValuePBin2 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[10,110,60,20],'CallBack',@getEditPBin2);
set(ValuePBin2,'String','0.1');
NamePBin3 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_Merge','Position',[10,80,60,20]);
ValuePBin3 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditPBin3);
set(ValuePBin3,'String','0.1');
NamePBin4 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_Birth','Position',[85,180,60,20]);
ValuePBin4 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[85,160,60,20],'CallBack',@getEditPBin4);
set(ValuePBin4,'String','0.1');
NamePBin5 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_Death','Position',[85,130,60,20]);
ValuePBin5 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[85,110,60,20],'CallBack',@getEditPBin5);
set(ValuePBin5,'String','0.1');
NamePBin6 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_Convert','Position',[10,30,60,20]);
ValuePBin6 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[10,10,60,20],'CallBack',@getEditPBin6);
set(ValuePBin6,'String','0.1');
NamePBin7 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_GSplit','Position',[85,80,60,20]);
ValuePBin7 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[85,60,60,20],'CallBack',@getEditPBin7);
set(ValuePBin7,'String','0.1');
NamePBin8 = uicontrol('Parent',hPBurninPanel,'Style','text','String','P_GMerge','Position',[85,30,60,20]);
ValuePBin8 = uicontrol('Parent',hPBurninPanel,'Style','edit','Position',[85,10,60,20],'CallBack',@getEditPBin8);
set(ValuePBin8,'String','0.1');

hPTrialPanel = uipanel('Parent',guiFig,'Title','P_Trial','Position',[0.622 0.185 0.16 0.45]);

NamePTin1 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_Jump','Position',[10,180,60,20]);
ValuePTin1 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[10,160,60,20],'CallBack',@getEditPTin1);
set(ValuePTin1,'String','0.4')
NamePTin2 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_Split','Position',[10,130,60,20]);
ValuePTin2 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[10,110,60,20],'CallBack',@getEditPTin2);
set(ValuePTin2,'String','0');
NamePTin3 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_Merge','Position',[10,80,60,20]);
ValuePTin3 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditPTin3);
set(ValuePTin3,'String','0');
NamePTin4 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_Birth','Position',[85,180,60,20]);
ValuePTin4 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[85,160,60,20],'CallBack',@getEditPTin4);
set(ValuePTin4,'String','0.05');
NamePTin5 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_Death','Position',[85,130,60,20]);
ValuePTin5 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[85,110,60,20],'CallBack',@getEditPTin5);
set(ValuePTin5,'String','0.05');
NamePTin6 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_Convert','Position',[10,30,60,20]);
ValuePTin6 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[10,10,60,20],'CallBack',@getEditPTin6);
set(ValuePTin6,'String','0.2');
NamePTin7 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_GSplit','Position',[85,80,60,20]);
ValuePTin7 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[85,60,60,20],'CallBack',@getEditPTin7);
set(ValuePTin7,'String','0.15');
NamePTin8 = uicontrol('Parent',hPTrialPanel,'Style','text','String','P_GMerge','Position',[85,30,60,20]);
ValuePTin8 = uicontrol('Parent',hPTrialPanel,'Style','edit','Position',[85,10,60,20],'CallBack',@getEditPTin8);
set(ValuePTin8,'String','0.15');

hSMAparametersPanel = uipanel('Parent',guiFig,'Title','SMAStruct','Position',[0.343 0.285 0.085 0.35]);

NameSMA1 = uicontrol('Parent',hSMAparametersPanel,'Style','text','String','MeanPhotons','Position',[10,130,70,20]);
ValueSMA1 = uicontrol('Parent',hSMAparametersPanel,'Style','edit','Position',[10,110,60,20],'CallBack',@getEditSMA1);
set(ValueSMA1,'String','800')
NameSMA2 = uicontrol('Parent',hSMAparametersPanel,'Style','text','String','MinPhotons','Position',[10,80,60,20]);
ValueSMA2 = uicontrol('Parent',hSMAparametersPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditSMA2);
set(ValueSMA2,'String','200');
NameSMA3 = uicontrol('Parent',hSMAparametersPanel,'Style','text','String','MinPvalue','Position',[10,30,60,20]);
ValueSMA3 = uicontrol('Parent',hSMAparametersPanel,'Style','edit','Position',[10,10,60,20],'CallBack',@getEditSMA3);
set(ValueSMA3,'String','0.01');

hParametersPanel = uipanel('Parent',guiFig,'Title','Parameters','Position',[0.8 0.185 0.177 0.45]);

NameP1 = uicontrol('Parent',hParametersPanel,'Style','text','String','Gain','Position',[10,180,60,20]);
ValueP1 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[10,160,60,20],'CallBack',@getEditP1);
set(ValueP1,'String','11')
NameP2 = uicontrol('Parent',hParametersPanel,'Style','text','String','Offset','Position',[10,130,60,20]);
ValueP2 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[10,110,60,20],'CallBack',@getEditP2);
set(ValueP2,'String','100');
NameP3 = uicontrol('Parent',hParametersPanel,'Style','text','String','BoxSize','Position',[10,80,60,20]);
ValueP3 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditP3);
set(ValueP3,'String','16');
NameP4 = uicontrol('Parent',hParametersPanel,'Style','text','String','Overlap','Position',[85,180,60,20]);
ValueP4 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[85,160,60,20],'CallBack',@getEditP4);
set(ValueP4,'String','3');
NameP5 = uicontrol('Parent',hParametersPanel,'Style','text','String','KeepChain(N/Y)','Position',[85,130,80,20]);
ValueP5 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[85,110,60,20],'CallBack',@getEditP5);
set(ValueP5,'String','N');
NameP6 = uicontrol('Parent',hParametersPanel,'Style','text','String','UsePPT(N/Y)','Position',[85,80,70,20]);
ValueP6 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[85,60,60,20],'CallBack',@getEditP6);
set(ValueP6,'String','Y');
NameP7 = uicontrol('Parent',hParametersPanel,'Style','text','String','FindPSF(N/Y)','Position',[85,30,70,20]);
ValueP7 = uicontrol('Parent',hParametersPanel,'Style','edit','Position',[85,10,60,20],'CallBack',@getEditP7);
set(ValueP7,'String','N');

hThreshPanel = uipanel('Parent',guiFig,'Title','Threshold','Position',[0.05 0.02 0.16 0.25]);
hThreshCheck = uicontrol('Style','checkbox','Value',1,'Position',[27 103 20 40],'Callback',@enableThresh);

NameThresh1 = uicontrol('Parent',hThreshPanel,'Style','text','String','Max XY-SE','Position',[10,80,60,20]);
ValueThresh1 = uicontrol('Parent',hThreshPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditThresh1);
set(ValueThresh1,'String','0.25')
NameThresh2 = uicontrol('Parent',hThreshPanel,'Style','text','String','Min XY-SE','Position',[10,30,60,20]);
ValueThresh2 = uicontrol('Parent',hThreshPanel,'Style','edit','Position',[10,10,60,20],'CallBack',@getEditThresh2);
set(ValueThresh2,'String','0.001');
NameThresh3 = uicontrol('Parent',hThreshPanel,'Style','text','String','MaxPhoton_SE','Position',[82,80,75,20]);
ValueThresh3 = uicontrol('Parent',hThreshPanel,'Style','edit','Position',[85,60,60,20],'CallBack',@getEditThresh3);
set(ValueThresh3,'String','50');
NameThresh4 = uicontrol('Parent',hThreshPanel,'Style','text','String','MinPhotons','Position',[85,30,60,20]);
ValueThresh4 = uicontrol('Parent',hThreshPanel,'Style','edit','Position',[85,10,60,20],'CallBack',@getEditThresh4);
set(ValueThresh4,'String','200');

hFConnectPanel = uipanel('Parent',guiFig,'Title','FrameConnect','Position',[0.25 0.02 0.16 0.25]);
hFCCheck = uicontrol('Style','checkbox','Value',0,'Position',[230 103 20 40],'Callback',@enableFC);

NameFC1 = uicontrol('Parent',hFConnectPanel,'Style','text','String','MinPValue','Position',[10,80,60,20]);
ValueFC1 = uicontrol('Parent',hFConnectPanel,'Style','edit','Position',[10,60,60,20],'CallBack',@getEditFC1);
set(ValueFC1,'String','0.01')
NameFC2 = uicontrol('Parent',hFConnectPanel,'Style','text','String','MaxFrame-gap','Position',[10,30,75,20]);
ValueFC2 = uicontrol('Parent',hFConnectPanel,'Style','edit','Position',[10,10,60,20],'CallBack',@getEditFC2);
set(ValueFC2,'String','10');
NameFC3 = uicontrol('Parent',hFConnectPanel,'Style','text','String','MaxDistance','Position',[85,80,70,20]);
ValueFC3 = uicontrol('Parent',hFConnectPanel,'Style','edit','Position',[85,60,60,20],'CallBack',@getEditFC3);
set(ValueFC3,'String','4');
NameFC4 = uicontrol('Parent',hFConnectPanel,'Style','text','String','LOS','Position',[85,30,60,20]);
ValueFC4 = uicontrol('Parent',hFConnectPanel,'Style','edit','Position',[85,10,60,20],'CallBack',@getEditFC4);
set(ValueFC4,'String','0.01');

hDriftCheck = uicontrol('Style','checkbox','Value',0,'Position',[450 30 20 40],'Callback',@enableDrift);
hDriftStatic = uicontrol('Style','text','string','Drift','Position',[465 17 30 40])

runButton1 = uicontrol('Style', 'pushbutton', 'String','run RJMCMC','Enable','on','Position', [850 30 120 40],'Callback',@runRJMCMC);

    function getFileList(~,~) % Choose data file(s) and list them in the File List box
        if isempty(obj.DataDir)
            [filename, pathname]=uigetfile('Y:\*.mat;*.ics;*.h5','MultiSelect','on');
        else
            % open the saved directory path
            dirpath=strcat(obj.DataDir,'*.mat;*.ics;*.h5');
            [filename, pathname]=uigetfile(dirpath,'MultiSelect','on');
        end
        obj.DataDir=pathname;
        if ~isa(filename,'cell')
            filename = {filename};
        end
        obj.FileList=filename';
        obj.SaveDir=fullfile(obj.DataDir,'RJ_Results');
        if ~isdir(obj.SaveDir)
            mkdir(obj.SaveDir)
        end
        set(fileList1,'String',filename)
    end

    function deleteFiles(~,~)
        currentItems=get(fileList1,'String');
        rowToDelete = get(fileList1, 'Value');
        newItems = currentItems;
        newItems(rowToDelete) = [];
        set(fileList1, 'String', newItems);  % display new filelist
        if ~isempty(newItems) 
            if length(newItems) >= 1 
                set(fileList1, 'Value', 1); 
            end 
        end
        obj.FileList=newItems;
    end

    function getDriftX(~,~)
        if isempty(obj.DataDir)
            [Filename, pathname]=uigetfile('Y:\*.mat;*.ics;*.h5','MultiSelect','on');
        else
            % open the saved directory path
            dirpath=strcat(obj.DataDir,'*.mat;*.ics;*.h5');
            [Filename, pathname]=uigetfile(dirpath,'MultiSelect','on');
        end
        if ~isa(Filename,'cell')
            Filename = {Filename};
        end
        obj.DriftXDir=pathname;
        obj.FileListDriftX=Filename;
        obj.DriftFlag = 0;
        set(hDriftCheck,'Value',0);
        set(fileList2,'String',Filename)
    end

    function deleteDriftX(~,~)
        currentItems=get(fileList2,'String');
        rowToDelete = get(fileList2, 'Value');
        newItems = currentItems;
        newItems(rowToDelete) = [];
        set(fileList2, 'String', newItems);  % display new filelist
        if ~isempty(newItems) 
            if length(newItems) >= 1 
                set(fileList2, 'Value', 1); 
            end 
        end
        obj.RJStruct.XDrift = [];
    end

    function getDriftY(~,~)
        if isempty(obj.DataDir)
            [Filename, pathname]=uigetfile('Y:\*.mat;*.ics;*.h5','MultiSelect','on');
        else
            % open the saved directory path
            dirpath=strcat(obj.DataDir,'*.mat;*.ics;*.h5');
            [Filename, pathname]=uigetfile(dirpath,'MultiSelect','on');
        end
        if ~isa(Filename,'cell')
            Filename = {Filename};
        end
        obj.DriftYDir=pathname;
        obj.FileListDriftY=Filename;
        obj.DriftFlag = 0;
        set(hDriftCheck,'Value',0);
        set(fileList3,'String',Filename)
    end

    function deleteDriftY(~,~)
        currentItems=get(fileList3,'String');
        rowToDelete = get(fileList3, 'Value');
        newItems = currentItems;
        newItems(rowToDelete) = [];
        set(fileList3, 'String', newItems);  % display new filelist
        if ~isempty(newItems) 
            if length(newItems) >= 1 
                set(fileList3, 'Value', 1); 
            end 
        end
        obj.RJStruct.YDrift = [];
    end

    function getEditRJ1(~,~) 
        obj.RJStruct.SRZoom = str2double(get(ValueRJ1,'String'));
    end

    function getEditRJ2(~,~) 
        obj.RJStruct.PSF_Sigma = str2double(get(ValueRJ2,'String'));
    end

    function getEditRJ3(~,~)
        obj.RJStruct.N_Trials = str2double(get(ValueRJ3,'String'));
    end

    function getEditRJ4(~,~)
        obj.RJStruct.N_Burnin = str2double(get(ValueRJ4,'String'));
    end

    function getEditRJ5(~,~)
        obj.RJStruct.N_TrialsMCMC = str2double(get(ValueRJ5,'String'));
    end

    function getEditRJ6(~,~)
        obj.RJStruct.N_BurninMCMC = str2double(get(ValueRJ6,'String'));
    end

    function getEditRJ7(~,~)
        obj.RJStruct.P_dP = str2double(get(ValueRJ7,'String'));
    end

    function getEditRJ8(~,~)
        obj.RJStruct.Rho = str2double(get(ValueRJ8,'String'));
    end

    function getEditRJ9(~,~)
        obj.RJStruct.XstdFg = str2double(get(ValueRJ9,'String'));
    end

    function getEditRJ10(~,~)
        obj.RJStruct.IstdFg = str2double(get(ValueRJ10,'String'));
    end

    function getEditRJ11(~,~)
        obj.RJStruct.BGstd = str2double(get(ValueRJ11,'String'));
    end

    function getEditRJ12(~,~)
        obj.RJStruct.BgPriorRange = str2double(get(ValueRJ12,'String'));
    end

    function getEditPBin1(~,~)
        obj.RJStruct.P_Burnin(1) = str2double(get(ValuePBin1,'String'));
    end

    function getEditPBin2(~,~)
        obj.RJStruct.P_Burnin(2) = str2double(get(ValuePBin2,'String'));
    end

    function getEditPBin3(~,~)
        obj.RJStruct.P_Burnin(3) = str2double(get(ValuePBin3,'String'));
    end

    function getEditPBin4(~,~)
        obj.RJStruct.P_Burnin(4) = str2double(get(ValuePBin4,'String'));
    end

    function getEditPBin5(~,~)
        obj.RJStruct.P_Burnin(5) = str2double(get(ValuePBin5,'String'));
    end

    function getEditPBin6(~,~)
        obj.RJStruct.P_Burnin(8) = str2double(get(ValuePBin6,'String'));
    end
    
    function getEditPBin7(~,~)
        obj.RJStruct.P_Burnin(6) = str2double(get(ValuePBin7,'String'));
    end

    function getEditPBin8(~,~)
        obj.RJStruct.P_Burnin(7) = str2double(get(ValuePBin8,'String'));
    end

    function getEditPTin1(~,~)
        obj.RJStruct.P_Trial(1) = str2double(get(ValuePTin1,'String'));
    end

    function getEditPTin2(~,~)
        obj.RJStruct.P_Trial(2) = str2double(get(ValuePTin2,'String'));
    end

    function getEditPTin3(~,~)
        obj.RJStruct.P_Trial(3) = str2double(get(ValuePTin3,'String'));
    end

    function getEditPTin4(~,~)
        obj.RJStruct.P_Trial(4) = str2double(get(ValuePTin4,'String'));
    end

    function getEditPTin5(~,~)
        obj.RJStruct.P_Trial(5) = str2double(get(ValuePTin5,'String'));
    end

    function getEditPTin6(~,~)
        obj.RJStruct.P_Trial(8) = str2double(get(ValuePTin6,'String'));
    end

    function getEditPTin7(~,~)
        obj.RJStruct.P_Trial(6) = str2double(get(ValuePTin7,'String'));
    end

    function getEditPTin8(~,~)
        obj.RJStruct.P_Trial(7) = str2double(get(ValuePTin8,'String'));
    end

    function getEditSMA1(~,~)
        obj.SMAStruct.MeanPhotons = str2double(get(ValueSMA1,'String'));
    end

    function getEditSMA2(~,~)
        obj.SMAStruct.MinPhotons = str2double(get(ValueSMA2,'String'));
    end

    function getEditSMA3(~,~)
        obj.SMAStruct.MinPValue = str2double(get(ValueSMA3,'String'));
    end

    function getEditP1(~,~)
        obj.Camera_Gain = str2double(get(ValueP1,'String'));
    end

    function getEditP2(~,~)
        obj.Camera_Offset = str2double(get(ValueP2,'String'));
    end

    function getEditP3(~,~)
        obj.BoxSize = str2double(get(ValueP3,'String'));
    end

    function getEditP4(~,~)
        obj.Overlap = str2double(get(ValueP4,'String'));
    end

    function getEditP5(~,~)
        S = get(ValueP5,'String');
        if strcmp(S,'Y')
            obj.KeepChain = 1;
        elseif strcmp(S,'N')
            obj.KeepChain = 0;
        else
           error('The string could be either N or Y.'); 
        end
    end

    function getEditP6(~,~)
        S = get(ValueP6,'String');
        if strcmp(S,'Y')
            obj.UsePPToolbox = 1;
        elseif strcmp(S,'N')
            obj.UsePPToolbox = 0;
        else
           error('The string could be either N or Y.'); 
        end
    end

    function getEditP7(~,~)
        S = get(ValueP7,'String');
        if strcmp(S,'Y')
            obj.FindPSF = 1;
        elseif strcmp(S,'N')
            obj.FindPSF = 0;
        else
           error('The string could be either N or Y.'); 
        end
    end
    function enableThresh(~,~)
        obj.ThreshFlag = get(hThreshCheck,'Value');
    end

    function getEditThresh1(~,~)
        obj.Threshold.Max_XY_SE = str2double(get(ValueThresh1,'String'));
    end

    function getEditThresh2(~,~)
        obj.Threshold.Min_XY_SE = str2double(get(ValueThresh2,'String'));
    end

    function getEditThresh3(~,~)
        obj.Threshold.Max_Photons_SE = str2double(get(ValueThresh3,'String'));
    end

    function getEditThresh4(~,~)
        obj.Threshold.MinPhotons = str2double(get(ValueThresh4,'String'));
    end

    function enableFC(~,~)
        obj.FCFlag = get(hFCCheck,'Value');
    end

    function getEditFC1(~,~)
        obj.FrameConnect.MinPValue = str2double(get(ValueFC1,'String'));
    end

    function getEditFC2(~,~)
        obj.FrameConnect.MaxFrame_gap = str2double(get(ValueFC2,'String'));
    end

    function getEditFC3(~,~)
        obj.FrameConnect.MaxDistance = str2double(get(ValueFC3,'String'));
    end

    function getEditFC4(~,~)
        obj.FrameConnect.LOS = str2double(get(ValueFC4,'String'));
    end

    function enableDrift(~,~)
        obj.DriftFlag = get(hThreshCheck,'Value');
    end

    function runRJMCMC(~,~)
        for TimeID = 1:length(obj.FileList)
            obj.processData(TimeID);
        end
        obj=obj.assembleResults();
        if obj.ThreshFlag
            obj = obj.thresholdRJ();
        end
        if obj.DriftFlag
            SMDin.X = obj.TotClust.X;
            SMDin.Y = obj.TotClust.Y;
            SMDin.Photons = obj.TotClust.Photons;
            SMDin.X_SE = obj.TotClust.X_SE;
            SMDin.Y_SE = obj.TotClust.Y_SE;
            SMDin.FrameNum = obj.TotClust.FrameNum;
            SMDin.Bg = zeros(size(SMDin.FrameNum));
            SMDin.DatasetNum = obj.TotClust.DatasetNum;
            SMDin.XSize = size(obj.Data,2); SMDin.YSize = size(obj.Data,1);      
            SMDin.Nframes=max(SMDin.FrameNum);
            SMDin.Ndatasets=max(SMDin.DatasetNum);
            OptParams.PDegree      = 1;
            OptParams.SRZoom       = obj.RJStruct.SRZoom;
            OptParams.TolFun_intra = 1e-2;
            OptParams.TolX_intra   = 1e-4;
            OptParams.TolFun_inter = 1e-2;
            OptParams.TolX_inter   = 1e-4;
            OptParams.Init_inter   = 1;
            [obj.SMD, ~] = SMA_Core.driftCorrectKNN(SMDin, OptParams);
        end
        if obj.FCFlag
            if isempty(obj.SMD)
                SMDin.X = obj.TotClust.X;
                SMDin.Y = obj.TotClust.Y;
                SMDin.Photons = obj.TotClust.Photons;
                SMDin.X_SE = obj.TotClust.X_SE;
                SMDin.Y_SE = obj.TotClust.Y_SE;
                SMDin.FrameNum = obj.TotClust.FrameNum;
                SMDin.Bg = zeros(size(SMDin.FrameNum));
                SMDin.DatasetNum = obj.TotClust.DatasetNum;
                SMDin.XSize = size(obj.Data,2); 
                SMDin.YSize = size(obj.Data,1);      
                SMDin.Nframes=max(SMDin.FrameNum);
                SMDin.Ndatasets=max(SMDin.DatasetNum);
            end
             SMDin.ThreshFlag = zeros(size(obj.SMD.X));
             SMDin.Photons_SE = obj.TotClust.Photons_SE;
             SMDin.Bg_SE = 0.1*ones(size(obj.SMD.X));
             SMDin.LogL = 0.1*ones(size(obj.SMD.X));
             SMDin.PSFSigma = ones(size(obj.SMD.X));
             SMDin.PSFSigma_SE = ones(size(obj.SMD.X));
             LOS=0.01;
             MaxDistance=2;
             MaxFrameGap=5;
             FitType='GaussianBasic';
             [obj.SMD,obj.SMD_combined]=SMA_Core.frameConnect(SMDin,LOS,MaxDistance,MaxFrameGap,FitType);
        end
        Im=0;
        if obj.FCFlag
            SMR = obj.SMD_combined;
            Struct = 'SMD_combined';
        elseif obj.FCFlag==0 && obj.DriftFlag
            SMR = obj.SMD;
            Struct = 'SMD';
        elseif obj.FCFlag==0 && obj.DriftFlag==0 && obj.ThreshFlag 
            SMR = obj.TotClust;
            Struct = 'TotStruct';
        end
        for nn = 1:length(obj.FileList)
            Ind = SMR.DatasetNum == nn;
            ClustInfoT.X = SMR.X(Ind);
            ClustInfoT.Y = SMR.Y(Ind);
            ClustInfoT.I = SMR.Photons(Ind);
            ClustInfoT.X_SE = SMR.X_SE(Ind);
            ClustInfoT.Y_SE = SMR.Y_SE(Ind);
            ClustInfoT.I_SE = ones(size(ClustInfoT.X));
            ClustInfoT.FrameNum = SMR.FrameNum(Ind);
            ClustInfoT.ROInum = ones(size(ClustInfoT.X));
            Im = Im+obj.makeGaussIm(ClustInfoT);
        end
        RJ.produceIm(Im,fullfile(obj.SaveDir,'MCMC'));
        save(fullfile(obj.SaveDir,'MCMC'),'SMR','Im','-v7.3');
        PlotFlag = 0;
        SaveFlag = 1;
        obj.makeHist(Struct,'X_SE',PlotFlag,SaveFlag);
        obj.makeHist(Struct,'Y_SE',PlotFlag,SaveFlag);
        obj.makeHist(Struct,'Photons',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','Photons_SE',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','BG',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','ABG',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','BBG',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','JumpRJ',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','SplitRJ',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','MergeRJ',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','BirthRJ',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','DeathRJ',PlotFlag,SaveFlag);
        obj.makeHist('TotClust','JumpMC',PlotFlag,SaveFlag);
    end
    
end
