function [Data] = loadData(FileName,varargin)
%loadData Loads raw microsope dataset from mat, ics, h5 or file
% INPUT
%    FileName - full path to datafile, must be mat, ics or h5
%    varargin - input parameters, different for each file type:
%       mat: MatVarName, name of matlab variable containing the data
%       ics: no parameters needed
%       h5: DatasetIdx, index of dataset to be loaded from h5 file
%           ChannelIdx (optional), index of channel to be loaded, default is 1
% OUTPUT
%    Data - data loaded from FileName, converted to type single
% REQUIRES
%    Dipimage (www.diplib.org)

% Marjolein Meddens, Lidke Lab, 2017
% Sandeep Pallikkuth, Lidke Lab, 2018

FileType = FileName(end-2:end);

switch FileType
    case 'mat'
        % check input
        if isempty(varargin)
            error('SMA_Core:loadData:noMatVarName','No Matlab variable name given for mat file, please give as input: Data = SMA_Core.loadData(FileName,MatVarName)')
        end
        MatVarName = varargin{1};
        if iscell(MatVarName)
            MatVarName=MatVarName{1};
        end
        % load data
        tmp = load(FileName);
        Data = single(tmp.(MatVarName));
    case 'ics'
        % load data
        tmp = readim(FileName);
        Data = single(tmp);
    case '.h5'
        if isempty(varargin)
            error('SMA_Core:loadData:noDatasetIndex','No dataset index given for h5 file, please give as input: Data = SMA_Core.loadData(FileName,DatasetIndex)')
        elseif nargin == 2
            ChannelIdx = 1;
        elseif nargin > 2
            ChannelIdx = varargin{2};
        end
        DatasetIdx = varargin{1};
        HD5Info = h5info(FileName);
       % setup directory into H5 file
        for ii = 1 : numel(HD5Info.Groups)
            if strcmp(HD5Info.Groups(ii).Name,'/Data')
                ChannelName = sprintf('Channel%02i',ChannelIdx);
                DataSetName = sprintf('Data%04i',DatasetIdx);
                DataName = sprintf('/Data/%s/%s',ChannelName,DataSetName);
            elseif strcmp(HD5Info.Groups(ii).Name,'/Channel01')
                ChannelName = sprintf('Zposition%03i',ChannelIdx);
                DataSetName = sprintf('Data%04i',DatasetIdx);
                DataName = sprintf('/Channel01/%s/%s',ChannelName,DataSetName);
            end
        end
        % check whether channel and dataset exist
        for ii = 1 : numel(HD5Info.Groups)
            if strcmp(HD5Info.Groups(ii).Name,'/Data')
                DataGroup = HD5Info.Groups(ii);
                break
            elseif strcmp(HD5Info.Groups(ii).Name,'/Channel01')
                DataGroup = HD5Info.Groups(ii);
                break
            end
        end
        % check channel input
        ChannelExists = 0;
        for ii = 1 : numel(DataGroup.Groups)
            if strcmp(DataGroup.Groups(ii).Name,['/Data/' ChannelName])
                ChannelGroup = DataGroup.Groups(ii);
                ChannelExists = 1;
                break
            elseif strcmp(DataGroup.Groups(ii).Name,['/Channel01/' ChannelName])
                ChannelGroup = DataGroup.Groups(ii);
                ChannelExists = 1;
                break
           end
        end
        if ~ChannelExists
            error('SMA_Core:loadData:unknownChannel','h5 file does not contain %s, cannot load data from that channel',ChannelName);
        end
        % check dataset input
        DataSetExists = 0;
        for ii = 1 : numel(ChannelGroup.Datasets)
            if strcmp(ChannelGroup.Datasets(ii).Name,DataSetName)
                DataSetExists = 1;
                break
            end
        end
        if ~DataSetExists
            error('SMA_Core:loadData:unknownDataset','h5 file does not contain dataset %s in %s, cannot load data',DataSetName,ChannelName);
        end
        % load data
        Data=single(h5read(FileName,DataName));
    otherwise
        % unknown file type
        error('SMA_Core:loadData:unknownFileType','Unknown type of file, datafile should be of type mat, ics or h5')
end

end



