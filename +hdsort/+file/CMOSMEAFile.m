classdef CMOSMEAFile < hdsort.file.SingleFileWrapper
    
    properties
        sourceFname
        dataSets
        filter_settings
        
        h5matrix_raw
        session_str
        size_buffer
        CL
        connected_channels
        
        lsb
    end
    
    methods
        
        %------------------------------------------------------------------
        %% CONSTRUCTOR
        function self = CMOSMEAFile(fileName, varargin)
            if nargin == 1
                warning('Not recommended to use this object directly. Better use multi-file supporting hdsort.file.CMOSMEA!')
            end
            
            session_str_ = '/Sessions/Session0/'; % This wrapper only supports one session per file!
            samplesPerSecond = hdf5read(fileName, [session_str_ 'sr']);   
            samplesPerSecond = samplesPerSecond(1);
            
            % build the multielectrode
            % CHECK IF THIS IS A NEW FILE VERSION WITH THE CHANNEL LIST
            % SPLIT UP
            try 
                % THIS IS MUCH FASTER
                connectedChannels = hdf5read(fileName, [session_str_ 'channel_connected']);
                nC_ = length(connectedChannels);
                connectedChannels = connectedChannels == 1;
                assert(any(connectedChannels), 'None of the channels was connected!')
                cpx = double(hdf5read(fileName, [session_str_ 'channel_posx']));
                cpy = double(hdf5read(fileName, [session_str_ 'channel_posy']));
                cnr = double(hdf5read(fileName, [session_str_ 'channel_nr']));
                cpx = cpx(connectedChannels==1);
                cpx = cpx(:);
                cpy = cpy(connectedChannels==1);
                cpy = cpy(:);
                ME = hdsort.file.MultiElectrode([cpx cpy], cnr);
            catch
                % Seems to be the old format
                % WHY EVER, BUT THIS READ CALL IS SLOW LIKE HELL
                CL__ = hdf5read(fileName, [session_str_ 'channel_list']);
                CL__ = get(CL__, 'Data');
                nC = size(CL__,1);
                nV = length(CL__{1});
                CL_ = zeros(nC, nV);
                for i=1:nC
                    CL_(i,:) = cellfun(@double, CL__{i});
                end
                connectedChannels = CL_(:,2)==1;
                CL_(~connectedChannels, :) = [];   % remove unconnected channels
                CL_(:,3:4) = CL_(:,3:4)/1000; % convert to micro meter
                ME = hdsort.file.MultiElectrode(CL_(:,3:4), CL_(:,5));
            end
            
            self = self@hdsort.file.SingleFileWrapper('CMOSMEAFile', samplesPerSecond, ME, fileName, varargin{:})
            
            self.fileName = fileName;
            self.session_str = session_str_;
            self.size_buffer = [];
            
            %find out whether sig is data of a string if a binary file name:
            fileInfo = h5info(self.fileName, [self.session_str 'sig']);
            is_int = strcmp(fileInfo.Datatype.Class, 'H5T_INTEGER');
            is_str = strcmp(fileInfo.Datatype.Class, 'H5T_STRING') ;
            assert(is_int || is_str, 'sig has to be in the format of an integer or a string');
            
            if is_str
                binFileName = cell2mat( h5read(self.fileName, [self.session_str 'sig']) );
                
                [pathstr,name,ext] = fileparts(self.fileName);
                binFile = fullfile(pathstr, binFileName);
                
                binDims = h5read(self.fileName, [self.session_str 'bin_dims']);
                assert( exist(binFile, 'file') == 2, ['Task aborted: binary file ' binFile ' not found!']);
                
                self.h5matrix_raw =  hdsort.file.util.BinaryFileMatrix(binFile, binDims);
            else
                self.h5matrix_raw = hdsort.file.hdf5.matrix(self.fileName, [self.session_str 'sig'], true);
            end
            
            try
                self.sourceFname = get(hdf5read(self.fileName, [self.session_str 'fileName']), 'data');
            catch
                self.sourceFname = [];
            end
            
            self.MultiElectrode.setDataSource(self);
            
            %% Get LSB:
            self.lsb = self.getLSB_();
            self.connected_channels = find(connectedChannels);
            
            %% Read filter settings:
            self.dataSets.filter_settings.name = [self.session_str 'filter'];
            self.readFilterSettings();
            
            self.info = ['This object is intended to wrap a single preprocessed CMOSMEA datafile. ' ...
                         'It can be used as a standalone object, but it is recommended to use it ' ...
                         'within hdsort.file.CMOSMEA, an object that provides the same functionality ' ...
                         'but for multiple files at the same time.'];
        end
        
        %------------------------------------------------------------------
        function [out] = tryReadingDataset(self, dsName, displayError)
            if nargin < 3
                displayError = true;
            end
            try
                out = h5read(self.fileName, dsName);
            catch
                if displayError
                    str = hdsort.util.buildLastErrString();
                    disp(str);
                end
                out = [];
            end
        end
        
        %------------------------------------------------------------------
        function readFilterSettings(self)
            try
                fsinfo = h5info(self.fileName, self.dataSets.filter_settings.name);
                fs_names = {fsinfo.Datasets.Name};
                
                for f_ = fs_names
                    fs_name = [self.dataSets.filter_settings.name '/' f_{1}];
                    self.filter_settings.(f_{1}) = self.tryReadingDataset(fs_name);
                end
            catch
                warning('Not filter information found!')
            end
        end
        
        %------------------------------------------------------------------
        function missingFrames = getMissingFrameNumbers(self)
            FR = hdsort.file.hdf5.recursiveLoad(self.fileName, [self.session_str 'frame_numbers']);
            if ~isfield(FR, 'missing_fns')
                warning('There is a problem with missing framenumbers in CMOSMEASession.getFrameNumbers()!')
                FR.missing_fns = -1;
            end
            
            missingFrames.first = FR.first_fn;
            missingFrames.last = FR.last_fn;
            if FR.missing_fns == -1
                missingFrames.begin = [];
                missingFrames.n = 0;
                missingFrames.length = [];
            else
                missingFrames.begin = FR.missing_fns(1, :);
                missingFrames.n = size(FR.missing_fns, 2);
                missingFrames.length = FR.missing_fns(2, :);
            end
        end
        
        %------------------------------------------------------------------
        function x = isBinaryFile(self)
            x = isa( self.h5matrix_raw, 'hdsort.file.util.BinaryFileMatrix');
        end
        
        %------------------------------------------------------------------
        function gain = getGain(self)
            gain = hdf5read(self.fileName, [self.session_str 'gain']);
            assert(length(gain) == 4, 'Gain must have 4 values!');
            try
                gainmultiplier = hdf5read(self.fileName, [self.session_str 'filter/gainmultiplier']);
            catch
                gainmultiplier = 1;
            end
            gain(1) = gain(1)*single(gainmultiplier);
        end
        
        %------------------------------------------------------------------
        function lsb = getLSB_(self)
            try
                ntk.adc_resolution = hdf5read(self.fileName, [self.session_str 'adc_resolution']);
            catch
                ntk.adc_resolution = 8;
            end
            try
                ntk.adc_range = hdf5read(self.fileName, [self.session_str 'adc_range']);
            catch
                ntk.adc_range = 2.9;
            end            
            gain = self.getGain();
            if ~isempty(gain) && any(gain~=1)
                lsb = ntk.adc_range/(2^ntk.adc_resolution-1)/gain(1)*1000000;
            else
                lsb = 1;
            end
        end
        
        %------------------------------------------------------------------
        function CL = getChannelList(self)
            if isempty(self.CL)
                CL_ = hdf5read(self.fileName, [self.session_str 'channel_list']);
                CL_ = get(CL_, 'Data');
                nC = size(CL_,1);
                nV = length(CL_{1});
                self.CL = zeros(nC, nV);
                for i=1:nC
                    self.CL(i,:) = cellfun(@double, CL_{i});
                end
            end
            CL = self.CL;
        end    
        
        %------------------------------------------------------------------
        function wf = getWaveform_(self, nCut, channelIndex, varargin)
            channelIndex = self.connected_channels(channelIndex);
            wf = self.lsb * double( self.h5matrix_raw.getWaveform_(nCut, channelIndex, varargin{:}) );
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, timeIndex, channelIndex, channels)
            if nargin < 4
                channels = [];
            end
            
            if isempty(channels) || strcmp(channels, 'connected')
                assert(max(channelIndex)<=length(self.connected_channels), 'channel index out of bounds')
                channelIndex = self.connected_channels(channelIndex);
            end
            
            % channel idx is first dimension !!!
            X = double( self.h5matrix_raw(timeIndex, channelIndex) ) * self.lsb;
        end
        
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = size(self.h5matrix_raw,1);
        end
        
        %------------------------------------------------------------------
        function gainMultiplier = getGainMultiplier(self)
            try
                gainMultiplier = double(h5read(self.fileName, '/Sessions/Session0/filter/gainmultiplier'));
            catch
                gainMultiplier = 1.0;
                warning('GainMultiplier not defined!')
            end
        end
        
    end
end
