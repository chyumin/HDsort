classdef MaxWellFile < hdsort.file.SingleFileWrapper
    % Wrapper for BEL-MEA recording file format.
    % This object does not support multi-file handling!
    % Its main purpose is to open recording files and to run preprocessing
    % such that the data is in a format that can be used directly for
    % spike-sorting (see hdsort.file.CMOSMEA).
    
    properties (Hidden)
        fileID
    
        nDims
        dims
        maxDims
        dimAddChans
        dataSets
    end
    
    properties
        P
        connectedChannel
        chip_version
        chip_number
        software_version
    end
    
    
    methods
        
        %------------------------------------------------------------------
        %% CONSTRUCTOR
        function self = MaxWellFile(fileName, varargin)
            
            P.mapping_ds_name = '/mapping';
            P.signal_ds_name = '/sig';
            P.dateTimeFormate = 'yyyy-MM-dd''T''HH:mm:ss';
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            assert(~iscell(fileName), 'Only one file can be used as input!')
            samplesPerSecond = 20000;
            
%             %% Extract Multielectrode:
%             mapping = h5read(fileName, P.mapping_ds_name);
%             
%             x = h5info(fileName, P.signal_ds_name);
%             nChannels = x.Dataspace.Size(2)-1;
%             allRecordedChannels = (1:nChannels)';
%             mappedChannels = double(mapping.channel)+1;
%             allUnsortedChannels = [mappedChannels; allRecordedChannels(~ismember(allRecordedChannels, mappedChannels))];
%             [allSortedChannels sortidx] = sort(allUnsortedChannels);
%             
%             xPos = [double(mapping.x); zeros(nChannels - size(mapping.x, 1), 1)-1 ];
%             yPos = [double(mapping.y); zeros(nChannels - size(mapping.y, 1), 1)-1 ];
%             elns = [double(mapping.electrode); zeros(nChannels - size(mapping.electrode, 1), 1)-1 ];
%             
%             ME = hdsort.file.MultiElectrode([xPos(sortidx) yPos(sortidx)], elns(sortidx));
            
            
            %% Extract Multielectrode:
            mapping = h5read(fileName, P.mapping_ds_name);
            
            % This block is for perfomance reasons:
            fileID = H5F.open(fileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            datasetID = H5D.open(fileID, P.signal_ds_name);
            dataspaceID = H5D.get_space(datasetID);
            [nDims, h5_dims, h5_maxDims] = H5S.get_simple_extent_dims(dataspaceID);
            nChannels = h5_dims(1) - 1; % ignore channel 1025, beause it contains only the command counter
            H5D.close(datasetID);
            H5F.close(fileID);
            
            allRecordedChannels = (1:nChannels)';
            mappedChannels = double(mapping.channel)+1;
            allUnsortedChannels = [mappedChannels; allRecordedChannels(~ismember(allRecordedChannels, mappedChannels))];
            [allSortedChannels, sortidx] = sort(allUnsortedChannels);
            
            xPos = [double(mapping.x); zeros(nChannels - size(mapping.x, 1), 1)-1 ];
            yPos = [double(mapping.y); zeros(nChannels - size(mapping.y, 1), 1)-1 ];
            elns = [double(mapping.electrode); zeros(nChannels - size(mapping.electrode, 1), 1)-1 ];
            
            ME = hdsort.file.MultiElectrode([xPos(sortidx) yPos(sortidx)], elns(sortidx));
            
            
            
            
            self = self@hdsort.file.SingleFileWrapper('MaxWellFile', samplesPerSecond, ME, fileName, varargin{:})
            
            %% Test the multielectrode:
            for ii = 1:size(self.MultiElectrode.electrodeNumbers, 1)
                nEl = self.MultiElectrode.electrodeNumbers(ii);
                
                if nEl < 0
                    continue
                end
                
                nCh = ii;
                idx = find(double(mapping.channel)+1 == nCh );
                
                if ~isempty(mapping.electrode(idx))
                    assert( mapping.electrode(idx) == nEl, 'Oops!');
                else
                    assert( ~any(mapping.electrode == nEl) , 'Oops!');
                end
            end
            
            %%
            self.P = P;
            self.dataSets.mapping.name = P.mapping_ds_name;
            self.dataSets.signal.name = P.signal_ds_name;
            self.dataSets.frame_numbers.name = '/ephys/frame_numbers';
            self.dataSets.bits.name = '/bits';
            self.dataSets.raw_spikes.name = '/proc0/spikeTimes';
            self.dataSets.gain.name = '/settings/gain';
            self.dataSets.version.name = '/version';
            self.dataSets.raw_spikes.name = '/proc0/spikeTimes';
            
            self.fileName = fileName;
            self.connectedChannel = find(self.MultiElectrode.electrodeNumbers > -1);
            self = self.openH5File();
            
            self.chip_version = 'MaxOne';
            self.software_version = self.tryReadingDataset(self.dataSets.version.name, false);
            
            self.info = ['This object is intended to wrap a single MaxWell-MEA recording file.'];
        end
        
        %------------------------------------------------------------------
        function self = openH5File(self)
            plist = 'H5P_DEFAULT';
            rmode = 'H5F_ACC_RDONLY';
            
            try
                self.fileID = H5F.open(self.fileName, rmode, plist);
            catch
                str = hdsort.util.buildLastErrString();
                disp(str);
                error('Could not open file %s!', self.fileName);
            end
            
            try
                self.dataSets.signal.datasetID = H5D.open(self.fileID, self.dataSets.signal.name);
                self.dataSets.signal.dataspaceID = H5D.get_space(self.dataSets.signal.datasetID);
            catch
                str = hdsort.util.buildLastErrString();
                disp(str);
                error('Could not find H5 Variable %s!', self.dataSets.signal.name);
            end
            self.getProps();
            
            %% Frame Numbers:
%             try
%                 self.dataSets.frame_numbers.datasetID = H5D.open(self.fileID, self.dataSets.frame_numbers.name);
%                 self.dataSets.frame_numbers.dataspaceID = H5D.get_space(self.dataSets.frame_numbers.datasetID);
%             catch
%                 str = hdsort.util.buildLastErrString();
%                 disp(str);
%                 error('Could not find H5 Variable %s!', self.dataSets.frame_numbers.name);
%             end
%             
%             %% Other info:
%             self.chip_version    = self.tryReadingDataset(self.dataSets.chip_version.name);
%             self.chip_number     = self.tryReadingDataset(self.dataSets.chip_number.name, false);
%             self.software_version = self.tryReadingDataset(self.dataSets.software_version.name);
%             
%             self.loadBits();
%             self.readHardwareSettings();
        end
        
        %------------------------------------------------------------------
        function h5info_ = getH5Info(self)
            if ~isfield(self.buffer, 'h5info') || isempty(self.buffer.h5info)
                self.buffer.h5info = h5info(self.fileName);
            end
            h5info_ = self.buffer.h5info;
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
        function getProps(self)
            [self.nDims, h5_dims, h5_maxDims] = H5S.get_simple_extent_dims(self.dataSets.signal.dataspaceID);
            self.dims = fliplr(h5_dims);
            self.maxDims = fliplr(h5_maxDims);
        end
        
        %------------------------------------------------------------------
        function [missingFrames] = getMissingFrameNumbers(self)
            if ~isfield(self.buffer, 'missingFrames') || isempty(self.buffer.missingFrames)
                self.buffer.missingFrames = hdsort.file.util.getMissingFrameNumbers(self.getFrameNumbers());
            end
            missingFrames = self.buffer.missingFrames;
        end
        
        %------------------------------------------------------------------
        function gain = getGain(self)
            if ~isfield(self.buffer, 'gain') || isempty(self.buffer.gain)
                self.buffer.gain = self.tryReadingDataset(self.dataSets.gain.name, false);
            end
            gain = self.buffer.gain;
        end
        
        %------------------------------------------------------------------
        function LSB_volts = getLSB_(self)
            gain = self.getGain();
            LSB_volts = 3.3 / (1024.0 * gain);
        end
        
        %------------------------------------------------------------------
        function rawSpikes = getRawSpikes(self)
            if ~isfield(self.buffer, 'rawSpikes') || isempty(self.buffer.rawSpikes)
                rawSpikes_ = self.tryReadingDataset(self.dataSets.raw_spikes.name, false);
                if ~isempty(rawSpikes_)            
                    self.buffer.rawSpikes.spiketimes = double(rawSpikes_.frameno);
                    self.buffer.rawSpikes.channel = double(rawSpikes_.channel) + 1;
                    self.buffer.rawSpikes.positions = self.MultiElectrode.electrodePositions(self.buffer.rawSpikes.channel, :);
                    self.buffer.rawSpikes.electrode = self.MultiElectrode.electrodeNumbers(self.buffer.rawSpikes.channel);
                    self.buffer.rawSpikes.amplitude = double(rawSpikes_.amplitude);
                else
                    self.buffer.rawSpikes = [];
                end
            end
            rawSpikes = self.buffer.rawSpikes;
        end
        
        %------------------------------------------------------------------
        function self = closeH5File(self)
            try
                H5D.close(self.datasetID);
            catch
                warning('could not close h5 dataset!');
            end
            try
                H5F.close(self.fileID);
            catch
                warning('could not close h5 file!');
            end
        end
        
        %------------------------------------------------------------------
        function disp(self)
            details(self)
            disp(['Channels: [' num2str(self.dims(1)) ' ' num2str(size(self.MultiElectrode.electrodePositions,1)) ']']);
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, varargin)
            sz = self.size;
            [bb relIdx] = hdsort.file.hdf5.getBoundingBoxFromIndexing(sz, varargin{:});
            if nargin == 2
                % if only one index was requested, only ask for one
                % dimension
                bb = bb(:,1);
            end
            %   assert(~any(bb(:)<=0) && ~any(any(bb>repmat(sz(1:size(bb,2)),size(bb,1),1))), 'Indexing out of bounds!');
            block = bb(2,:) - bb(1,:) +1;
            numEl = prod(block);
            assert(numEl < 10^9, sprintf('Loading that much data (%d = %s) at once is not recommended!', numEl, sprintf('%d x ', block)));
            offset = bb(1,:) - 1;
            % read outer bounding box of requested data, this is usually faster
            % than reading individual rows or indices
            X = double(hdsort.file.hdf5.read_dset(self.dataSets.signal.datasetID, block, offset));
            % select actually requested data
            if ~isempty(relIdx)
                X = X(relIdx{:});
            end
        end
        
        %------------------------------------------------------------------
        function frameNo=getFrameNoAt(self, index)
            X = h5read( self.fileName , self.dataSets.signal.name, [index 1027], [1, 2]);
            frameNo = bitor( bitshift( double(X(:,2)) , 16 ) , double(X(:,1)) );
        end
        
        %------------------------------------------------------------------
        function [frameNo] = getFrameNumbers(self, varargin)
            if ~isfield(self.buffer, 'frame_numbers') || isempty(self.buffer.frame_numbers)
                X = h5read( self.fileName , self.dataSets.signal.name, [1 1027], [self.size(1), 2]);
                self.buffer.frame_numbers = bitor( bitshift( double(X(:,2)) , 16 ) , double(X(:,1)) )';
            end
            frameNo = self.buffer.frame_numbers;
        end
        
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = self.dims(1);
        end
        
        %------------------------------------------------------------------
        function fileName = getFullFileName(self)
            fileName = self.fileName;
        end
        
        %------------------------------------------------------------------
        function rootRileName = getRootFileName(self)
            [path_ name_ ext_] = fileparts(self.fileName);
            [path_ rootRileName ext_] = fileparts(name_);
        end
        
        %------------------------------------------------------------------
        function restrictToConnectedChannels(self)
            self.restrictToChannels();
            self.restrictToChannels(self.connectedChannel);
        end
        
    end
end