classdef BELMEAFile < hdsort.file.SingleFileWrapper
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
        function self = BELMEAFile(fileName, varargin)
            P.mapping_ds_name = '/ephys/mapping';
            P.sampling_rate_ds_name = '/ephys/frame_rate';
            P.signal_ds_name = '/ephys/signal';
            
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            assert(~iscell(fileName), 'Only one file can be used as input!')
            samplesPerSecond = double(h5read(fileName, P.sampling_rate_ds_name));
            
            
            %% Extract Multielectrode:
            mapping = h5read(fileName, P.mapping_ds_name);
            
            % This block is for perfomance reasons:
            fileID = H5F.open(fileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            datasetID = H5D.open(fileID, P.signal_ds_name);
            dataspaceID = H5D.get_space(datasetID);
            [~, h5_dims, ~] = H5S.get_simple_extent_dims(dataspaceID);
            nChannels = h5_dims(1); % ignore channel 1025, beause it contains only the command counter
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
            
            self = self@hdsort.file.SingleFileWrapper('BELMEAFile', samplesPerSecond, ME, fileName, varargin{:})
            
            %% Test the multielectrode:
            for ii = 1:size(self.MultiElectrode.electrodeNumbers, 1)
                nEl = self.MultiElectrode.electrodeNumbers(ii);
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
            
            self.dataSets.sampling_rate.name = P.sampling_rate_ds_name;
            self.dataSets.mapping.name = P.mapping_ds_name;
            self.dataSets.signal.name = P.signal_ds_name;
            self.dataSets.frame_numbers.name = '/ephys/frame_numbers';
            self.dataSets.additional_channels.name = '/ephys/additional_channels';
            
            self.dataSets.chip_version.name = '/chipinformation/chip_version';
            self.dataSets.chip_number.name = '/chipinformation/chip_number';
            self.dataSets.software_version.name = '/chipinformation/software_version';
            self.dataSets.hardware_settings.name = '/hardware_settings';
            self.dataSets.raw_spikes.name = '/spikesorting/spikeTimes';
            
            self.dataSets.status.name = '/status';
            self.dataSets.bits.name = '/bits';
            
            self.fileName = fileName;
            self.connectedChannel = find(self.MultiElectrode.electrodeNumbers > -1);
            
            self = self.openH5File();
            
            self.info = ['This object is intended to wrap a single BEL-MEA recording file.'];
        end
        
        %------------------------------------------------------------------
        function self = openH5File(self)%, bReadOnly)
            plist = 'H5P_DEFAULT';
            %if nargin == 1 || ~bReadOnly
            %    rmode = 'H5F_ACC_RDWR';
            %else
            rmode = 'H5F_ACC_RDONLY';
            %end
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
            try
                self.dataSets.frame_numbers.datasetID = H5D.open(self.fileID, self.dataSets.frame_numbers.name);
                self.dataSets.frame_numbers.dataspaceID = H5D.get_space(self.dataSets.frame_numbers.datasetID);
            catch
                str = hdsort.util.buildLastErrString();
                disp(str);
                error('Could not find H5 Variable %s!', self.dataSets.frame_numbers.name);
            end
            
            %% Additional Channels:
            try
                self.dataSets.additional_channels.datasetID = H5D.open(self.fileID, self.dataSets.additional_channels.name);
                self.dataSets.additional_channels.dataspaceID = H5D.get_space(self.dataSets.additional_channels.datasetID);
                
                [nDims h5_dims h5_maxDims] = H5S.get_simple_extent_dims(self.dataSets.additional_channels.dataspaceID);
                self.dimAddChans = fliplr(h5_dims);
            catch
                %str = hdsort.util.buildLastErrString();
                %disp(str);
                self.dimAddChans = [0,0];
            end
            
            %% Other info:
            self.chip_version = self.tryReadingDataset(self.dataSets.chip_version.name);
            self.chip_number  = self.tryReadingDataset(self.dataSets.chip_number.name, false);
            self.software_version = self.tryReadingDataset(self.dataSets.software_version.name);
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
            
            if iscell(out) && numel(out) == 1
                out = out{1};
            end
        end
        
        %------------------------------------------------------------------
        function x = getH5Info(self)
            if isempty(self.buffer.h5info)
                self.buffer.h5info = h5info(self.fileName);
            end
            x = self.buffer.h5info;
        end
        
        %------------------------------------------------------------------
        function getProps(self)
            [self.nDims, h5_dims, h5_maxDims] = H5S.get_simple_extent_dims(self.dataSets.signal.dataspaceID);
            self.dims = fliplr(h5_dims);
            self.maxDims = fliplr(h5_maxDims);
        end
        
        %------------------------------------------------------------------
        function bits = getBits(self)
            if ~isfield(self.buffer, 'bits') || isempty(self.buffer.bits)
                bits = self.tryReadingDataset(self.dataSets.bits.name, false);
                if ~isempty(bits) 
                self.buffer.bits = [bits.x0x2Fephys0x2Fframe_numbers, bits.x0x2Fbits];
                end
            end
            bits = self.buffer.bits;
        end
        
        %------------------------------------------------------------------
        function status = getStatus(self)
            if ~isfield(self.buffer, 'status') || isempty(self.buffer.status)
                status_ = self.tryReadingDataset(self.dataSets.status.name, false);
                if ~isempty(status_)
                    self.buffer.status = double([status_.x0x2Fephys0x2Fframe_numbers, status_.x0x2Fstatus]);
                end
            end
            status = self.buffer.status;
        end
        
        %------------------------------------------------------------------
        function hws = getHardwareSettings(self)
            if ~isfield(self.buffer, 'hardware_settings') || isempty(self.buffer.hardware_settings)
                try
                    hwinfo = h5info(self.fileName, self.dataSets.hardware_settings.name);
                    hws_names = {hwinfo.Datasets.Name};
                    
                    for h_ = hws_names
                        hw_name = [self.dataSets.hardware_settings.name '/' h_{1}];
                        self.buffer.hardware_settings.(h_{1}) = self.tryReadingDataset(hw_name);
                    end
                catch
                    warning('No hardware info found!')
                    self.buffer.hardware_settings = [];
                end
            end
            hws = self.buffer.hardware_settings;
        end
        
        %------------------------------------------------------------------
        function gain = getGain_(self)
            if strcmp(self.chip_version, 'Mea1k')
                hws = self.getHardwareSettings();
                gain = hdsort.file.BELMEAFile.getGainMea1k(hws);
            else
                warning(['Unknown chip version: ' self.chip_version])
                gain = 1.0;
            end
        end
        
        %------------------------------------------------------------------
        function LSB_volts = getLSB_(self)
            if ~isfield(self.buffer, 'lsb') || isempty(self.buffer.lsb)
                if strcmp(self.chip_version, 'Mea1k')
                    hws = self.getHardwareSettings();
                    self.buffer.lsb = hdsort.file.BELMEAFile.getLSBMea1k(hws);
                else
                    warning(['Unknown chip version: ' self.chip_version])
                    self.buffer.lsb = 1;
                end
            end
            LSB_volts = self.buffer.lsb;
        end
        
        %------------------------------------------------------------------
        function rawSpikes = getRawSpikes(self)
            if ~isfield(self.buffer, 'rawSpikes') || isempty(self.buffer.rawSpikes)
                rawSpikes_ = self.tryReadingDataset(self.dataSets.raw_spikes.name, false);
                if ~isempty(rawSpikes_)
                    self.buffer.rawSpikes.spiketimes = double(rawSpikes_.x0x2Fephys0x2Fframe_numbers);
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
            disp(['Restriction: [' num2str(self.dims(1)) ' ' num2str(size(self.MultiElectrode.electrodePositions,1)) ']']);
            disp(['Add. chans.: [' num2str(self.dimAddChans(1)) ' ' num2str(self.dimAddChans(2)) ']']);
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
            assert(numEl < 10^10, sprintf('Loading that much data (%d = %s) at once is not recommended!', numEl, sprintf('%d x ', block)));
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
        function X = getAdditionalChannels(self, varargin)
            if nargin == 2
                error('Please provide frame index and channel number!')
            end
            
            sz = self.dimAddChans;
            [bb relIdx] = hdsort.file.hdf5.getBoundingBoxFromIndexing(sz, varargin{:});
            
            block = bb(2,:) - bb(1,:) +1;
            numEl = prod(block);
            assert(numEl < 10^9, sprintf('Loading that much data (%d = %s) at once is not recommended!', numEl, sprintf('%d x ', block)));
            offset = bb(1,:) - 1;
            % read outer bounding box of requested data, this is usually faster
            % than reading individual rows or indices
            X = double(hdsort.file.hdf5.read_dset(self.dataSets.additional_channels.datasetID, block, offset));
            % select actually requested data
            if ~isempty(relIdx)
                X = X(relIdx{:});
            end
        end
        
        %------------------------------------------------------------------
        function [frameNo] = getFrameNumbers(self, varargin)
            sz = self.size; sz(2) = 1;
            
            [bb relIdx] = hdsort.file.hdf5.getBoundingBoxFromIndexing(sz, varargin{:});
            assert(~any(bb(:)<=0) && ~any(any(bb>repmat(sz(1:size(bb,2)),size(bb,1),1))), 'Indexing out of bounds!');
            block = bb(2,1) - bb(1,1) +1;
            offset = bb(1,1) - 1;
            
            try
            frameNo = double(hdsort.file.hdf5.read_dset(self.dataSets.frame_numbers.datasetID, block, offset));
            catch
            frameNo = double(hdsort.file.hdf5.read_dset(self.dataSets.frame_numbers.datasetID, [block, 1], [offset 0]));
            end
            frameNo = frameNo(:)';
            
            % select actually requested data
            if ~isempty(relIdx)
                frameNo = frameNo(relIdx{:});
            end
        end
        
        %------------------------------------------------------------------
        function [missingFrames] = getMissingFrameNumbers(self)
            if ~isfield(self.buffer, 'missingFrames') || isempty(self.buffer.missingFrames)
                self.buffer.missingFrames = hdsort.file.util.getMissingFrameNumbers(self.getFrameNumbers());
            end
            missingFrames = self.buffer.missingFrames;
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
    
    methods (Static)
        
        %------------------------------------------------------------------
        function gain = getGainMea1k(hws)
            gain_ = stage1(hws.Stage1Gain, hws.Stage1Bypass) * stage2(hws.Stage2Gain, hws.Stage2Bypass) * stage3(hws.Stage3Gain, hws.Stage3Bypass);
            gain = double(gain_);
            function ret = stage1(s, sbypass)
                if s < 0 || sbypass
                    ret = 1;
                    return
                end
                ret = (s*16 + (1-s)*7);
            end
            
            function ret = stage2(s, sbypass)
                if s < 0 || sbypass
                    ret = 1;
                    return
                end
                ret = 2^(s-1);
            end
            
            function ret = stage3(s, sbypass)
                if s < 0 || sbypass
                    ret = 1;
                    return
                end
                ret = 2^(s+1);
            end
        end
        
        %------------------------------------------------------------------
        function LSB_volts = getLSBMea1k(hws)
            gain = hdsort.file.BELMEAFile.getGainMea1k(hws);
            LSB_volts = 3.3 / (1024.0 * gain);
        end
        
    end
    
end