classdef BELMEAFile < hdsort.filewrapper.SingleFileWrapper
    % Wrapper for BEL-MEA recording file format.
    % This object does not support multi-file handling!
    % Its main purpose is to open recording files and to run preprocessing
    % such that the data is in a format that can be used directly for
    % spike-sorting (see hdsort.filewrapper.CMOSMEA).
    
    properties
        P
        h5info
        
        fileID
        %className
        
        
        nDims
        dims
        maxDims
        connectedChannel
        
        bits
        missingFrames
        
        chip_version
        chip_number
        software_version
        
        hardware_settings
        
        dimAddChans
        
        dataSets
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
            
            x = h5info(fileName, P.signal_ds_name);
            nChannels = x.Dataspace.Size(2)-1;
            allRecordedChannels = (1:nChannels)';
            mappedChannels = double(mapping.channel)+1;
            allUnsortedChannels = [mappedChannels; allRecordedChannels(~ismember(allRecordedChannels, mappedChannels))];
            [allSortedChannels sortidx] = sort(allUnsortedChannels);
            
            xPos = [double(mapping.x); zeros(nChannels - size(mapping.x, 1), 1)-1 ];
            yPos = [double(mapping.y); zeros(nChannels - size(mapping.y, 1), 1)-1 ];
            elns = [double(mapping.electrode); zeros(nChannels - size(mapping.electrode, 1), 1)-1 ];
            
            ME = hdsort.filewrapper.MultiElectrode([xPos(sortidx) yPos(sortidx)], elns(sortidx));
            
            self = self@hdsort.filewrapper.SingleFileWrapper('BELMEAFile', samplesPerSecond, ME, fileName, varargin{:})
            
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
            
            self.dataSets.bits.name = '/bits';
            
            self.fileName = fileName;
            %self.className = P.className;
            self.h5info = h5info(fileName);
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
            self.chip_version    = self.tryReadingDataset(self.dataSets.chip_version.name);
            self.chip_number     = self.tryReadingDataset(self.dataSets.chip_number.name, false);
            self.software_version = self.tryReadingDataset(self.dataSets.software_version.name);
            
            self.loadBits();
            
            self.readHardwareSettings();
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
        function loadBits(self)
            bits = self.tryReadingDataset(self.dataSets.bits.name, false);
            if isempty(bits) return; end
            self.bits = [bits.x0x2Fephys0x2Fframe_numbers, bits.x0x2Fbits];
        end
        
        %------------------------------------------------------------------
        function hs = readHardwareSettings(self)
            try
                hwinfo = h5info(self.fileName, self.dataSets.hardware_settings.name);
                hws_names = {hwinfo.Datasets.Name};
            
                for h_ = hws_names
                    hw_name = [self.dataSets.hardware_settings.name '/' h_{1}];
                    self.hardware_settings.(h_{1}) = self.tryReadingDataset(hw_name);
                end
            catch
                disp('No hardware info found!')
                self.hardware_settings = [];
            end
            hs = self.hardware_settings;
        end
        
        %------------------------------------------------------------------
        function getProps(self)
            [self.nDims h5_dims h5_maxDims] = H5S.get_simple_extent_dims(self.dataSets.signal.dataspaceID);
            self.dims = fliplr(h5_dims);
            self.maxDims = fliplr(h5_maxDims);
        end
        
        %------------------------------------------------------------------
        function bits = getBits(self)
            bits = self.bits;
        end
        
        %------------------------------------------------------------------
        function hws = getHardwareSettings(self)
            hws = self.hardware_settings;
        end
        
        %------------------------------------------------------------------
        function hws = getGain(self)
            error('Not implemented yet!')
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
        function display(self)
            disp(self)
            disp(['Restriction: [' num2str(self.dims(1)) ' ' num2str(size(self.MultiElectrode.electrodePositions,1)) ']']);
            disp(['Add. chans.: [' num2str(self.dimAddChans(1)) ' ' num2str(self.dimAddChans(2)) ']']);
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, varargin)
            sz = self.size;
            [bb relIdx] = hdsort.filewrapper.hdf5.getBoundingBoxFromIndexing(sz, varargin{:});
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
            X = double(hdsort.filewrapper.hdf5.read_dset(self.dataSets.signal.datasetID, block, offset));
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
            [bb relIdx] = hdsort.filewrapper.hdf5.getBoundingBoxFromIndexing(sz, varargin{:});
            
            block = bb(2,:) - bb(1,:) +1;
            numEl = prod(block);
            assert(numEl < 10^9, sprintf('Loading that much data (%d = %s) at once is not recommended!', numEl, sprintf('%d x ', block)));
            offset = bb(1,:) - 1;
            % read outer bounding box of requested data, this is usually faster
            % than reading individual rows or indices
            X = double(hdsort.filewrapper.hdf5.read_dset(self.dataSets.additional_channels.datasetID, block, offset));
            % select actually requested data
            if ~isempty(relIdx)
                X = X(relIdx{:});
            end
        end
        
        %------------------------------------------------------------------
        function [frameNo] = getFrameNumbers(self, varargin)
            sz = self.size; sz(2) = 1;
            
            [bb relIdx] = hdsort.filewrapper.hdf5.getBoundingBoxFromIndexing(sz, varargin{:});
            assert(~any(bb(:)<=0) && ~any(any(bb>repmat(sz(1:size(bb,2)),size(bb,1),1))), 'Indexing out of bounds!');
            block = bb(2,1) - bb(1,1) +1;
            offset = bb(1,1) - 1;
            
            try
            frameNo = double(hdsort.filewrapper.hdf5.read_dset(self.dataSets.frame_numbers.datasetID, block, offset));
            catch
            frameNo = double(hdsort.filewrapper.hdf5.read_dset(self.dataSets.frame_numbers.datasetID, [block, 1], [offset 0]));
            end
            frameNo = frameNo(:)';
            
            % select actually requested data
            if ~isempty(relIdx)
                frameNo = frameNo(relIdx{:});
            end
        end
        
        %------------------------------------------------------------------
        function [missingFrames] = getMissingFrameNumbers(self)
            if ~isstruct(self.missingFrames)
                self.missingFrames = hdsort.filewrapper.util.getMissingFrameNumbers(self.getFrameNumbers());
            end
            missingFrames = self.missingFrames;
        end
        
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = self.dims(1);
        end
        
        %------------------------------------------------------------------
        function fileName = getFullFileName(self)
            fileName = self.fileName;
        end
        
        function X = getFilteredData(self, varargin)
            hpf = 300;
            lpf = 7000;
            fir_filterOrder = 110;
            b  = hdsort.util.filter_design_fir(hpf, lpf, self.getSampleRate(), fir_filterOrder);
            
            X = self.getData(varargin{:});
            
            X = double(X);
            X = bsxfun(@minus, X, mean(X,1));
            
            
            X = conv2(X, b(:), 'same');
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
        
        
        %------------------------------------------------------------------
        function [outFile, P] = preprocessFile(self, path, varargin)
            % Input: path where the preprocessed file is to be stored.
            
            P.outFileName = [];
            P.override = 0;
            P.hpf = 300;
            P.lpf = 7000;
            P.fir_filterOrder = 110;
            P.deflation = 1;
            P.chunkLen = [];
            P.downSample = 1;
            P.chunkIndividualChannels = 0;
            P.subtractMeanOverAllChannels = true;
            P.restrictToTimePeriod = [];
            P.sessionName = '/Sessions/Session0';
            P.save_as_binary = false; %true;
            P.restrictToChannels = [];
            P.chunkSize = 250000;
            P.chunkOverlap = 111;
            P.progressDisplay = 'console';
            P = hdsort.util.parseInputs(P, varargin);
            
            if isempty(P.outFileName)
                P.outFileName = [self.getRootFileName() '.h5'];
            end
            outFile = fullfile(path, P.outFileName);
            
            if P.override
                delete(outFile)
            end
            
            %assert(~exist(outFile, 'file'), ['Output File does already exist! ' outFile]);
            if exist(outFile, 'file')
                warning(['Output File does already exist! ' outFile]);
                return;
            end
            
            assert(int32(P.downSample)==P.downSample, 'downSample must be an integer!');
            assert(P.downSample==1 || P.lpf >0, 'You must high pass filter before dowsampling, otherwise you will get artefacts at the endpoints!');
            
            if isempty(P.restrictToChannels)
                self.restrictToConnectedChannels();
            else
                self.restrictToChannels();
                self.restrictToChannels(P.restrictToChannels);
            end
            
            fileInfo = struct();
            if ~isempty(self.chip_number)
                fileInfo.chipid = self.chip_number;
            else
                fileInfo.chipid = -1;
            end
            
            fileInfo.gain1 = 1;
            fileInfo.gain2 = 1;
            fileInfo.gain3 = 1;
            fileInfo.adc_resolution = 13.8;
            fileInfo.adc_range = 1;
            
            nSamples = self.getNSamples;
            
            if ~isempty(P.restrictToTimePeriod)
                error('Not Implemented At the moment');
                assert(length(P.restrictToTimePeriod) == 2, 'P.restrictToTimePeriod must be a time period in samples [start end]!');
                assert(~any(round(P.restrictToTimePeriod) ~= P.restrictToTimePeriod), 'P.restrictToTimePeriod must be in intergers!');
                assert(P.restrictToTimePeriod(1) > 0 && P.restrictToTimePeriod(2) < nSamples, 'P.restrictToTimePeriod out of bounds!');
                assert(P.restrictToTimePeriod(1) < P.restrictToTimePeriod(2), 'P.restrictToTimePeriod must be increasing');
            end
            
            map = struct();
            map.noChans = numel(self.MultiElectrode.parentElectrodeIndex);
            
            nC = map.noChans;
            assert(length(nC) == 1, 'map.noChans must be a number!')
            assert(nC>0, 'map.noChans must be greater 0!')
            
            % Make sure the map vectors are all in the right format
            map.elNo  = self.MultiElectrode.electrodeNumbers;
            map.chans = self.MultiElectrode.parentElectrodeIndex;
            map.mposx  = self.MultiElectrode.electrodePositions(:, 1);
            map.mposy  = self.MultiElectrode.electrodePositions(:, 2);
            
            map.elNo  = map.elNo(:)';
            map.chans = map.chans(:)';
            map.mposx = map.mposx(:)';
            map.mposy = map.mposy(:)';
            
            fileInfo.sr = self.samplesPerSecond;
            fileInfo.version = -1;
            warning('here!')
            
            filterType = 'None      ';
            
            if P.downSample > 1
                fileInfo.sr = fileInfo.sr/P.downSample;
                assert(P.lpf <= .5*fileInfo.sr, 'The lowpass filter must be below the Nyquist Frequency (AFTER downsampling, if enabled)!');
            end
            
            % Prefiltered
            gainmultiplier = 256;
            h5Type = 'H5T_NATIVE_SHORT';
            dfun = @(x) int16(gainmultiplier*x);
            
            % create filter object
            d = P.downSample;
            
            % INIT FIR Filter
            b  = hdsort.util.filter_design_fir(P.hpf, P.lpf, fileInfo.sr, P.fir_filterOrder);
            filterType = 'FIR fircls';
            
            h = P.hpf;
            l = P.lpf;
            fo = P.fir_filterOrder;
            
            
            % Set File as being in process
            proc = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
            proc(1,1) = int32(1);
            
            % Prepare the spike detection if switched on
            %spikeDetMatrix = prepareSpikeDetection();
            
            % Set the other variables in the file
            prepareFilterVariables();
            prepareMeaVariables();
            
            % Set the two channel lists
            isConnected = prepareChannelList();
            nC_effective = sum(isConnected);
            assert(nC_effective == count, 'Must be identical');
            prepareNewChannelList(isConnected);
            
            if isempty(P.chunkLen) || (~isempty(P.chunkLen) && P.chunkLen == 0)
                % This disables chunking and deflation
                chunkDims = [];
            else
                if P.chunkIndividualChannels
                    chunkDims = [P.chunkLen 1];
                else
                    chunkDims = [P.chunkLen nC_effective];
                end
            end
            
            isConnectedIdx = find(isConnected);
            
            if isempty(P.restrictToTimePeriod)
                lastSample = 0;
            else
                nSamples = P.restrictToTimePeriod(2) - P.restrictToTimePeriod(1) +1;
                lastSample = P.restrictToTimePeriod(1)-1;
            end
            
            dims = [nSamples nC_effective];
            maxDims = [nSamples nC_effective];
            
            if ~P.save_as_binary
                sig = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/sig'], ...
                    dims, maxDims, h5Type, chunkDims, P.deflation);
            else
                
                [pathstr,name,ext] = fileparts(outFile);
                P.binFile = fullfile(pathstr, [name, '.dat']);
                
                % create a binary file where the data is stored
                sig = hdsort.filewrapper.util.BinaryFileMatrix(P.binFile, [1 nC_effective], 'writable', true);
                
                % Save a link to the binary file into /sig:
                binFileName = [name, '.dat'];
                sig_link = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/sig'], [1 length(binFileName) ], [1 length(binFileName) ], 'H5T_C_S1');
                sig_link(1,1:length(binFileName)) = binFileName;
                clear sig_link;
                
                bin_dims = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/bin_dims'], [1 2], [1 2], 'H5T_NATIVE_LONG');
                bin_dims(1, :) = [nSamples nC_effective];
                clear bin_dims;
            end
            
            chunker = hdsort.util.Chunker(nSamples, 'chunkSize', P.chunkSize, ...
                'chunkOverlap', P.chunkOverlap, 'one_sided_overlap', false, 'progressDisplay', P.progressDisplay);
            
            while chunker.hasNextChunk()
                [chunkS] = chunker.getNextChunk();
                
                % Indices of the whole chunks (with the artifacts):
                firstOverlapIdx = chunker.chunk_start_overlap + lastSample;
                lastOverlapIdx  = chunker.chunk_end_overlap + lastSample;
                
                % Indices of kept chunks (without the artifacts):
                firstFullIdx = chunker.chunk_start + lastSample;
                lastFullIdx  = chunker.chunk_end + lastSample;
                
                % Indices within X:
                firstXdiff = chunker.chunk_start - chunker.chunk_start_overlap;
                lastXdiff = chunker.chunk_end_overlap - chunker.chunk_end;
                
                X = self.getData(firstOverlapIdx:lastOverlapIdx, isConnectedIdx);
                
                if P.downSample > 1
                    X = double(X);
                    X = bsxfun(@minus, X, mean(X,1));
                    if P.subtractMeanOverAllChannels
                        X = X-repmat(mean(X,2), 1, size(X,2));
                    end
                    X = resample(X, 1, P.downSample);
                end
                
                X = double(X);
                X = bsxfun(@minus, X, mean(X,1));
                if P.subtractMeanOverAllChannels
                    X = X-repmat(mean(X,2), 1, size(X,2));
                end
                
                X = conv2(X, b(:), 'same');
                
                X = dfun(X);
                sig(firstFullIdx:lastFullIdx,:) = X((1+firstXdiff):(end-lastXdiff), :);
            end
            
            % FRAME NUMBERS
            % write the missing framenumbers to the file. We write the missing
            % framenumbers as a matrix with two columns.
            % The first column is the FIRST INDEX INTO THE DATA for which the
            % framenumbers must be corrected and the second column is the number of
            % missing frames AT THIS POSITION that need to be added
            missingFrames = self.getMissingFrameNumbers();
            
            ffn = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/frame_numbers/first_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ffn(1,1) = int32(missingFrames.first);
            clear ffn
            ffn = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/frame_numbers/last_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ffn(1,1) = int32(missingFrames.last);
            clear ffn
            ffn = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/frame_numbers/dataDims'], [1 2], [1 2], 'H5T_NATIVE_INT');
            ffn(1,:) = int32(self.dims);
            clear ffn
            if missingFrames.n > 0
                lL = numel(missingFrames.begin);
                ffn = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/frame_numbers/missing_fns'], [lL 2], [lL 2], 'H5T_NATIVE_INT');
                ffn(:,:) = [missingFrames.begin(:) missingFrames.length(:)];
            else
                ffn = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/frame_numbers/missing_fns'], [1 1], [1 1], 'H5T_NATIVE_INT');
                ffn(1,1) = int32(-1);
            end
            clear ffn
            
            disp('Done preprocessing.')
            clear sig
            
            % Set File as being done
            proc(1,1) = int32(0);
            clear proc
            
            % ---------------------------------------------------------------------
            function prepareFilterVariables()
                pref = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/prefiltered'], [1 1], [1 1], 'H5T_NATIVE_INT');
                pref(1,1) = int32(1);
                clear pref
                high = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/highpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
                high(1,1) = int32(h);
                clear high
                low = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/lowpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
                low(1,1) = int32(l);
                clear low
                down = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/downsamplefactor'], [1 1], [1 1], 'H5T_NATIVE_INT');
                down(1,1) = int32(d);
                clear down
                ord = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/order'], [1 1], [1 1], 'H5T_NATIVE_INT');
                ord(1,1) = int32(fo);
                clear type
                ord = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/type'], [1 20], [1 20], 'H5T_C_S1');
                ord(1,1:length(filterType)) = filterType;
                clear ord
                gd = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/filter/gainmultiplier'], [1 1], [1 1], 'H5T_NATIVE_INT');
                gd(1,1) = int32(gainmultiplier);
                clear gd
            end
            % ---------------------------------------------------------------------
            function prepareMeaVariables()
                % SAVE THE SOURCEFILES
                ord = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/source_files/raw_h5'], [1 length(self.fileName)], [1 length(self.fileName)], 'H5T_C_S1');
                ord(1,1:length(self.fileName)) = self.fileName;
                clear ord
                
                % CHIP ID
                chipid = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/chipid'], [1 1], [1 1], 'H5T_NATIVE_INT');
                chipid(1,1) = int32(fileInfo.chipid);
                clear chipid
                
                % GAIN
                gain = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/gain'], [1 4], [1 4], 'H5T_NATIVE_DOUBLE');
                gain(1,2:4) = [fileInfo.gain1 fileInfo.gain2 fileInfo.gain3];
                gain(1,1) = prod(gain(1,2:4)); % total gain
                clear gain
                
                % ADC range and resolution
                gain = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/adc_resolution'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
                gain(1,1) = fileInfo.adc_resolution;
                gain = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/adc_range'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
                gain(1,1) = fileInfo.adc_range;
                clear gain;
                
                % SR
                sr_ = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/sr'], [1 1], [1 1], 'H5T_NATIVE_INT');
                sr_(1,1) = int32(fileInfo.sr);
                clear sr_
                
                % VERSION
                version = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/version'], [1 1], [1 1], 'H5T_NATIVE_INT');
                version(1,1) = int32(fileInfo.version);
                clear version
            end
            
            % ---------------------------------------------------------------------
            function isConnected = prepareChannelList()
                % CHANNEL LIST
                names = {'channel_nr', 'connected', 'x', 'y', 'idx', 'dummy', 'damaged'};
                type_id = H5T.create('H5T_COMPOUND',length(names)*32);
                for ii=1:length(names)
                    H5T.insert(type_id, names{ii}, (ii-1)*32, 'H5T_NATIVE_INT');
                end
                cl = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/channel_list'], nC, nC, type_id);
                % get one object
                a = cl(1);
                isConnected = zeros(nC,1);
                count = 0;
                for ii=1:nC
                    % set values into object
                    a.channel_nr = int32(map.elNo(ii));
                    if a.channel_nr > 1
                        a.connected = int32(1);
                        isConnected(ii) = 1;
                        a.x = int32(map.mposx(ii));
                        a.y = int32(map.mposy(ii));
                        a.idx = int32(map.chans(ii));
                        a.dummy = int32(0);
                        a.damaged = int32(0);
                        
                        % set correct index in file to values of object
                        count = count+1;
                        cl(count)=a;
                    else
                        a.connected = int32(0);
                        a.x = int32(0);
                        a.y = int32(0);
                        a.idx = int32(0);
                        a.dummy = int32(0);
                        a.damaged = int32(0);
                    end
                end
                clear a cl
            end
            
            % ---------------------------------------------------------------------
            function prepareNewChannelList(isConnected)
                % CHANNEL LIST IN New CL Format
                idx = isConnected==1;
                idx = idx(:);
                x = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/channel_nr'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
                x(1,1:nC_effective) = map.elNo(idx);
                clear x
                x = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/channel_posx'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
                x(1,1:nC_effective) = map.mposx(idx);
                clear x
                x = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/channel_posy'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
                x(1,1:nC_effective) = map.mposy(idx);
                clear x
                x = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, [P.sessionName '/channel_connected'], [1 nC_effective], [1 nC_effective], 'H5T_NATIVE_DOUBLE');
                x(1,1:nC_effective) = ones(1,nC_effective);
                clear x
            end
        end
        
        
    end
end