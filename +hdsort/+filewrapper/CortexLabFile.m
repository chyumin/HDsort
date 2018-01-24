classdef CortexLabFile < hdsort.filewrapper.SingleFileWrapper
    % Wrapper for files in http://phy.cortexlab.net/data/sortingComparison/
    
    properties
        P
        prmFile
        prbFile
        
        fileInfo
        
        memmap
        className
        
        
        nDims
        dims
        connectedChannel
        
        missingFrames
        
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = CortexLabFile(fileName, prmFile, varargin)
            assert(exist(fileName)==2, ['File ' fileName ' does not exist!'])
            assert(exist(prmFile)==2, ['Prm file ' prmFile ' does not exist!'])
            
            %% Read Parameter file:
            prm = hdsort.filewrapper.CortexLabFile.readPRM(prmFile);
             
            samplingRate_ = prm.sample_rate;
            nChannels = prm.n_channels;
            
            
            %% Construct MultiElectrode:
            [path] = fileparts(fileName);
            [~, prbBaseName] = fileparts(prm.prb_file)
            prbFile = fullfile(path, [prbBaseName '.prb.mat']);
            prb = load(prbFile);
            if 0
                allRecordedChannels = (1:nChannels)';
                mappedChannels = double(prb.channels(:))+1;
                allUnsortedChannels = [mappedChannels; allRecordedChannels(~ismember(allRecordedChannels, mappedChannels))];
                [allSortedChannels, sortidx] = sort(allUnsortedChannels);
                
                xPos = [double(prb.geometry(:,2)); zeros(nChannels - size(prb.geometry, 1), 1)-1 ];
                yPos = [double(prb.geometry(:,3)); zeros(nChannels - size(prb.geometry, 1), 1)-1 ];
                elns = [double(prb.geometry(:,1)); zeros(nChannels - size(prb.geometry, 1), 1)-1 ];
                
                multiElectrode_ = hdsort.filewrapper.MultiElectrode([xPos(sortidx) yPos(sortidx)], elns(sortidx));
            else
                % Create a vector for each channel, starting at 1 (Matlab).
                % Then, map each channel from prb to this vector (+1
                % because it starts at 0).
                % Finally sort this vector and use the sorting index to
                % sort the electrode positions. The unknown channels are at
                % [-1,-1]
                allRecordedChannels = (1:nChannels)';
                mappedChannels = prb.geometry(:,1)+1;
                allUnsortedChannels = [mappedChannels; allRecordedChannels(~ismember(allRecordedChannels, mappedChannels))];
                [allSortedChannels, sortidx] = sort(allUnsortedChannels);
                
                xPos = [double(prb.geometry(:,2)); zeros(nChannels - size(prb.geometry, 1), 1)-1 ];
                yPos = [double(prb.geometry(:,3)); zeros(nChannels - size(prb.geometry, 1), 1)-1 ];
                elns = [double(prb.geometry(:,1)); zeros(nChannels - size(prb.geometry, 1), 1)-1 ];
                multiElectrode_ = hdsort.filewrapper.MultiElectrode([xPos(sortidx) yPos(sortidx)], elns(sortidx));
            end
            
            
            self = self@hdsort.filewrapper.SingleFileWrapper('CortexLabFile', samplingRate_, multiElectrode_, fileName, varargin{:});
             
            %%
            self.fileName = fileName;
            self.prbFile = prbFile;
            self.fileInfo = dir(self.fileName);
            
            
            self.dims = [nChannels, (self.fileInfo.bytes/nChannels/2)];
            self.memmap = memmapfile(self.fileName, 'Format', {'int16', self.dims, 'X'});
            
            self.connectedChannel = find(self.MultiElectrode.electrodeNumbers > -1);
            
        end
         
        %------------------------------------------------------------------
        function display(self)
            disp(self)
            disp(['Restriction: [' num2str(self.dims(2)) ' ' num2str(size(self.MultiElectrode.electrodePositions,1)) ']']);
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, varargin)
            if length(varargin) == 2
                c = {varargin{[2 1]}};
            else
                c = varargin;
            end
            X =  self.memmap.Data.X(c{:})';
        end
        
        %------------------------------------------------------------------
        function [frameNo] = getFrameNumbers(self, varargin)
            frameNo = 1:self.getNSamples();
        end
       
        %------------------------------------------------------------------
        function [missingFrames] = getMissingFrameNumbers(self, bufferFile)
            if nargin == 2 && exist(bufferFile)
                load(bufferFile)
            else
                if ~isstruct(self.missingFrames)
                    missingFrames = hdsort.filewrapper.util.getMissingFrameNumbers(self.getFrameNumbers());
                end
                if nargin == 2
                    save(bufferFile, 'missingFrames');
                end
            end
            self.missingFrames = missingFrames;
        end
        
        %------------------------------------------------------------------
        function L = getNSamples_(self)
            L = self.dims(2);
        end
        
        %------------------------------------------------------------------
        function fileName = getFullFileName(self)
            fileName = self.fileName;
        end
        
        function X = getFilteredData(self, varargin)
            hpf = 300;
            lpf = 7000;
            fir_filterOrder = 110;
            fir_filter  = hdsort.util.filter_design_fir(hpf, lpf, self.getSampleRate(), fir_filterOrder);
            
            X = self.getData(varargin{:});
            X = double(X);
            X = bsxfun(@minus, X, mean(X,1));
            X = conv2(X, fir_filter(:), 'same');
        end
        
        %------------------------------------------------------------------
        function rootFileName = getRootFileName(self)
            [path_ name_ ext_] = fileparts(self.fileName);
            [path_ rootFileName ext_] = fileparts(name_);
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
            P.bUseFPGA_IIR_Filter = false;
            P.deflation = 1;
            P.chunkLen = [];
            P.downSample = 1;
            P.chunkIndividualChannels = 0;
            P.subtractMeanOverAllChannels = true;
            P.restrictToTimePeriod = [];
            P.sessionName = '/Sessions/Session0';
            P.save_as_binary = true;
            P.restrictToChannels = [];
            P.chunkSize = 250000;
            P.chunkOverlap = 111;
            P.progressDisplay = 'console';
            P = hdsort.util.parseInputs(P, varargin);
            
            if isempty(P.outFileName)
                P.outFileName = ['pre_' self.getRootFileName() '.h5'];
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
            fileInfo.chipid = -1;
            
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
            gainmultiplier = 1;
            dfun = @(x) x;
            
            % create filter object
            d = P.downSample;
            
            % INIT FIR Filter
            fir_filter  = hdsort.util.filter_design_fir(P.hpf, P.lpf, fileInfo.sr, P.fir_filterOrder);
            filterType = 'FIR fircls';
            
            h = P.hpf;
            l = P.lpf;
            fo = P.fir_filterOrder;
            
            % Set File as being in process
            proc = hdsort.filewrapper.hdf5.createVariableAndOrFile(outFile, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
            proc(1,1) = int32(1);
            
            % Set the other variables in the file
            prepareFilterVariables();
            prepareMeaVariables();
            
            % Set the two channel lists
            [isConnected, count] = prepareChannelList();
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
                sig = hdsort.filewrapper.util.binaryFileMatrix(P.binFile, [1 nC_effective], 'writable', true);
                
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
                X = conv2(X, fir_filter(:), 'same');
                
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
            
            if 0
                %% TEST SCENARIO FOR FRAMENUMBER
                %FRAMES = [1 2 3 4 5 6 7 8 15 16 17 19 20 21 23]
                FRAMES = 1:50; FRAMES(30) = []; FRAMES(33) = [];
                
                missingFrames = hdsort.filewrapper.util.getMissingFrameNumbers(FRAMES);
                ffn = [missingFrames.begin(:) missingFrames.length(:)];
                FRAMES_RESTORED = hdsort.filewrapper.util.getFrameNumbersFromMissing(missingFrames.first, missingFrames.last, ffn')
            end
            
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
            function [isConnected, count] = prepareChannelList()
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
    
    %% STATIC FUNCTIONS ###################################################
    methods(Static)
        
        %------------------------------------------------------------------
        function [out] = readPRM(prmFile)
            
            text = fileread(prmFile);
            text = strrep(text, '', '');
            text = strrep(text, char(10), '');
            text = strrep(text, char(32), '');
            text = strrep(text, char(9), '');
            
            %%
            tsplit = strsplit(text, {'''', '=', ','});
            out = struct();
            for ii = 1:numel(tsplit)
                
                
                if strcmp(strtrim(tsplit{ii}), 'sample_rate') ...
                        || strcmp(strtrim(tsplit{ii}), 'n_channels') ...
                        
                    
                    key_str = strtrim(tsplit{ii});
                    ii = ii + 1;
                    
                    %%
                    value_str = strtrim(tsplit{ii});
                    value_str = strrep(value_str,':', '');
                    value = str2num(value_str);
                    out.(key_str) = value;
                    
                elseif strcmp(strtrim(tsplit{ii}), 'prb_file') && ~isfield(out, 'prb_file')
                    
                    key_str = strtrim(tsplit{ii});
                    ii = ii + 1;
                    
                    %%
                    value_str = strtrim(tsplit{ii});
                    value_str = strrep(value_str,':', '');
                    out.(key_str) = value_str;
                    
                end
            end
        end
        
      
        
    end
    
    
end