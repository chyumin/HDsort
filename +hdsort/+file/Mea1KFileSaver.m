classdef Mea1KFileSaver < handle
    properties
        outFilePath
        outFileNameH5
        outFileNameBinary
        dims
        P
        
        h5infos
        handles
        nextIdx
        
        srate
        gainmultiplier
        elX
        elY
        elNo
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Mea1KFileSaver(path, fname, ...
                                             h5type, dims, maxDims, chunkDims, deflation, ...
                                             gainmultiplier, srate, elX, elY, elNo, varargin)
            P.deflation = 1;
            P.forceFileDeletionIfExists = false;
            P.sessionName = '/Sessions/Session0';
            P.h5path = []; % default is '/sig' (see below)
            P.save_as_binary = true;
            P.self = [];
            P.amplifierGains = [1 1 1];
            P.adc_resolution = 13.8;
            P.adc_range = 1;            
            P.fileVersion = -1;
            P = hdsort.util.parseInputs(P, varargin);
            
            if isempty(P.h5path)
                P.h5path = '/sig';
            end
            
            self.P = P;
            self.outFilePath = path;
            if ~exist(path, 'file') && ~isempty(path)
                mkdir(path);
            end
            self.outFileNameH5 = fullfile(path, fname);
            
            [path, name, ext] = fileparts(self.outFileNameH5);
            if isempty(ext)
                self.outFileNameH5 = [self.outFileNameH5 '.h5'];
            end
            self.outFileNameBinary = fullfile(path, [name '.dat']);
            
            if P.forceFileDeletionIfExists
                if exist(self.outFileNameBinary, 'file')
                    delete(self.outFileNameBinary);
                end
                if exist(self.outFileNameH5, 'file')
                    delete(self.outFileNameH5);
                end
            end
            
            self.dims = dims(:)';
            
            assert(numel(dims) == 2, 'dims must have 2 dimensions!');
            assert(~any(dims<1), 'dims must be >0');
            assert(~any(round(dims) ~= dims), 'dims must integer');
            assert(~exist(self.outFileNameH5, 'file'), ['Output File does already exist! ' self.outFileNameH5]);
            assert(~exist(self.outFileNameBinary, 'file'), ['Output File does already exist! ' self.outFileNameBinary]);
            
            self.srate = srate;
            self.elX = elX;
            self.elY = elY;
            self.elNo = elNo;
            self.gainmultiplier = gainmultiplier;
            
            % Make sure the self vectors are all in the right format
            self.elNo  = self.elNo(:)';
%             self.chans = self.chans(:)';
            self.elX = self.elX(:)';
            self.elY = self.elY(:)';
            
            assert(~any(round(maxDims) ~= maxDims), 'maxDims must be integer!');
            
            self.h5infos.h5type = h5type;
            self.h5infos.dims = dims;
            self.h5infos.maxDims = maxDims;
            self.h5infos.chunkDims = chunkDims;
            self.h5infos.deflation = deflation;   
            
            self.nextIdx = 1;
            self.initFiles();
        end
        
        % -----------------------------------------------------------------
        function initFiles(self)
            % Set File as being in process
            self.handles.proc = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, '/bFileIsInProcess', [1 1], [1 1], 'H5T_NATIVE_INT');
            self.handles.proc(1,1) = int32(1);     
            
            if ~self.P.save_as_binary
                self.handles.sig = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/sig'], ...
                    self.h5infos.dims, self.h5infos.maxDims, self.h5infos.h5type, self.h5infos.chunkDims, self.h5infos.deflation);
            else
                % create a binary file where the data is stored
                self.handles.sig = hdsort.file.util.BinaryFileMatrix(self.outFileNameBinary, [1, self.h5infos.maxDims(2)], 'writable', true);
                
                % Save a link to the binary file into /sig:
                [path, name, ext] = fileparts(self.outFileNameBinary);
                fileNameBinary = [name ext];
                
                sig_link = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/sig'], [1 length(fileNameBinary) ], [1 length(fileNameBinary) ], 'H5T_C_S1');
                sig_link(1,1:length(fileNameBinary)) = fileNameBinary;
                clear sig_link;
                bin_dims = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/bin_dims'], [1 2], [1 2], 'H5T_NATIVE_LONG');
                bin_dims(1, :) = int32(self.h5infos.maxDims);
                clear bin_dims;                
            end
            
            % SR
            sr_ = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/sr'], [1 1], [1 1], 'H5T_NATIVE_INT');
            sr_(1,1) = int32(self.srate);
            clear sr_    
            
            % VERSION
            version = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/version'], [1 1], [1 1], 'H5T_NATIVE_INT');
            version(1,1) = int32(self.P.fileVersion);
            clear version
            
            % GAIN
            gain = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/gain'], [1 4], [1 4], 'H5T_NATIVE_DOUBLE');
            gain(1,2:4) = self.P.amplifierGains;
            gain(1,1) = prod(self.P.amplifierGains); % total gain
            clear gain

            % ADC range and resolution
            gain = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/adc_resolution'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
            gain(1,1) = self.P.adc_resolution;
            gain = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/adc_range'], [1 1], [1 1], 'H5T_NATIVE_DOUBLE');
            gain(1,1) = self.P.adc_range;
            clear gain;         
            
            self.saveNewChannelList();
        end
        
        % -----------------------------------------------------------------
        function finalizeFile(self)
            % Set File as being done
            self.handles.proc(1,1) = int32(0);
            % clear handles (not really necessary)
            self.handles = [];
        end
        
        % -----------------------------------------------------------------
        function saveChunk(self, X)
            if isempty(self.handles)
                error('File was already finalized!');
            end
            self.handles.sig(self.nextIdx:self.nextIdx+size(X,1)-1,:) = self.gainmultiplier*X;
            self.nextIdx = self.nextIdx+size(X,1);     
            if self.nextIdx == self.h5infos.maxDims(1)+1
                self.finalizeFile();
            end
        end
        
        % -----------------------------------------------------------------
        function saveChunkTransposed(self, Xt)
            % This function is faster with X is really large than the
            % non-transposed one. Use this if you can transpose X before
            % calling this function by calling X=X';
            if isempty(self.handles)
                error('File was already finalized!');
            end
            self.handles.sig.appendDataTransposed(self.gainmultiplier*Xt);
            self.nextIdx = self.nextIdx + size(Xt,2);     
            if self.nextIdx == self.h5infos.maxDims(1)+1
                self.finalizeFile();
            end
        end        
        
        % -----------------------------------------------------------------
        function saveFramenumbers(self, firstFN, lastFN, missingIDX, missingLength)
            % FRAME NUMBERS
            % write the missing framenumbers to the file. We write the missing
            % framenumbers as a matrix with two columns.
            % The first column is the FIRST INDEX INTO THE DATA for which the
            % framenumbers must be corrected and the second column is the number of
            % missing frames AT THIS POSITION that need to be added
            missingLength = missingLength(:);
            missingIDX = missingIDX(:);
            L = length(missingIDX);
            assert(L == length(missingLength), 'Must be same!');
           
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/first_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ffn(1,1) = int32(firstFN);
            clear ffn
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/last_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ffn(1,1) = int32(lastFN);
            clear ffn
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/dataDims'], [1 2], [1 2], 'H5T_NATIVE_INT');
            ffn(1,:) = int32(self.h5infos.maxDims);
            clear ffn
            if ~isempty(missingIDX)
                ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/missing_fns'], [L 2], [L 2], 'H5T_NATIVE_INT');
                ffn(:,:) = [missingIDX(:) missingLength(:)];
            else
                ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/missing_fns'], [1 1], [1 1], 'H5T_NATIVE_INT');
                ffn(1,1) = int32(-1);
            end
            clear ffn
        end
        
        
        % ---------------------------------------------------------------------
        function saveFilterSettings(self, filter)
            % INPUT FORMAT:
            % filter.prefiltered
            % filter.highpass
            % filter.lowpass
            % filter.downsamplefactor
            % filter.order
            % filter.gainmultiplier;
            % filter.type (string)
            
            %% Set filter info:
            pref = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/prefiltered'], [1 1], [1 1], 'H5T_NATIVE_INT');
            pref(1,1) = int32(filter.prefiltered);
            clear pref
            high = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/highpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
            high(1,1) = int32(filter.highpass);
            clear high
            low = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/lowpass'], [1 1], [1 1], 'H5T_NATIVE_INT');
            low(1,1) = int32(filter.lowpass);
            clear low
            down = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/downsamplefactor'], [1 1], [1 1], 'H5T_NATIVE_INT');
            down(1,1) = int32(filter.downsamplefactor);
            clear down
            ord = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/order'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ord(1,1) = int32(filter.order);
            clear type
            
            L = length(filter.type);
            ord = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/type'], [1 L], [1 L], 'H5T_C_S1');
            ord(1,1:L) = filter.type;
            clear ord
            
            gd = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/filter/gainmultiplier'], [1 1], [1 1], 'H5T_NATIVE_INT');
            gd(1,1) = int32(self.gainmultiplier);
            clear gd
            
        end
        
        % -----------------------------------------------------------------
        function saveDetectedSpikes(self, spikeDetection)
            gdf_spikesDetectedUp = hdsort.spiketrain.toGdf( spikeDetection.spikesDetectedUp, 1:size(spikeDetection.pks_down, 1));
            gdf_pks_up = hdsort.spiketrain.toGdf( spikeDetection.pks_up, 1:size(spikeDetection.pks_down, 1));
            gdf_spikesDetectedDown = hdsort.spiketrain.toGdf( spikeDetection.spikesDetectedDown, 1:size(spikeDetection.pks_down, 1));
            gdf_pks_down = hdsort.spiketrain.toGdf( spikeDetection.pks_down, 1:size(spikeDetection.pks_down, 1));
            
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/spikesDetectedUp'], size(gdf_spikesDetectedUp), size(gdf_spikesDetectedUp), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = gdf_spikesDetectedUp;
            clear ffn
            
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/pks_up'], size(gdf_pks_up), size(gdf_pks_up), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = gdf_pks_up;
            clear ffn
            
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/spikesDetectedDown'], size(gdf_spikesDetectedDown), size(gdf_spikesDetectedDown), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = gdf_spikesDetectedDown;
            clear ffn
            
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/pks_down'], size(gdf_pks_down), size(gdf_pks_down), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = gdf_pks_down;
            clear ffn
        
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/minDist'], size(spikeDetection.minDist), size(spikeDetection.minDist), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = spikeDetection.minDist;
            clear ffn
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/threshold'], size(spikeDetection.threshold), size(spikeDetection.threshold), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = spikeDetection.threshold;
            clear ffn
        end
        
        % ---------------------------------------------------------------------
        function saveNoiseStd(self, noiseStd)
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/spike_detection/noiseStd'], size(noiseStd), size(noiseStd), 'H5T_NATIVE_DOUBLE');
            ffn(:,:) = noiseStd;
            clear ffn
        end

        % ---------------------------------------------------------------------
        function saveRawH5File(self, ffile)
            ord = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/source_files/raw_h5'], [1 length(ffile)], [1 length(ffile)], 'H5T_C_S1');
            ord(1,1:length(ffile)) = ffile;
            clear ord
        end
        
        % ---------------------------------------------------------------------    
        function saveselfpingFile(self, mfile)
            ord = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/source_files/selfping_file'], [1 length(mfile)], [1 length(mfile)], 'H5T_C_S1');
            ord(1,1:length(mfile)) = mfile;
            clear ord            
        end
        
        % ---------------------------------------------------------------------    
        function saveChipID(self, chipid)
            % CHIP ID
            h = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/chipid'], [1 1], [1 1], 'H5T_NATIVE_INT');
            h(1,1) = int32(chipid);
            clear h   
        end        

        % -----------------------------------------------------------------
        function saveMissingFrames(self, missingFrames)
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/first_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ffn(1,1) = int32(missingFrames.first);
            clear ffn
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/last_fn'], [1 1], [1 1], 'H5T_NATIVE_INT');
            ffn(1,1) = int32(missingFrames.last);
            clear ffn
            ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/dataDims'], [1 2], [1 2], 'H5T_NATIVE_INT');
            ffn(1,:) = int32(self.dims);
            clear ffn
            if missingFrames.n > 0
                lL = numel(missingFrames.begin);
                ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/missing_fns'], [lL 2], [lL 2], 'H5T_NATIVE_INT');
                ffn(:,:) = [missingFrames.begin(:) missingFrames.length(:)];
            else
                ffn = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/frame_numbers/missing_fns'], [1 1], [1 1], 'H5T_NATIVE_INT');
                ffn(1,1) = int32(-1);
            end
            clear ffn
        end
        
        % ---------------------------------------------------------------------
        function saveNewChannelList(self)
            nC = self.h5infos.maxDims(2);
            % CHANNEL LIST IN New CL Format
            x = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/channel_nr'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
            x(1,:) = self.elNo;
            clear x
            x = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/channel_posx'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
            x(1,:) = self.elX;
            clear x
            x = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/channel_posy'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
            x(1,:) = self.elY;
            clear x
            x = hdsort.file.hdf5.createVariableAndOrFile(self.outFileNameH5, [self.P.sessionName '/channel_connected'], [1 nC], [1 nC], 'H5T_NATIVE_DOUBLE');
            x(1,:) = ones(1, nC);
            clear x
        end
    end
end