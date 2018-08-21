classdef WaveFormFileMat < hdsort.file.FileWrapperInterface
   
    properties (Access = protected) 
        binFile
        wfs
        nextIdx
        nRowsGdf
        writable
    end
        
    properties
        P
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = WaveFormFileMat(fileName, varargin) % cutLeft, Tf, nC, samplesPerSecond
            
            P.debug = false;
            P.Tf = [];
            P.cutLeft = [];
            P.nC = [];
            P.samplesPerSecond = 20000;
            P.MultiElectrode = [];
            P.precision = 'double';
            P.writable = false;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            self = self@ hdsort.file.FileWrapperInterface('WaveFormFileMat', P.samplesPerSecond, P.MultiElectrode);
            
            self.P = P; clear P
            
            self.wfs = [];
            self.nRowsGdf = 4;
            
            self.name = fileName;
            
            if exist(self.name, 'file') ~= 2
                self.writable = true;
                
                self.buffer.Tf = self.P.Tf;
                self.buffer.cutLeft = self.P.cutLeft;
                self.buffer.nC = self.P.nC;
                self.buffer.samplesPerSecond = self.P.samplesPerSecond;
                self.buffer.precision = self.P.precision;
                
                if self.P.debug
                    disp('Create new empty waveform-file')
                    disp(self.name)
                end
                
                [pathstr,name,ext] = fileparts(self.name);
                self.buffer.binFileName = [name, '.dat'];
                self.binFile = fullfile(pathstr, self.buffer.binFileName);
                
                if exist(self.binFile) == 2
                    delete(self.binFile)
                end
                
                % create a binary file where the data is stored
                self.wfs = hdsort.file.util.BinaryFileMatrix(self.binFile, [1 self.buffer.nC*self.buffer.Tf], ...
                    'writable', true, 'precision', self.buffer.precision);
                
                self.buffer.gdf = [];
                self.buffer.binDims = [0, int32(self.buffer.nC*self.buffer.Tf)];
                
                self.nextIdx = 1;
         
            else
                self.writable = self.P.writable;
                
                self.buffer = load(self.name);
                
                pathstr = fileparts(self.name);
                self.binFile = fullfile(pathstr, self.buffer.binFileName);
                assert( exist(self.binFile, 'file') == 2, ['Binary file missing! ' self.name])
                
                binDims = self.buffer.binDims;
                
                % Open the binary file matrix:
                if binDims(1)
                    self.wfs = hdsort.file.util.BinaryFileMatrix(self.binFile, binDims,...
                        'writable', self.writable, 'precision', self.buffer.precision);
                else
                    self.wfs = hdsort.file.util.BinaryFileMatrix(self.binFile, [1, self.buffer.nC*self.buffer.Tf], ...
                        'writable', self.writable, 'precision', self.buffer.precision);
                end
                self.nextIdx = binDims(1)+1;
            end
            
        end
        
        % DESTRUCTOR
        %------------------------------------------------------------------
        function delete(self)
            self.saveBufferToMatFile();
            self.wfs = [];
        end
        
        %------------------------------------------------------------------
        function saveBufferToMatFile(self)
            buffer = self.buffer;
            save(self.name, '-struct', 'buffer')
        end
        
        %------------------------------------------------------------------
        function display(self)
           disp(self)
        end
        
        %------------------------------------------------------------------
        function w = iswritable(self)
           	w = logical(self.writable);
        end
        
        %------------------------------------------------------------------
        function X = getData_(self, waveformidx, tfchidx)
            binDims = self.buffer.binDims;
            assert( max(waveformidx) <= binDims(1), 'Index exceeding number of waveforms!') 
            assert( max(tfchidx) <= binDims(2), 'Index exceeding length of waveforms!') 
            X = self.wfs(waveformidx, tfchidx);
        end
        
        function X = getData(self, tfidx, chidx, waveformidx)
            binDims = self.buffer.binDims;
            if ~binDims(1)
                X = [];
                return;
            end
            
            %% If exactly two indices are used, give the data back in the vector format:
            if nargin == 3
                waveformidx = tfidx;
                tfchidx = chidx;
                
                if ischar(waveformidx) && strcmp(waveformidx, ':')
                    waveformidx = 1:binDims(1);
                end
                if ischar(tfchidx) && strcmp(tfchidx, ':')
                    tfchidx = 1:(self.buffer.Tf*self.buffer.nC);
                end
                
                X = self.getData_(waveformidx, tfchidx);
                return;
            end
            
            if nargin < 2 tfidx = 1:self.buffer.Tf; end
            if nargin < 3 chidx = 1:self.buffer.nC; end
            if nargin < 4 waveformidx = 1:binDims(1); end
            
            if ischar(tfidx) && strcmp(tfidx, ':')
                tfidx = 1:self.buffer.Tf;
            end
            if ischar(chidx) && strcmp(chidx, ':')
                chidx = 1:self.buffer.nC;
            end
            if ischar(waveformidx) && strcmp(waveformidx, ':')
                waveformidx = 1:binDims(1);
            end
            
            idx2 = false(self.buffer.Tf, self.buffer.nC);
            idx2(tfidx, chidx) = true;
            tfchidx = find(idx2(:)');
            vX = self.getData_(waveformidx, tfchidx);
            
            [nS, TfnC] = size(vX);
            X = reshape(vX', TfnC/numel(chidx), numel(chidx), nS);
        end
        
        %------------------------------------------------------------------
        function gdf = getGdf(self, idx)
            binDims = self.buffer.binDims;
            
            if ~binDims(1)
                gdf = [];
                return
            end
            if nargin < 2 || ( ischar(idx) && strcmp(idx, ':') )
                idx = 1:binDims(1); 
            end
            assert( max(idx) <= binDims(1), 'Index exceeding number of waveforms!')
            
            gdf = self.buffer.gdf(idx, :);
        end
        
        %------------------------------------------------------------------
        function dims = getDims(self)
            dims = [self.getNSamples(), self.getTf() * self.getNChannels()];
        end
        
        %------------------------------------------------------------------
        function X = getWaveform_(self, waveformidx, tfchidx)
            X = self.getData_(waveformidx, tfchidx)
        end
        
        %------------------------------------------------------------------
        function N = getNSamples_(self)
            binDims = self.buffer.binDims;
            N = double(binDims(1));
        end
        
        %------------------------------------------------------------------
        function n = getNChannels(self)
            n = self.buffer.nC;
        end
        
        %------------------------------------------------------------------
        function n = getTf(self)
            n = self.buffer.Tf;
        end
        
        %------------------------------------------------------------------
        function n = getCutLeft(self)
            n = self.buffer.cutLeft;
        end
        
        %------------------------------------------------------------------
        
        function addWaveforms(self, wfs, gdf)
            assert( self.iswritable(), 'File is writable!')
            
            % Reshape input if necessary:
            if ndims(wfs) == 3
                vWFS = hdsort.waveforms.t2v(wfs);
            else
                assert(ndims(wfs) == 2, 'Input waveforms must either be a vector or tensor!')
                vWFS = wfs;
            end
            clear wfs
            
            % Check that input size corresponds to the expected values:
            [N, C] = size(vWFS);
            assert( C == self.buffer.nC * self.buffer.Tf, 'The input waveforms must conform to the filesize!')
            
            % Chack and if necessary reshape the gdf:
            if nargin == 2
                gdf = ones(N, self.nRowsGdf);
            end
            if size(gdf, 2) < self.nRowsGdf 
                gdf = [gdf, ones(N, self.nRowsGdf-size(gdf, 2))]; 
            end
            assert(size(gdf, 1) == N, 'The input gdf must have the same number of spikes as the waveforms!')
            
            a0 = self.nextIdx;
            a1 = self.nextIdx+N-1;
            
            % Save data:
            self.wfs(a0:a1, :) = vWFS;
            self.buffer.binDims = int32([a1, self.buffer.nC*self.buffer.Tf]);
            self.buffer.gdf = [self.buffer.gdf; gdf]; 
            
            assert( isequal(size(self.wfs), self.buffer.binDims), 'Saved binDims do not correspond to size of wfs!')
            self.nextIdx = self.nextIdx + N;
        end
        
    end
    
end
