classdef WaveFormFile < hdsort.file.FileWrapperInterface
    % Wrapper for BEL-MEA-file format.
    
    properties (Access = protected) 
        wfs
        gdf
        binDims
        nextIdx
        
        nRowsGdf
        
        Tf
        cutLeft
        nC
        R
        precision
        
        writable
    end
        
    properties
        P
        fileName
        
        binFile
        binFileName
        
        datasets
    end
    
    methods
        %% CONSTRUCTOR
        %------------------------------------------------------------------
        function self = WaveFormFile(fileName, varargin) % cutLeft, Tf, nC, samplesPerSecond
            
            P.debug = false;
            P.Tf = [];
            P.cutLeft = [];
            P.nC = [];
            P.samplesPerSecond = 20000;
            P.MultiElectrode = [];
            P.precision = 'double';
            P.writable = false;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            self = self@hdsort.file.FileWrapperInterface('WaveFormFile', P.samplesPerSecond, P.MultiElectrode);
            
            self.wfs = [];
            self.datasets.binFile = '/binFile';
            self.datasets.binDims = '/binDims';
            self.datasets.gdf = '/gdf';
            self.datasets.nC = '/nC';
            self.datasets.cutLeft = '/cutLeft';
            self.datasets.samplesPerSecond = '/samplesPerSecond';
            self.datasets.precision = '/precision';
            self.nRowsGdf = 4;
            
            self.fileName = fileName;
            
            if exist(self.fileName) ~= 2
                self.writable = true;
                
                self.Tf = P.Tf;
                self.cutLeft = P.cutLeft;
                self.nC = P.nC;
                self.samplesPerSecond = P.samplesPerSecond;
                self.precision = P.precision;
                
                if P.debug
                    disp('Create new empty waveform-file')
                    disp(self.fileName)
                end
                
                [pathstr,name,ext] = fileparts(self.fileName);
                self.binFileName = [name, '.dat'];
                self.binFile = fullfile(pathstr, self.binFileName);
                
                if exist(self.binFile) == 2
                    delete(self.binFile)
                end
                
                % create a binary file where the data is stored
                self.wfs = hdsort.file.util.BinaryFileMatrix(self.binFile, [1 self.nC*self.Tf], ...
                    'writable', true, 'precision', P.precision);
                
                self.gdf = hdsort.file.hdf5.matrixCreate(self.fileName, self.datasets.gdf, ...
                    [0 self.nRowsGdf], [-1, self.nRowsGdf], 'H5T_NATIVE_DOUBLE', [1, 1], 0);
                
                self.binDims = hdsort.file.hdf5.matrixCreate(self.fileName, self.datasets.binDims, ...
                    [1 2], [-1 2], 'H5T_NATIVE_LONG', [1, 1], 0);
                self.binDims(1, :) = [0, int32(self.nC*self.Tf)];
                
                self.nextIdx = 1;
                
                % Save a link to the binary file into /binFile:
                binFileName = hdsort.file.hdf5.createVariableAndOrFile(self.fileName, self.datasets.binFile, [1 length(self.binFileName) ], [1 length(self.binFileName) ], 'H5T_C_S1');
                binFileName(1,1:length(self.binFileName)) = self.binFileName;
                clear binFileName;
                
                precision_ = hdsort.file.hdf5.createVariableAndOrFile(self.fileName, self.datasets.precision, [1 length(self.precision) ], [1 length(self.precision) ], 'H5T_C_S1');
                precision_(1,1:length(self.precision)) = self.precision;
                clear precision_;
                
                dims = [1 1];
                maxDims = [1 1];
                h5type = 'H5T_NATIVE_INT';
                chunkDims = [1 1];
                deflation = 0;
                
                nC_ = hdsort.file.hdf5.matrixCreate(self.fileName, self.datasets.nC, ...
                    dims, maxDims, h5type, chunkDims, deflation);
                nC_(1,1) = int32(self.nC);
                clear nC_
                
                cutLeft_ = hdsort.file.hdf5.matrixCreate(self.fileName, self.datasets.cutLeft, ...
                    dims, maxDims, h5type, chunkDims, deflation);
                cutLeft_(1,1) = int32(self.cutLeft);
                clear cutLeft_
                
                samplesPerSecond_ = hdsort.file.hdf5.matrixCreate(self.fileName, self.datasets.samplesPerSecond, ...
                    dims, maxDims, h5type, chunkDims, deflation);
                samplesPerSecond_(1,1) = int32(self.samplesPerSecond);
                clear samplesPerSecond_
            else
                self.writable = P.writable;
                
                % Extract necessary information from the h5 file to
                % construct the binaryfile handle:
                self.binFileName = cell2mat( h5read(self.fileName, self.datasets.binFile) );
                [pathstr,name,ext] = fileparts(self.fileName);
                self.binFile = fullfile(pathstr, self.binFileName);
                assert( exist(self.binFile, 'file') == 2, ['Task aborted: binary file ' self.binFile ' not found!']);
                
                % Open the writable h5-matrices:
                self.gdf = hdsort.file.hdf5.matrix(self.fileName, self.datasets.gdf, ~self.writable);
                self.binDims = hdsort.file.hdf5.matrix(self.fileName, self.datasets.binDims, ~self.writable);
                binDims = self.binDims(:,:);
                
                % Read out other parameters of the file:
                self.precision = cell2mat( h5read(self.fileName, self.datasets.precision));
                self.nC = double( h5read(self.fileName, self.datasets.nC));
                self.Tf = double(binDims(2)/self.nC);
                self.cutLeft = double( h5read(self.fileName, self.datasets.cutLeft));
                self.samplesPerSecond = double(h5read(self.fileName, self.datasets.samplesPerSecond));
                
                % Open the binary file matrix:
                if binDims(1)
                    self.wfs = hdsort.file.util.BinaryFileMatrix(self.binFile, binDims,...
                        'writable', self.writable, 'precision', self.precision);
                else
                    self.wfs = hdsort.file.util.BinaryFileMatrix(self.binFile, [1 self.nC*self.Tf], ...
                        'writable', self.writable, 'precision', P.precision);
                end
                self.nextIdx = binDims(1)+1;
                
            end
            
            self.P = P;
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
            binDims = self.binDims(:,:);
            
            assert( max(waveformidx) <= binDims(1), 'Index exceeding number of waveforms!') 
            assert( max(tfchidx) <= binDims(2), 'Index exceeding length of waveforms!') 
            
            X = self.wfs(waveformidx, tfchidx);
        end
        
        function X = getData(self, tfidx, chidx, waveformidx)
            
            binDims = self.binDims(:,:);
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
                    tfchidx = 1:(self.Tf*self.nC);
                end
                
                X = self.getData_(waveformidx, tfchidx);
                return;
            end
            
            if nargin < 2 tfidx = 1:self.Tf; end
            if nargin < 3 chidx = 1:self.nC; end
            if nargin < 4 waveformidx = 1:binDims(1); end
            
            if ischar(tfidx) && strcmp(tfidx, ':')
                tfidx = 1:self.Tf;
            end
            if ischar(chidx) && strcmp(chidx, ':')
                chidx = 1:self.nC;
            end
            if ischar(waveformidx) && strcmp(waveformidx, ':')
                waveformidx = 1:binDims(1);
            end
            
            idx2 = false(self.Tf, self.nC);
            idx2(tfidx, chidx) = true;
            tfchidx = find(idx2(:)');
            vX = self.getData_(waveformidx, tfchidx);
            
            [nS, TfnC] = size(vX);
            X = reshape(vX', TfnC/numel(chidx), numel(chidx), nS);
            
        end
        
        %------------------------------------------------------------------
        function X = getWaveform_(self, tfidx, chidx, waveformidx)
            X = self.getData(tfidx,chidx,waveformidx);
        end
        
        %------------------------------------------------------------------
        function gdf = getGdf(self, idx)
            binDims = self.binDims(:,:);
            
            if ~binDims(1)
                gdf = [];
                return
            end
            if nargin < 2 || ( ischar(idx) && strcmp(idx, ':') )
                idx = 1:binDims(1); 
            end
            assert( max(idx) <= binDims(1), 'Index exceeding number of waveforms!') 
            gdf = double(self.gdf(idx,:));
        end
        
        %------------------------------------------------------------------
        function dims = getDims(self)
            %dims = [self.getTf(), self.getNChannels(), self.getNSamples()];
            dims = [self.getNSamples(), self.getTf() * self.getNChannels()];
        end
        
        function N = getNSamples_(self)
            binDims = self.binDims(:,:);
            N = double(binDims(1));
        end
        
        function n = getNChannels(self)
            n = self.nC;
        end
        
        function n = getTf(self)
            n = self.Tf;
        end
        
        function n = getCutLeft(self)
            n = self.cutLeft;
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
            assert( C == self.nC * self.Tf, 'The input waveforms must conform to the filesize!')
            
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
            %disp(['Before: ' num2str(size(self.wfs)) ' - ' num2str(self.binDims(:,:))])
            self.wfs(a0:a1, :) = vWFS;
            %disp(['Middle: ' num2str(size(self.wfs)) ' - ' num2str(self.binDims(:,:))])
            self.binDims(:, :) = int32([a1, self.nC*self.Tf]);
            %disp(['After: ' num2str(size(self.wfs)) ' - ' num2str(self.binDims(:,:))])
            
            self.gdf( a0:a1, :) = gdf;
            assert( isequal(size(self.wfs), self.binDims(:,:)), 'Saved binDims do not correspond to size of wfs!')
            self.nextIdx = self.nextIdx + N;
            
        end
        
        
    end
    
end