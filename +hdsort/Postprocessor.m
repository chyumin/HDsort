classdef Postprocessor < handle
    % This class gets the sorting results of different LEGs an input,
    % removes units that did not have their template maximum in the LEG in
    % which they were found, finds and resolves duplicates, and then saves
    % and returns the final results.
    
    properties
        P
        LEGs
        name
        
        files
        folders
        
        buffer
    
        summary
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Postprocessor(LEGs, legFolders, sortingName, outputLocation, varargin)
            
            P.debug = false;
            P.useparfor = false;
            
            P.gdfFileNames = '.gdf.mat';
            P.templateFileName = '_templates.mat';
            P.covFileName = '.060cov.mat';
            P.parametersFileName = '.P.mat';
            
            self.P = hdsort.util.parseInputs(P, varargin);
            
            assert(LEGs.N == numel(legFolders), 'You must provide as many legFolders as there are LEGs!')
            
            self.LEGs = LEGs;
            self.folders.legs = legFolders;
            self.name = sortingName;
            
           % self.folders.
            self.files.GStruct = fullfile(outputLocation, 'G_struct.mat');
            self.summary.durations = zeros(3,1);
        end
        
        function results = run(self)
            self.createGStruct();
            self.merge();
            results = self.getFinalResults();
        end
        
        
        % -----------------------------------------------------------------
        function [G, meanNoiseStd] = createGStruct(self)
            
            try
                self.buffer = load(self.files.GStruct);
                disp('G structure already created.')
                return;
            catch
                disp('Create new GStruct...')
            end
            
            t1 = tic;
            
            G = struct('sortedElectrodes', {}, 'gdf', {}, 'units', {});
            nUnitsPerLocalSorting = zeros(1, self.LEGs.N);
            Ns = cell(self.LEGs.N,1);
            nT = 0;
            
            for ii = 1:self.LEGs.N
                
                %% Reshape the GDFs:
                load( fullfile( self.folders.legs{ii}, [self.name, self.P.gdfFileNames]) );
                units = unique(gdf(:,1));
                G(ii).sortedElectrodes = self.LEGs.groupsidx{ii};
                G(ii).gdf = gdf;
                G(ii).units = units;
                nUnitsPerLocalSorting(ii) = length(units);
                
                %% Load Templates
                G(ii).templates = load( fullfile( self.folders.legs{ii}, [self.name, self.P.templateFileName]) );
                nT = nT + size(G(ii).templates.wfs,3);
            
                %% Load Noise
                %S = load( fullfile( self.folders.legs{ii}, [self.name, self.P.covFileName]) );
                %Ns{ii} = S.noise.CestS.CCol;
                
                %% Load Parameters
                parameters = load( fullfile( self.folders.legs{ii}, [self.name, self.P.parametersFileName]) );
                G(ii).P = parameters.S.P;
                
                if ii > 1
                    assert(isequal(G(ii).P, G(1).P), 'The different LEGs were sorted with different parameters!')
                end
            end
            
            disp(['Total number of units found: ' num2str(sum(nUnitsPerLocalSorting)) ])
            disp(['Total number of templates found: ' num2str(nT) ])
            
            buffer.G = G;
            buffer.is_merged = false;
            %meanNoiseStd = mean(sqrt(diag(Ns{1})));
            %buffer.meanNoiseStd = meanNoiseStd;
            self.summary.durations(1) = toc(t1);
            
            save(self.files.GStruct, '-struct', 'buffer', '-v7.3');
            
            self.buffer = buffer;
        end
        
        % -----------------------------------------------------------------
        function G = merge(self)
            
            if ~isempty(self.buffer.G) && self.buffer.is_merged
                G = self.buffer.G;
                return;
            end
            
            assert( exist(self.files.GStruct) == 2, 'GStruct does not exist yet!')
            
            try
                self.buffer = load(self.files.GStruct);
                assert( self.buffer.is_merged , '!');
                disp('G structure already merged.')
                return;
            catch
                disp('Merge...')
            end
            
            t2 = tic;
            
            %% Check the G structure:
            for ii = 1:length(self.buffer.G)
                nUnitsInGDF = length(unique(self.buffer.G(ii).gdf(:,1)));
                
                if ~isempty(self.buffer.G(ii).templates.wfs)
                    nUnitsInTempates = size(self.buffer.G(ii).templates.wfs,3);
                else
                    nUnitsInTempates = 0;
                end
                str_ = sprintf('i: %d, nGdf: %d, nWfs: %d, rn: %s', ii, nUnitsInGDF, nUnitsInTempates, self.name);
                assert(nUnitsInGDF == nUnitsInTempates, ['must be identical - ' str_])
            end
            
            %% Merge double Templates
            fprintf('Merging double Templates...');
            self.buffer.G = hdsort.leg.mergeLocalSortings(self.buffer.G);
            
            for ii = 1:length(self.buffer.G)
                if ~isempty(self.buffer.G(ii).gdf(:,1))
                    assert(length(unique(self.buffer.G(ii).gdf(:,1))) == size(self.buffer.G(ii).templates.wfs,3), ['must be identical - ' str_])
                else
                    assert(isempty(self.buffer.G(ii).templates.wfs), 'If there is no spike, there cannot be a waveform!')
                end
            end
            
            self.buffer.is_merged = true;
            
            self.summary.durations(2) = toc(t2);
            
            buffer = self.buffer;
            save(self.files.GStruct, '-struct', 'buffer', '-v7.3');
            disp('Done.')
        end
        
        % -----------------------------------------------------------------
        function results = getFinalResults(self)
            
            if isfield(self.buffer, 'results') && ~isempty(self.buffer.results) 
                results = self.buffer.results;
                return;
            end
            
            if isempty(self.buffer.G)
                assert( exist(self.files.GStruct) == 2, 'GStruct does not exist yet!')
                self.buffer = load(self.files.GStruct);
            end
            assert( self.buffer.is_merged , '!');
            
            t3 = tic;
            
            %% Extract final gdf and templates
            disp('Extracting final gdf and templates...');
            results = hdsort.leg.computeFinalGDFFromMergedLocalSortings(self.buffer.G);
            nT = size(results.T_merged,3);
            disp(['Number of templates after merging: ' num2str(nT) ])
            
            results.gdf_merged = sortrows(results.gdf_merged,2);
            if ~isempty(results.gdf_discarded)
                results.gdf_discarded = sortrows(results.gdf_discarded, 2);
            end
            results.G = self.buffer.G;
            results.sortingParameters = self.buffer.G(1).P;
            
            self.summary.durations(3) = toc(t3);
            
            self.buffer.results = results;
            disp('Final results computed.')
        end
        
        % -----------------------------------------------------------------
        function summary = summarizeReports(self)
            self.summary.totalDuration = sum(self.summary.durations);
            summary = self.summary; 
        end
        
        
    end
end

