classdef Sorting < handle
    properties
        
        rawDS
        LEGs
        name
        
        folders
        files
        buffer
        
        sortingParameters
        
        groupFilesList
        all_tasks_complete
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Sorting(rawDSList, mainFolder, sortingName, varargin)
            P.maxElPerGroup = 9;
            P.legsFile = '';
            P.rawResultFile = '';
            P.resultFile = '';
            P.preprocessFile = '';
            P.sortingParameters = hdsort.defaultParameters();
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            %% Initialize sortingParameters
            self.sortingParameters = P.sortingParameters;
            
            %% Initialize the raw DS:
            if ~iscell(rawDSList)
                rawDSList = {rawDSList};
            end
            self.rawDS = rawDSList;
            MultiElectrode = self.rawDS{1}.MultiElectrode;
            for ii = 2:numel(self.rawDS)
                assert(isequal(MultiElectrode, self.rawDS{ii}.MultiElectrode), 'All raw DS must have the same MultiElectrode!')
            end
            
            %% Initialize file names:
            self.name = sortingName;
            self.folders.main = fullfile(mainFolder, self.name);
            [dir_exists,mess,messid] = mkdir( self.folders.main );
            assert(dir_exists, 'Main directory could not be created!');
            
            %self.files.sortingclass = fullfile(self.folders.main, [ self.name '_sorting.mat']);
            self.setFile(fullfile(self.folders.main, [ self.name '_sorting.mat']), 'sortingclass');
            
            if ~exist(self.files.sortingclass)
                save(self.files.sortingclass, 'sortingName');
            end
            
            if isempty(P.legsFile)
                self.setFile(fullfile(self.folders.main, 'legsFile.mat'), 'legs');
                %self.files.legs = fullfile(self.folders.main, 'legsFile.mat');
            else
                self.setFile(P.legsFile, 'legs');
                %self.files.legs = P.legsFile;
            end
            
            if isempty(P.preprocessFile)
                self.setFile(fullfile(self.folders.main, [self.name '_preprocessor.mat']), 'preprocessor');
                %self.files.preprocessor = fullfile(self.folders.main, [self.name '_preprocessor.mat']);
            else
                self.setFile(P.preprocessFile, 'preprocessor');
                %self.files.preprocessor = P.preprocessFile;
            end
            
            if isempty(P.rawResultFile)
                self.setFile(fullfile(self.folders.main, [self.name '_rawres.mat']), 'rawResults');
                %self.files.rawResults = fullfile(self.folders.main, [self.name '_rawres.mat']);
            else
                self.setFile(P.rawResultFile, 'rawResults');
                %self.files.rawResults = P.rawResultFile;
            end
            
            if isempty(P.resultFile)
                self.setFile(fullfile(self.folders.main, [self.name '_results.mat']), 'results');
                %self.files.results = fullfile(self.folders.main, [self.name '_results.mat']);
            else
                self.setFile(P.resultFile, 'results');
                %self.files.results = P.resultFile;
            end
            
            %% Create LEGs and group folders
            self.folders.groups = fullfile(self.folders.main, 'groups');
            [dir_exists,mess,messid] = mkdir( self.folders.groups );
            assert(dir_exists, 'Output directory could not be created!');
            
            self.LEGs = self.createLocalElectrodeGroups(MultiElectrode, P.maxElPerGroup);
            
            for ii = 1:self.LEGs.N
                self.folders.legs{ii} = fullfile(self.folders.groups, self.LEGs.name{ii});
                [dir_exists, mess, messid] = mkdir( self.folders.legs{ii}  );
                assert(dir_exists, 'Output directory could not be created!');
            end
            
            %% Keep track of all files in each group folder:
            keySet = ...
                {'spikes_cut',
                'spikes_aligned',
                'cov',
                'spikes_prewhitened',
                'spikes_features',
                'clusters_meanshift',
                'botm_matching',
                'clusters_meanshift_merged',
                'P',
                'gdf',
                'templates'};
            valueSet = ...
                {'.040spikes_cut.h5',
                '.050spikes_aligned.mat',
                '.060cov.mat',
                '.070spikes_prewhitened.mat',
                '.080spikes_features.mat',
                '.090clusters_meanshift.mat',
                '.100botm_matching.mat',
                '.110clusters_meanshift_merged.mat',
                '.P.mat',
                '.gdf.mat',
                '_templates.mat'};
            
            self.groupFilesList = containers.Map(keySet,valueSet);
            
            self.buffer.progress = self.loadProgress();
        end
        
        function display(self)
            disp(self)
            disp('##############################################')
            disp(['Sorting Name: ' self.name])
            
            disp('----------- Sorting file names ---------------')
            disp(self.groupFilesList.keys')
            
            if ~isempty(self.LEGs)
                disp(['Number of LEGs: ' num2str(self.LEGs.N)])
                % Add here some info on the LEGs
            end
            
            if isfield(self.buffer, 'progress')
                if self.buffer.progress.postprocessed
                    disp('Post-processing finished.')
                elseif sum(self.buffer.progress.sortedLEGs)
                    disp([num2str(sum(self.buffer.progress.sortedLEGs)) ' LEGs have been sorted.'])
                elseif self.buffer.progress.preprocessed
                    disp('Pre-processing finished.')
                end
            end
            
            if ~isfield(self.buffer, 'result')
                disp(' ')
                if exist(self.files.results)
                    disp('----------------------------------------------')
                    disp('Result file already exists');
                else
                    disp('-----------------------')
                    disp('Result file not yet created');
                end
            else
                disp(' ')
                disp('----------------------------------------------')
                disp('Results already loaded:');
            end
            
        end
        
        % -----------------------------------------------------------------
        function LEGs = createLocalElectrodeGroups(self, MultiElectrode, maxElPerGroup)
            try
                LEGs = load(self.files.legs)
                if ~isfield(LEGs, 'N') || ~isfield(LEGs, 'name')
                    LEGs.N = numel(LEGs.groups);
                    for ii= 1:LEGs.N
                        LEGs.name{ii} = sprintf('group%04d', ii);
                    end
                end
                disp('Groups already exist')
            catch
                LEGs.electrodeNumbers = MultiElectrode.electrodeNumbers;
                LEGs.electrodePositions = MultiElectrode.electrodePositions;
                LEGs.maxElPerGroup = maxElPerGroup;
                
                disp('Create groups...');
                [LEGs.groupsidx, LEGs.nGroupsPerElectrode] = hdsort.leg.constructLocalElectrodeGroups(...
                    LEGs.electrodePositions(:,1), LEGs.electrodePositions(:,2), ...
                    'maxElPerGroup', maxElPerGroup);
                disp(['Number of groups created: ' num2str(length(LEGs.groupsidx))])
                
                LEGs.groups = {};
                for ii= 1:length(LEGs.groupsidx)
                    LEGs.groups{ii} = LEGs.electrodeNumbers(LEGs.groupsidx{ii});
                    LEGs.name{ii} = sprintf('group%04d', ii);
                end
                
                LEGs.N = numel(LEGs.groups);
                save(self.files.legs, '-struct', 'LEGs')
            end
        end
        
        % -----------------------------------------------------------------
        function R = loadRawResults(self)
            if isempty(self.buffer.rawResults)
                try
                    load(self.files.rawResults);
                    self.buffer.rawResults = R;
                catch
                    error('Could not load sorting result, run sorter first!');
                end
            end
            R = self.buffer.rawResults;
        end
        
        % -----------------------------------------------------------------
        function progress = loadProgress(self)
            try
                load(self.files.sortingclass, 'progress');
                assert( isfield(progress, 'preprocessed'), '!')
            catch
                progress.preprocessed = false;
                progress.preprocess_summary = [];
                progress.preprocessedFiles = {};
                
                progress.sortedLEGs = false(self.LEGs.N, 1);
                progress.sort_summary = [];
                
                progress.postprocessed = false;
                progress.postprocess_summary = [];
                
                save(self.files.sortingclass, 'progress', '-append');
            end
            self.buffer.progress = progress;
        end
        function progress = saveProgress(self, key, value)
            progress = self.loadProgress();
            progress.(key) = value;
            save(self.files.sortingclass, 'progress', '-append');
        end

        % -----------------------------------------------------------------
        function preprocess_summary = preprocess(self, varargin)
            P = self.sortingParameters;
            P.saveRawH5FileNameList = {};
            P = hdsort.util.parseInputs(P, varargin, 'merge');
            
            progress = self.loadProgress();
            
            if isempty(P.saveRawH5FileNameList)
                if isa(self.rawDS{1}, 'hdsort.file.CMOSMEA')
                    fileList = self.rawDS{1}.fileNames;
                else
                    try
                        fileList = cellfun(@(x) x.fileName, self.rawDS, 'UniformOutput', 0);
                    catch
                        fileList = arrayfun(@(x) num2str(x,'%04d.h5'), 1:numel(self.rawDS), 'un',0)
                    end
                end
                
                for ii = 1:numel(fileList)
                    try
                        [~, name_, ~] = fileparts(fileList{ii});
                        [~, name,  ~] = fileparts(name_);
                    catch
                        name = sprintf('%04d', ii);
                    end
                    P.saveRawH5FileNameList{ii} = [name, '.h5'];
                end
            end
            
            if ~progress.preprocessed
                self.buffer.preprocessor = hdsort.Preprocessor( ...
                        self.rawDS, self.LEGs.groupsidx, self.name, P);
                    
                disp('Preprocessing...')
                forceExecution = false;
                preprocess_summary = self.buffer.preprocessor.preprocess(...
                    self.folders.groups, self.files.preprocessor, forceExecution);
                
                
                if ~isempty(preprocess_summary.preprocessedFiles)
                    %self.files.preprocessed = hdsort.util.convertPathToOS(preprocess_summary.preprocessedFiles);
                    self.setFile(preprocess_summary.preprocessedFiles, 'preprocessed');
                else
                    try
                        %self.files.preprocessed = hdsort.util.convertPathToOS(self.rawDS{1}.fileNames);
                        self.setFile(self.rawDS{1}.fileNames, 'preprocessed');
                    catch
                        self.files.preprocessed = {};
                    end
                end
                
                % Save summary:
                self.saveProgress('preprocessed', true);
                self.saveProgress('preprocessedFiles', self.files.preprocessed);
                self.saveProgress('preprocess_summary', preprocess_summary);
                
                disp('Preprocessing finished.')
            else
                preprocess_summary = progress.preprocess_summary;
                %self.files.preprocessed = hdsort.util.convertPathToOS(progress.preprocessedFiles);
                self.setFile(progress.preprocessedFiles, 'preprocessed');
                disp('Preprocessing already completed.')
            end
        end
        
        function file_exists = setFile(self, fileName, name)
            self.files.(name) = hdsort.util.convertPathToOS(fileName);
            
            if ~iscell(self.files.(name))
                file_exists = exist(self.files.(name));
            else
                file_exists = false;
                for f = self.files.(name)
                    if ~exist(f{1}) return; end
                end
                file_exists = true;
            end
        end
        
        % -----------------------------------------------------------------
        function sort_summary = sort(self, varargin)
            
            P.sortingMode = 'grid'; % as opposed to 'local' , 'local_parfor'
            P.forceExecution = false;
            P.queue = 'regular';
            P.gridType = 'QSUB';
            P.sortjobFolder = self.folders.main;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            if ~isfield(self.buffer, 'sortjob') || isempty(self.buffer.sortjob)
            
                progress = self.loadProgress();
                preprocessedFiles = progress.preprocessedFiles;
                
                self.buffer.sortjob = hdsort.grid.SortJob( self.name, P.sortjobFolder, ...
                    self.folders.groups, self.LEGs, preprocessedFiles, self.sortingParameters, ...
                    'gridType', P.gridType, 'queue', P.queue)
                
                self.buffer.sortjob.setTaskParameters();
            end
            
            progress = self.loadProgress();
            if P.forceExecution || any(~progress.sortedLEGs)
                
                if strcmp(P.sortingMode, 'grid')
                    
                    self.buffer.sortjob.constructJobSH();
                    
                    [nCompleted, tasksNotCompleted, nErrors, tasksWithErrors] = self.buffer.sortjob.summarySnapshot(true);
                    if self.buffer.sortjob.nTasks > nCompleted
                        self.buffer.sortjob.createAutoSubmitToken();
                    else
                        disp('All tasks seem to be completed...')
                    end
                    all_tasks_completed = self.buffer.sortjob.waitForTasksToFinish(60);
                    assert(all_tasks_completed, 'Error in at least a few tasks!')
                    
                elseif strcmp(P.sortingMode, 'local')
                    
                    for t = 1:self.buffer.sortjob.nTasks()
                        self.buffer.sortjob.runTaskLocally(t);
                    end
                    
                    % This could be used to avoid hdsort.grid.SortJob, but
                    % currently, we want to use the summarizeReports
                    % function.
%                     for t = 1:self.buffer.sortjob.nTasks()
%                         groupFolder = self.folders.legs{t};
%                         wfsFile = fullfile(groupFolder, [taskParameters.name '.040spikes_cut.h5']);
%                         covFile = fullfile(groupFolder, [taskParameters.name '.060cov.mat']);
%                         groupidx = self.LEGs.groupsidx{t};
%                         
%                         hdsort.leg.sortOneLEG(self.name, wfsFile, covFile, ...
%                             groupFolder, self.sortingParameters, groupidx, ...
%                             self.files.preprocessed);
%                     end
                    
                elseif strcmp(P.sortingMode, 'local_parfor')
                    
                    parfor t = 1:self.buffer.sortjob.nTasks()
                        self.buffer.sortjob.runTaskLocally(t);
                    end
                    
                else
                    error('sortingMode unknown - use either grid, local or local_parfor!')
                end
            end
            
            %% Save summary:
            sort_summary = self.buffer.sortjob.summarizeReports();
            self.saveProgress('sortedLEGs', sort_summary.completedTasks);
            self.saveProgress('sort_summary', sort_summary);
        end
        
        
        % -----------------------------------------------------------------
        function sort_summary = reuptake(self, varargin)
            
            P.sortingMode = 'grid'; % as opposed to 'local' , 'local_parfor'
            P.queue = 'regular';
            P.gridType = 'QSUB';
            P.sortjobFolder = self.folders.main;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            if ~isfield(self.buffer, 'sortjob')
                
                progress = self.loadProgress();
                preprocessedFiles = progress.preprocessedFiles;
                
                self.buffer.sortjob = hdsort.grid.SortJob( self.name, P.sortjobFolder, ...
                    self.folders.groups, self.LEGs, preprocessedFiles, self.sortingParameters, ...
                    'gridType', P.gridType, 'queue', P.queue)
                
                self.buffer.sortjob.setTaskParameters();
            end
            
            all_tasks_completed = self.buffer.sortjob.waitForTasksToFinish(60);
            assert(all_tasks_completed, 'Error in at least a few tasks!')
            
            %% Save summary:
            sort_summary = self.buffer.sortjob.summarizeReports();
            self.saveProgress('sortedLEGs', sort_summary.completedTasks);
            self.saveProgress('sort_summary', sort_summary);
        end
        
        % -----------------------------------------------------------------
        function rawResults = postprocess(self, varargin)
            
            P.forceExecution = false;
            P = hdsort.util.parseInputs(P, varargin, 'merge');
            
            if exist(self.files.rawResults) && ~P.forceExecution
                rawResults = load(self.files.rawResults);
                self.saveProgress('postprocess', true);
                return;
            end
            
            %% Run Postprocessor:
            disp('Start Postprocessor...')
            tic
            self.buffer.postprocessor = hdsort.Postprocessor(...
                self.LEGs, self.folders.legs, self.name, self.folders.groups, ...
                'gdfFileNames', self.groupFilesList('gdf'),...
                'templateFileName', self.groupFilesList('templates'), ...
                'covFileName', self.groupFilesList('cov') )
            
            rawResults = self.buffer.postprocessor.run();
            postprocess_summary = self.buffer.postprocessor.summarizeReports();
            
            %% Check output:
            disp('Check and save data...');
            units = unique(rawResults.gdf_merged(:,1));
            nU = length(units);
            assert(length(rawResults.localSorting) == nU, 'must be identical');
            assert(length(rawResults.localSortingID) == nU, 'must be identical');
            assert(size(rawResults.T_merged,3) == nU, 'must be identical');
            
            %% Save summary:
            self.saveProgress('postprocess_summary', postprocess_summary);
            self.saveProgress('postprocessed', true);
            
            %% Save rawResults:
            disp('Saving postprocessing rawResults...')
            progress = self.loadProgress();
            
            rawResults.summary.preprocess  = progress.preprocess_summary;
            rawResults.summary.sort        = progress.sort_summary;
            rawResults.summary.postprocess = progress.postprocess_summary;
            
            save(self.files.rawResults, '-struct', 'rawResults', '-v7.3');
            disp('Done.')
            
            self.buffer.rawResults = rawResults;
        end

        % -----------------------------------------------------------------
        %function fileName = sortingResultFileName(self, outputLocation)
        %    fileName = fullfile(outputLocation, [ self.name '_results.mat']);
        %end
        %function srf_exists = sortingResultFileExists(self, outputLocation)
        %    sortingResultFileName = self.sortingResultFileName(outputLocation)
        %    srf_exists = exist(sortingResultFileName) > 0;
        %end
        
        % -----------------------------------------------------------------
        function [SortingResults, SortingResults_discarded] = createSpikeSortingResult(self, outputLocation)
            
            if isfield(self.buffer, 'SortingResults') && ~isempty(self.buffer.SortingResults) && ~isempty(self.buffer.SortingResults_discarded)
                SortingResults = self.buffer.SortingResults;
                SortingResults_discarded = self.buffer.SortingResults_discarded;
            else
                
                try
                    %assert(self.sortingResultFileExists(outputLocation), 'Create file...');
                    assert( exist(self.files.results) > 0, 'Create results file...')
                    disp('Loading SortingResults file...')
                    load(self.files.results);
                    disp(['SortingResults ' SortingResults.name ' file loaded.'])
                catch
                    disp('Create SortingResults SpikeSorting structure...')
                    
                    if ~isfield(self.buffer, 'rawResults') || isempty(self.buffer.rawResults)
                        self.buffer.rawResults = load(self.files.rawResults);
                    end
                    
                    preprocessor = load(self.files.preprocessor)
                    noiseStd = preprocessor.smadPerEl;
                    
                    [SortingResults, SortingResults_discarded] = hdsort.results.createSpikeSorting(...
                        self.name, self.buffer.rawResults, self.rawDS, noiseStd);
                    
                    SortingResults.filePath = outputLocation;
                    SortingResults_discarded.filePath = outputLocation;
                    
                    save(self.files.results, 'SortingResults', 'SortingResults_discarded', '-v7.3');
                    disp('New SortingResults SpikeSorting saved!')
                end
                
                self.buffer.SortingResults = SortingResults;
                self.buffer.SortingResults_discarded = SortingResults_discarded;
            end
            
        end
        
        
        
        % -----------------------------------------------------------------
        % ---------------------- DEBUGGING FUNCTIONS ----------------------
        % -----------------------------------------------------------------
        %% GROUP ORGANISATION
        function keys = listFilesInGroup(self)
            keys = self.groupFilesList.keys;
            disp(keys)
        end
        
        function [fname, does_exist] = getFileNameInGroup(self, legNr, key)
            foldername = self.folders.legs{legNr};
            fname = fullfile(foldername, [self.name self.groupFilesList(key)]);
            does_exist = exist(fname);
        end
        
        function out = loadFileInGroup(self, legNr, key, out)
            if nargin < 4
                out = struct();
            end
            fname = self.getFileNameInGroup(legNr, key);
            
            if strcmp(key, 'spikes_cut')
                spikes_cut = hdsort.file.WaveFormFile(fname);
                x.spikes_cut = spikes_cut;
            else
                x = load(fname);
            end
            
            out = util.mergeStructs(out, x);
        end
        
        % -----------------------------------------------------------------
        function out = deleteFilesInGroups(self, key)
            out = true;
            for legNr = 1:self.LEGs.N
                fname = self.getFileNameInGroup(legNr, key);
                
                try
                    delete(fname);
                catch
                    out = false;
                end
            end
        end
        
        % -----------------------------------------------------------------
        function G = loadGroupStruct(self)
            %self.files.GroupStruct   = fullfile(self.folders.groups, 'G_struct.mat');
            self.setFile(fullfile(self.folders.groups, 'G_struct.mat'), 'GroupStruct');
            if isempty(self.buffer.GroupStruct)
                self.buffer.GroupStruct = load(self.files.GroupStruct);
            end
            G = self.buffer.GroupStruct;
        end
        
        % -----------------------------------------------------------------
        function gdfs = loadGroupGDFs(self)
            G = self.loadGroupStruct();
            gdfs = {G.G.gdf};
        end
        
        % -----------------------------------------------------------------
        % ------------------------- PLOT FUNCTIONS ------------------------
        % -----------------------------------------------------------------
        function P = plotLEGs(self, varargin)
            P.fh = [];
            P.ah = [];
            P = util.parseInputs(P, varargin, 'error');
            
            [col, markerSet] = bel_plots.PlotInterface.vectorColor(1:self.LEGs.N);
            for ii = 1:self.LEGs.N
                ep = self.LEGs.electrodePositions(self.LEGs.groupsidx{ii},:);
                P = bel_plots.Gscatter(ep(:,1), ep(:,2), [], 'ah', P.ah, 'fh', P.fh, 'color', col(ii,:), 'Marker', markerSet{ii});%, 'MarkerSize', 16);
            end
        end
        
    end
end