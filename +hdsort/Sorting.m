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
                self.setFile(fullfile(self.folders.main, [self.name '_results_discarded.mat']), 'results_discarded');
            else
                self.setFile(P.resultFile, 'results');
                % Todo: Fix this self.setFile(P.resultDiscardedFile, 'results_discarded');
            end
            
            %% Create LEGs and group folders
            self.folders.groups = fullfile(self.folders.main, 'groups');
            [dir_exists,mess,messid] = mkdir( self.folders.groups );
            assert(dir_exists, 'Output directory could not be created!');
            
            self.LEGs = self.createLocalElectrodeGroups(MultiElectrode);
            
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
                {'.040spikes_cut.mat',
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
        
        % -----------------------------------------------------------------
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
        function LEGs = createLocalElectrodeGroups(self, MultiElectrode, varargin)
            P = self.sortingParameters.legs;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
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
                LEGs.P = P;
                
                disp('Create groups...');
                [LEGs.groupsidx, LEGs.nGroupsPerElectrode] = hdsort.leg.constructLocalElectrodeGroups(...
                    LEGs.electrodePositions(:,1), LEGs.electrodePositions(:,2), LEGs.P);
                
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
        function preprocess_summary = preprocess(self, varargin)
            P = self.sortingParameters;
            P.saveRawH5FileNameList = {};
            P = hdsort.util.parseInputs(P, varargin, 'merge');
            
            progress = self.loadProgress();
            
            if isempty(P.saveRawH5FileNameList)
                if isa(self.rawDS{1}, 'hdsort.file.CMOSMEA')
                    try
                        fileList = cellfun(@(x) x.fileNames, self.rawDS, 'UniformOutput', 0);
                    catch
                        fileList = arrayfun(@(x) num2str(x,'%04d.h5'), 1:numel(self.rawDS), 'un',0)
                    end
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
                        try
                            [~, name_, ~] = fileparts(fileList{ii}{1});
                            [~, name,  ~] = fileparts(name_);
                        catch
                            name = sprintf('%04d', ii);
                        end
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
                    self.setFile(preprocess_summary.preprocessedFiles, 'preprocessed');
                else
                    try
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
                self.setFile(progress.preprocessedFiles, 'preprocessed');
                disp('Preprocessing already completed.')
            end
        end
         
        % -----------------------------------------------------------------
        function sort_summary = sort(self, varargin)
            % sortingMode Options: 
            % 'local_parfor' (default) - processes LEGs in a parfor loop
            % 'local'                  - processes LEGs in a for loop
            % 'grid'                   - launches a batch process on a computer grid
            % 'local_parfeval'         - processes LEGs using parfeval (mostly for debugging)
            % 'local_parfor_chunked'   - processes LEGs using a parfor loop but step by step (mostly for debugging)
            
            P.sortingMode = 'local_parfor';
            P.forceExecution = false;
            P.queue = 'regular';
            P.gridType = 'QSUB';
            P.sortjobFolder = self.folders.main;
            P.debug = false;
            P.numWorkers = 8;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            if ~isfield(self.files, 'preprocessed')
                progress = self.loadProgress();
                self.files.preprocessed = progress.preprocessedFiles;
            end
            
            if ~isfield(self.buffer, 'sortjob') || isempty(self.buffer.sortjob)
                
                self.buffer.sortjob = hdsort.grid.SortJob( self.name, P.sortjobFolder, ...
                    self.folders.groups, self.LEGs, self.files.preprocessed, ...
                    self.sortingParameters, 'gridType', P.gridType, 'queue', P.queue)
                
                self.buffer.sortjob.setTaskParameters();
            end
            
            progress = self.loadProgress();
            if P.forceExecution || any(~progress.sortedLEGs)
                
                nCompleted = self.buffer.sortjob.summarySnapshot(true);
                if self.buffer.sortjob.nTasks > nCompleted
                    
                    if strcmp(P.sortingMode, 'grid')
                        self.buffer.sortjob.constructJobSH();
                        self.buffer.sortjob.createAutoSubmitToken();
                    else
                        
                        debugFlag = P.debug;
                        [job_id_str, taskFile] = self.buffer.sortjob.prepareLocalTask();
                        if strcmp(P.sortingMode, 'local')
                            for t = 1:self.buffer.sortjob.nTasks()
                                hdsort.grid.GridJob.runTask(taskFile, num2str(t), job_id_str, debugFlag);
                            end
                            
                        elseif strcmp(P.sortingMode, 'local_parfeval')
                            
                            ppool = gcp();
                            for t = 1:self.buffer.sortjob.nTasks()
                                parfeval(ppool, @hdsort.grid.GridJob.runTask, 0, taskFile, num2str(t), job_id_str, debugFlag);
                            end
                            
                        elseif strcmp(P.sortingMode, 'local_parfor')
                            w = warning('off','all');
                            
                            N = self.buffer.sortjob.nTasks();
                            parfor t = 1:N
                                hdsort.grid.GridJob.runTask(taskFile, num2str(t), job_id_str, debugFlag);
                            end
                            
                        elseif strcmp(P.sortingMode, 'local_parfor_chunked')
                            w = warning('off','all');
                            
                            N = self.buffer.sortjob.nTasks();
                            chunksize = P.numWorkers;
                            C = util.Chunker(N, 'chunkSize', chunksize, ...
                                'minChunkSize', 4, 'progressDisplay', 'console');
                            while C.hasNextChunk()
                                chunks = C.getNextChunk();
                                chunkIdx = chunks(1):chunks(2);
                                
                                N_chunk = numel(chunkIdx);
                                parfor t = 1:N_chunk
                                    task_str = num2str( chunkIdx(t) );
                                    hdsort.grid.GridJob.runTask(taskFile, task_str, job_id_str, debugFlag);
                                end
                            end
                        end
                    end
                    
                else
                    disp('All tasks seem to be completed...')
                end
                all_tasks_completed = self.buffer.sortjob.waitForTasksToFinish(60);
                assert(all_tasks_completed, 'Error in at least a few tasks!')
            end
            
            % Save summary:
            sort_summary = self.buffer.sortjob.summarizeReports();
            self.saveProgress('sortedLEGs', sort_summary.completedTasks);
            self.saveProgress('sort_summary', sort_summary);
            
        end
        
        % -----------------------------------------------------------------
        function sort_summary = reuptake(self, varargin)
            
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
        function [SortedPopulation, SortedPopulation_discarded] = createSortedPopulation(self, outputLocation)
            
            if isfield(self.buffer, 'SortedPopulation') && ~isempty(self.buffer.SortedPopulation) && ~isempty(self.buffer.SortedPopulation_discarded)
                SortedPopulation = self.buffer.SortedPopulation;
                SortedPopulation_discarded = self.buffer.SortedPopulation_discarded;
            else
                
                try
                    %assert(self.sortingResultFileExists(outputLocation), 'Create file...');
                    assert( exist(self.files.results) > 0, 'Create results file...')
                    disp('Loading SortedPopulation file...')
                    load(self.files.results);
                    disp(['SortedPopulation ' SortedPopulation.name ' file loaded.'])
                catch
                    disp('Create SortedPopulation Population structure...')
                    
                    if ~isfield(self.buffer, 'rawResults') || isempty(self.buffer.rawResults)
                        self.buffer.rawResults = load(self.files.rawResults);
                    end
                    
                    preprocessor = load(self.files.preprocessor)
                    noiseStd = preprocessor.smadPerEl;
                    
                    [SortedPopulation, SortedPopulation_discarded] = hdsort.results.createPopulation(...
                        self.name, self.buffer.rawResults, self.rawDS, noiseStd);
                    
                    SortedPopulation.fileLocation = outputLocation;
                    SortedPopulation_discarded.fileLocation = outputLocation;
                    
                    %save(self.files.results, 'SortedPopulation', 'SortedPopulation_discarded', '-v7.3');
                    SortedPopulation.save(self.files.results)
                    SortedPopulation_discarded.save(self.files.results_discarded);
                    disp('New hdsort.results.Population saved!')
                    
                end
                
                self.buffer.SortedPopulation = SortedPopulation;
                self.buffer.SortedPopulation_discarded = SortedPopulation_discarded;
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
                spikes_cut = hdsort.file.WaveFormFileMat(fname);
                x.spikes_cut = spikes_cut;
            else
                x = load(fname);
            end
            
            out = hdsort.util.mergeStructs(out, x);
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
        
    end
end