%Example:
% sj = grid.SortJob(jobName, dataFiles, destFolder);
% sj.createBOTMGroups();
% sj.copyDataFiles();
% sj.setTaskParameters();
% sj.prepareTasks();
%
% --> Submit to grid
%
% all_tasks_completed = sj.waitForTasksToFinish();
% if all_tasks_completed
%     sj.copyBackResults();
%     sj.postprocessBOTMSorting();
%     sj.runQC()
% end

classdef SortJob < grid.GridJob
    properties (SetAccess=private)
    end
    
    properties
        sortJobP
        groupsidx
        MES
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = SortJob(jobName, rootFolder, dataFiles, varargin)
            self = self@grid.GridJob(['sort_' jobName], rootFolder, varargin{:});
            self.taskType = 'SortJob';
            
            p = struct;
            p.useFilter = 0;
            p.nGroups = [];
            p.groupFile = '';
            
            p = hdsort.util.parseInputs(p, self.P_untreated, 'error');
            self.sortJobP = p;
            
            %% Set the data files:
            if ~iscell(dataFiles) dataFiles = {dataFiles}; end
            self.files.data = dataFiles;
            %self.folders.root = rootPath;
            
            try
                load(self.sortJobP.groupFile)
                self.groupsidx = groupsidx;
                self.startIndex = 1;
                self.endIndex = length(self.groupsidx);
                self.taskIDs = self.startIndex:self.endIndex;
                
                %% Create SortJob specific output folder (on scratch):
                self.folders.groups = fullfile(self.folders.main, 'groups');
                %[dir_exists,mess,messid] = mkdir( self.folders.groups );
                assert(dir_exists, 'Output directory could not be created!');
                
            catch
                self.createBOTMGroups();    
            end
            
        end
        
        % -----------------------------------------------------------------
        function createBOTMGroups(self)
            warning('grid.SortJob.createBOTMGroups() deprecated!')
            
            disp('Create DataStructure object...');
            DSFull = hdsort.filewrapper.CMOSMEA(self.files.data, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);%'PREFILT');
            self.MES = DSFull.MultiElectrode.toStruct();
            
            %% Create groups based on the electrode positions and save them to a file:
            self.files.groupFile = fullfile(self.folders.main, 'groupFile.mat');
            %if ~isempty(self.P.dataPath)
            %    self.destinationlocation.files.groupFile = self.files.groupFile;
            %end
            
            try
                load(self.files.groupFile)
                disp('Groups already exist')
            catch
                disp('Create groups...');
                electrodePositions = DSFull.MultiElectrode.electrodePositions;
                electrodeNumbers   = DSFull.MultiElectrode.electrodeNumbers;
                [groupsidx nGroupsPerElectrode] = hdsort.leg.constructLocalElectrodeGroups(electrodePositions(:,1), electrodePositions(:,2), 'maxElPerGroup', 9);
                
                %% Limit number of groups if necessary:
                if ~isempty(self.sortJobP.nGroups)
                    groupsidx = {groupsidx{1:self.sortJobP.nGroups}};
                end
                
                disp(['Number of groups created: ' num2str(length(groupsidx))])
                
                groups = {};
                for ii= 1:length(groupsidx)
                    groups{ii} = electrodeNumbers(groupsidx{ii});
                end
                
                save(self.files.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
            end
            
            self.groupsidx = groupsidx;
            self.startIndex = 1;
            self.endIndex = length(self.groupsidx);
            self.taskIDs = self.startIndex:self.endIndex;
            
            
            %% Create SortJob specific output folder (on scratch):
            self.folders.groups = fullfile(self.folders.main, 'groups');
            [dir_exists,mess,messid] = mkdir( self.folders.groups );
            assert(dir_exists, 'Output directory could not be created!');
            
            
        end
        
        function setTaskParameters(self)
            taskParameters = struct;
            
            %% Sorting parameters form ana.startHDSorting:
            taskParameters.sortingParameters = struct;
            taskParameters.sortingParameters.spikeDetection.method = '-';
            taskParameters.sortingParameters.spikeDetection.thr = 4.2;
            taskParameters.sortingParameters.artefactDetection.use = 0;
            taskParameters.sortingParameters.botm.run = 0;
            taskParameters.sortingParameters.spikeCutting.maxSpikes = 200000000000; % Set this to basically inf
            
            taskParameters.sortingParameters.noiseEstimation.minDistFromSpikes = 80;
            
            taskParameters.sortingParameters.spikeAlignment.initAlignment = '-';
            taskParameters.sortingParameters.spikeAlignment.maxSpikes = 50000;     % so many spikes will be clustered
            taskParameters.sortingParameters.clustering.maxSpikes = taskParameters.sortingParameters.spikeAlignment.maxSpikes;  % dont align spikes you dont cluster...
            taskParameters.sortingParameters.clustering.meanShiftBandWidth = sqrt(1.8*6);
            taskParameters.sortingParameters.mergeTemplates.merge = 1;
            taskParameters.sortingParameters.mergeTemplates.upsampleFactor = 3;
            taskParameters.sortingParameters.mergeTemplates.atCorrelation = .93; % DONT SET THIS TOO LOW! USE OTHER ELECTRODES ON FULL FOOTPRINT TO MERGE
            taskParameters.sortingParameters.mergeTemplates.ifMaxRelDistSmallerPercent = 30;
            %% ---------------------------------------------
            
            taskParameters.dataFiles = self.files.data;
            taskParameters.runName = self.jobName;
            taskParameters.outputPath = self.folders.groups;
            
            if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                %warning('Sorting not started from a linux machine might cause problems!')
                taskParameters.outputPath = grid.GridJob.convertToLinux(taskParameters.outputPath);
                taskParameters.dataFiles = grid.GridJob.convertToLinux(taskParameters.dataFiles);
            end
            
            disp('Create task files from groups...');
            %% Create cell variable allTaskParamters:
            for ii = 1:length(self.taskIDs)
                
                self.files.groupPaths{ii} = ['group' sprintf('%04d', self.taskIDs(ii))];
                
                taskParameters.groupidx = self.groupsidx{ii};
                taskParameters.taskID = self.taskIDs(ii);
                taskParameters.groupPath = self.files.groupPaths{ii};
                
                self.allTaskParameters{ii} = taskParameters;
            end
            
            self.prepareTasks();
        end
        
        % -----------------------------------------------------------------
        function postprocessBOTMSorting(self, newPostProcFunc)
            warning('This function is deprecated and should not be used anymore. Use hdsort.Sorting()')
            
            % Process all local sortings into a final sorting.
            
            % newPostProcFunc is a string that specifies which
            % prostprocessing function should be used (merging of
            % templates).
            % Default is '', i.e. no string is passed on.
            if nargin < 2
                newPostProcFunc = '';
            end
            
            disp('Start postprocessBOTMSorting...')
            self.destinationlocation.files.results = fullfile(self.folders.main, [self.jobName newPostProcFunc '_results.mat']);
            try
                load(self.destinationlocation.files.results)
                disp('Postprocessing already finished.')
            catch
                disp('Load group information...');
                GF = load(self.destinationlocation.files.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
                
                disp('Start postprocessing...');
                %                 [gdf_merged, T_merged, localSorting, localSortingID, G] =
                [R, P] = hdsort.leg.processLocalSortings(...
                    self.folders.groups,...
                    self.jobName, GF.groups, GF.groupsidx, ...
                    'groupPaths', self.destinationlocation.files.groupPaths, ...
                    'newPostProcFunc', newPostProcFunc);
                gdf_merged     = R.gdf_merged;
                gdf_discarded  = R.gdf_discarded;
                T_merged       = R.T_merged;
                localSorting   = R.localSorting;
                localSortingID = R.localSortingID;
                G              = R.G;
                disp('Check and save data...');
                units = unique(gdf_merged(:,1));
                nU = length(units);
                assert(length(localSorting) == nU, 'must be identical');
                assert(length(localSortingID) == nU, 'must be identical');
                assert(size(T_merged,3) == nU, 'must be identical');
                
                disp('Saving postprocessing results...')
                save(self.destinationlocation.files.results, 'gdf_merged', 'gdf_discarded', 'T_merged', 'localSorting', 'localSortingID', 'G', '-v7.3');
            end
            
        end
        
        
        %         % -----------------------------------------------------------------
        %         function createSummaryFile(self)
        %
        %             self.destinationlocation.files.summary = fullfile(self.folders.main, 'summary.mat');
        %             try
        %                 load(self.destinationlocation.files.summary)
        %             catch
        %                 try
        %                     res = load(self.destinationlocation.files.results)
        %
        %                     units = unique( res.gdf_merged(:,1) );
        %                     nUnits = length(units);
        %                     %[dir_exists,mess,messid] = mkdir(self.folders.main, 'qchdsort.plot.');
        %                     %self.folders.qchdsort.plot. = fullfile( self.folders.main, 'qchdsort.plot.');
        %
        %                     %DSFull = hdsort.filewrapper.CMOSMEA(self.files.data);%, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);
        %                     %MES = DSFull.MultiElectrode.toStruct();
        %                     save(self.destinationlocation.files.summary, 'units', 'nUnits');
        %
        %                 catch
        %                     disp('Creation of summary file failed!');
        %                     return;
        %                 end
        %             end
        %         end
        
        
        % -----------------------------------------------------------------
        function runQC(self)
            try
                res = load(self.destinationlocation.files.results)
                
                [dir_exists,mess,messid] = mkdir(self.folders.main, 'qchdsort.plot.');
                self.folders.qchdsort.plot. = fullfile( self.folders.main, 'qchdsort.plot.');
                
                DSFull = hdsort.filewrapper.CMOSMEA(self.files.data);%, 'useFilter', self.sortJobP.useFilter, 'name', self.jobName);
                MES = DSFull.MultiElectrode.toStruct();
                
            catch
                disp('Quality control failed!');
                return;
            end
            
            figures = struct;
            
            %% ISIH:
            self.destinationlocation.files.isih = fullfile( self.folders.qchdsort.plot., 'isih');
            if ~exist(self.destinationlocation.files.isih, 'file')
                F = mysort.hdsort.plot.isi( res.gdf_merged )
                mysort.hdsort.plot.savefig(F.figureHandle, self.destinationlocation.files.isih)
            end
            
            %% Footprints whole
            self.destinationlocation.files.footprints = fullfile( self.folders.qchdsort.plot., 'footprints_whole');
            if ~exist(self.destinationlocation.files.footprints, 'file')
                F.figureHandle = figure();
                mysort.hdsort.plot.hdsort.waveforms.D(res.T_merged, MES.electrodePositions, 'IDs', res.localSortingID);
                mysort.hdsort.plot.savefig(F.figureHandle, self.destinationlocation.files.footprints)
            end
            
            %% Footprints localized
            self.destinationlocation.files.footprints2 = fullfile( self.folders.qchdsort.plot., 'footprints_localized');
            if ~exist(self.destinationlocation.files.footprints2, 'file')
                F.figureHandle = figure(); P.AxesHandle = []
                for i = 1:length(res.localSorting)
                    id = res.localSortingID(i);
                    lu = res.localSorting(i);
                    P = mysort.hdsort.plot.hdsort.waveforms.D(0.1*res.T_merged(:,:,i), MES.electrodePositions, 'IDs', (1000*id+lu), 'maxNumberOfChannels', 10, 'AxesHandle', P.AxesHandle, 'hdsort.plot.rgs', {'color', hdsort.plot.PlotInterface.vectorColor(i)});
                end
                mysort.hdsort.plot.savefig(F.figureHandle, self.destinationlocation.files.footprints2)
            end
            
        end
        
        % -----------------------------------------------------------------
        function copyBackResults(self, destinationlocation, copyEverything)
            % Copy results from the grid folder back to the destination folder
            
            if nargin < 3 copyEverything = false; end
            
            disp('Copy back result file(s)...');
            
            [dir_exists,mess,messid] = mkdir(destinationlocation);
            assert(dir_exists, 'Destination directory could not be created!');
            
            self.files.groupFile = self.copyFile(self.files.groupFile, destinationlocation);
            self.files.summary = self.copyFile(self.files.summary, destinationlocation);
            
            newGroupsFolder = fullfile(destinationlocation, 'groups');
            [dir_exists,mess,messid] = mkdir(newGroupsFolder); %self.folders.destination);
            assert(dir_exists, 'Destination groups directory could not be created!');
            
            if copyEverything
                disp('Copy back entire sorting output...')
                self.folders.groups = self.copyFolderContent(self.folders.groups, newGroupsFolder);
                
                disp('Copy back data files...')
                self.folders.data = self.copyFolder(self.folders.data, destinationlocation);
                disp('Done copying.')
            else
                %% Only copy back necessary (and small) output files:
                filesToCopy = {'_templates.mat', ...
                    '.P.mat',...
                    '.030spikes_det_merged.mat', ...
                    '.060cov.mat', ...
                    '.090clusters_meanshift.mat',...
                    '.100botm_matching.mat', ...
                    '.110clusters_meanshift_merged.mat'};
                
                groupFolders = self.findSubFolders(self.folders.groups)
                
                for i = 1:length(groupFolders)
                    groupFolder = groupFolders{i};
                    [pathstr,name,ext] = fileparts(groupFolder);
                    newGroupFolder = fullfile(newGroupsFolder, name);
                    
                    [dir_exists,mess,messid] = mkdir(newGroupFolder);
                    assert(dir_exists, 'Destination group folder could not be created!');
                    
                    for j = 1:length(filesToCopy)
                        fullFile = fullfile(self.folders.groups, groupFolder, [self.jobName filesToCopy{j}])
                        self.copyFile(fullFile, newGroupFolder);
                    end
                    
                    self.files.groupPaths{i} = groupFolder; %newGroupFolder;
                end
            end
            
            disp('Sorting results have been copied to:')
            self.folders.main = destinationlocation;
            disp(self.folders.main)
            
        end
        
        
        % -----------------------------------------------------------------
        function copyDataFilesToBinary(self)
            error('Not implemented properly!')
            
            % This function copies the original h5 file to a binary file and a
            % h5 metadata file in the job folder. From now on, meta-file will
            % be treated as the original data file.
            % The binary files are created directly on the scratch.
            newDataFileLocations = {};
            
            disp('Copy data file(s) to binary...');
            for i = 1:length(self.files.data)
                
                [pathstr,name,ext] = fileparts(self.files.data{i});
                
                fileNameH5 = fullfile(pathstr, [name '_meta' ext]);
                fileNameBin = fullfile(pathstr, [name '.dat']);
                
                if exist(fileNameBin, 'file') ~= 2
                    mysort.mea.copyH5toBinary(self.files.data{i}, fileNameH5, fileNameBin);
                end
                
                %% Remember the location of the datafiles on the scratch now:
                newDataFileLocations{i} = fileNameH5;
                
                disp([ num2str(i) ' of ' num2str(length(self.files.data)) ' copied.']);
            end
            
            self.files.data = newDataFileLocations;
        end
        
        % -----------------------------------------------------------------
        function copyDataFiles(self, destinationLocation)
            
            newDataFileLocations = {};
            
            disp('Copy data file(s)...');
            for i = 1:length(self.files.data)
                
                %% Check if the data (the first file at least) is binary:
                m_test = hdsort.filewrapper.CMOSMEA(self.files.data{i});
                
                if ~m_test.isBinaryFile()
                    [pathstr,name,ext] = fileparts(self.files.data{i});
                    fileNameH5 = fullfile(destinationLocation, [name '.h5']);
                    fileNameBin= fullfile(destinationLocation, [name '.dat']);
                    
                    %% Create a binary file on the scratch:
                    if exist(fileNameBin, 'file') ~= 2
                        mysort.mea.copyH5toBinary(self.files.data{i}, fileNameH5, fileNameBin);
                    end
                    newDataFileLocations{i} = fileNameH5;
                else
                    [pathstr,name,ext] = fileparts(self.files.data{i});
                    fileNameH5 = self.files.data{i}; %fullfile(self.startlocation.folders.data, [name '_meta' ext]);
                    fileNameBin = fullfile(pathstr, [name '.dat']);
                    
                    %% Copy meta and binary file:
                    newDataFileLocations{i} = self.copyFile(fileNameH5, destinationLocation) ;
                    self.copyFile(fileNameBin,destinationLocation);
                end
                
                disp([ num2str(i) ' of ' num2str(length(self.files.data)) ' copied.']);
            end
            self.files.data = newDataFileLocations;
            
        end
        
        %% Auxilary function for finding non-hidden subfolders
        function list = findSubFolders(self, folder)
            d = dir(folder);
            isub = [d(:).isdir]; %# returns logical vector
            list = {d(isub).name}';
            list(ismember(list,{'.','..'})) = [];
        end
    end
    
    methods(Static)
        %------------------------------------------------------------------
        %------------------------- RUN  FUNCTION ---------------------------
        %------------------------------------------------------------------
        function run(taskFile, debugFlag)
            if nargin < 2
                debugFlag = false;
            end
            
            %% List of default parameters:
            sortP = struct;
            sortP.botm.Tf = 75;
            sortP.botm.cutLeft = 20;
            sortP.botm.run = 0;
            sortP.spikeCutting.blockwise = false;
            
            taskP = struct;
            taskP.runName = 'BOTMdefault';
            
            %% Load taskFile:
            T = load(taskFile);
            
            %% Write "taskParameters" to struct "taskP" and "sortP":
            sortP = hdsort.util.mergeStructs(sortP, T.taskParameters.sortingParameters);
            taskParameters = rmfield(T.taskParameters,'sortingParameters');
            
            taskP = hdsort.util.mergeStructs(taskP, T.taskParameters);
            clear T;
            
            %% Check necessary parameters:
            assert( isfield(taskP, 'dataFiles'), 'Task aborted: field taskParameters.dataFiles not specified!');
            for f = 1:length(taskP.dataFiles)
                % Note: exist(...) == 2: name is a full pathname to a file
                assert( exist(taskP.dataFiles{f}, 'file') == 2, ['Task aborted: dataFile ' taskP.dataFiles{f} ' not found!']);
            end
            assert( exist(taskP.outputPath, 'dir') == 7, ['Task aborted: no valid outputPath specified! Path given:' taskP.outputPath]);
            assert( isfield(taskP, 'groupidx'), 'Task aborted: field taskParameters.groupidx not specified!');
            assert( isfield(taskP, 'taskID'), 'Task aborted: field taskParameters.taskID not specified!');
            
            %% (Re-)Set reporting file:
            rep = hdsort.filewrapper.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
            rep(:,:) = [0 0];
            
            if ~debugFlag
                try
                    mainBlock();
                catch ME
                    errorHandling(ME);
                end
            else
                mainBlock();
            end
            
            function mainBlock()
                disp('Creating DS...')
                %% Sort:
                DS = hdsort.filewrapper.CMOSMEA(taskP.dataFiles, 'useFilter', 0, 'name', 'PREFILT');
                MES = DS.MultiElectrode.toStruct();
                
                DS.restrictToChannels(taskP.groupidx);
                disp('Start sorting...')
                [S P_] = mysort.sorters.sort(DS, fullfile(taskP.outputPath, taskP.groupPath ), taskP.runName, sortP);
                
                if S.STOP_ME_BECAUSE_I_AM_SLOW return; end
                
                %% RELEASE CHANNEL RESTRICTIONS FOR TEMPLATE ESTIMATION
                DS.restrictToChannels();
                disp('Starting template estimation...')
                hdsort.leg.startHDSortingTemplateEstimation(taskP.outputPath, taskP.groupPath, taskP.runName, sortP.botm.Tf, sortP.botm.cutLeft, DS, taskP.groupidx, MES);
                
                %% Write to reporter file:
                disp('Writing results...')
                rep = hdsort.filewrapper.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [1 0];
            end
            
            function errorHandling(ME)
                
                disp('Catch error...')
                errStr = hdsort.util.buildLastErrString(ME);
                disp(errStr)
                
                rep = hdsort.filewrapper.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [0 1];
                rethrow(ME)
            end
        end
        
    end
end
