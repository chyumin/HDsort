% Things necessary before you can run this script:
% TO DO ONCE:
% A) Become a Grid User -> talk to IT
% B) Make a checkout of the matlab SVN folder on your new Grid Home
% B-1) For this you will need to set up your SSH key on the submit host to
%      be able to access the SVN
% C) copy the startup.m from the matlab folder to your Grid Home and set
%    the correct paths
% D) Test this by logging into the submit system, loading the matlab
%    module, starting matlab and see if the hdsort-package is loaded
% E) ...
%
%
% TO DO EACH TIME YOU RUN A NEW JOB
% A) ...
classdef GridJob < handle
    properties (SetAccess=private)
    end
    
    properties
        gridConfig
        jobName
        allTaskParameters
        
        logFolder
        
        files
        folders
        
        m_nTasksCompletedLastTime
        m_tTimeLastChange
        m_tTimeStart
        taskIDs
        taskType
        
        memoryusage
        runtime
        startIndex
        endIndex
        gridType
        nWorkersPerNode
        
        P
        P_untreated
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = GridJob(job_name, rootFolder, varargin)
            % Starts a grid job by running a number of parallel matlab
            % calls on a series of config files. Each individual matlab task
            % will load a config file which is unique for this task and run
            % the same script as all other tasks as well.
            self.gridConfig = hdsort.grid.config();
            
            is_letter = isstrprop(job_name, 'alpha');
            assert(is_letter(1), 'Job name must begin with a letter!')
            
            p.runLocallyOn = [];
            p.queue = 'regular';
            p.resume = false;
            p.runtime_hours = 0;
            p.runtime_minutes = 10;
            p.memoryusage = [];
            p.gridType = 'QSUB';
            
            [p p2] = hdsort.util.parseInputs(p, varargin, 'split');
            self.P_untreated = p2; self.P = p;
            
            % For debugging on a local machine:
            if ~isempty(self.P.runLocallyOn) self.gridConfig.home = self.P.runLocallyOn; end
            
            self.jobName = job_name;
            self.folders.root = rootFolder;
            self.runtime.hours = self.P.runtime_hours;
            self.runtime.minutes = self.P.runtime_minutes;
            self.memoryusage = self.P.memoryusage;
            self.gridType = self.P.gridType;
            self.nWorkersPerNode = 1; % default that can be overwritten by derivative classes
            
            self.createFolders();
            
        end
        
        % -----------------------------------------------------------------
        function createFolders(self)
            [dir_exists,mess,messid] = mkdir(self.folders.root, self.jobName);
            assert(dir_exists, 'Job directory could not be created!');
            self.folders.main = fullfile(self.folders.root, self.jobName);
            
            % LOG folder:
            self.folders.log = fullfile(self.folders.main, 'log');
            [dir_exists,mess,messid] = mkdir(self.folders.log);
            assert(dir_exists, 'Job directory could not be created!');
            
            % Task files folder:
            self.folders.tasks = fullfile(self.folders.main, 'taskFiles');
            [dir_exists,mess,messid] = mkdir(self.folders.tasks);
            assert(dir_exists, 'Task file directory could not be created!');
            
            % Report files folder:
            self.folders.report = fullfile( self.folders.main, 'reportFiles');
            [dir_exists,mess,messid] = mkdir(self.folders.report);
            assert(dir_exists, 'Report file directory could not be created!');
        end
        
        % -----------------------------------------------------------------
        function prepareTasks(self)
            self.createReportFiles();
            self.createTaskFiles();
            self.constructJobSH();
            
            % Set resume value to true
            self.setResume(true);
        end
        
        % -----------------------------------------------------------------
        function setResume(self, r)
            self.P.resume = r;
        end
        
        % -----------------------------------------------------------------
        function fileName = saveToFile(self)
            fileName = fullfile(self.folders.main,'GridJob.mat');
            save(fileName,'self');
        end
        
        % -----------------------------------------------------------------
        % This function creates a file for each task, containing only one
        % structure named taskParameters. The class element
        % allTaskParameters{:} has to be specified before.
        function createTaskFiles(self)
            
            self.files.tasks = fullfile( self.folders.tasks, '/taskFile');
            
            assert( self.nTasks == length(self.allTaskParameters), 'Error: number of tasks does not correspond to size of allTaskParameters!');
            for ii = 1:self.nTasks
                taskParameters = self.allTaskParameters{ii};
                
                taskParameters.reportFile = hdsort.grid.GridJob.convertToLinux(self.files.report{ii});
                taskParameters.reportFolder = hdsort.grid.GridJob.convertToLinux(self.folders.report);
                
                taskType = self.taskType;
                save( [self.files.tasks num2str(taskParameters.taskID)], 'taskParameters', 'taskType');
            end
        end
        
        % -----------------------------------------------------------------
        function createReportFiles(self)
            for ii = 1:self.nTasks
                self.files.report{ii} = fullfile( self.folders.report, ['task' num2str(self.taskIDs(ii)) '.dat']);
                
                % Only create new report files when it doesn't exist yet:
                if exist(self.files.report{ii}, 'file') ~= 2
                    rep = hdsort.filewrapper.util.BinaryFileMatrix(self.files.report{ii}, [1 2], 'writable', true);
                    rep(:,:) = [0 0];
                end
            end
        end
        
        % -----------------------------------------------------------------
        function out = waitForTasksToFinish(self, pause_duration)
            %Wait for all tasks to finish
            % OUT = waitForTasksToFinish(self, PAUSE_DURATION=10)
            % This function regularly checks if the tasks have finished and
            % displays the number of successfully finished tasks and the
            % number of errors.
            %  - returns 'true' when all tasks have completed.
            %  - returns 'false' when there has been at least one error.
            %
            %  PAUSE_DURATION: waiting time between loops in seconds.
            
            if nargin < 2 pause_duration = 10; end
            disp([self.jobName ': Waiting for tasks to finish (wait for ' num2str(pause_duration) ' seconds)...']);
            
            while true
                [nCompleted, tasksNotCompleted, nErrors, tasksWithErrors] = self.summarySnapshot(true);
                diff_completed = nCompleted - self.m_nTasksCompletedLastTime;
                self.m_nTasksCompletedLastTime = nCompleted;
                try
                    dt = toc(self.m_tTimeLastChange);
                catch
                    dt = 0;
                end
                minutes = floor(dt/60);
                sec = round(dt-minutes*60);
                dtstr = sprintf('%d min %d sec', minutes, sec);
                if diff_completed>0
                    self.m_tTimeLastChange = tic;
                end
                fprintf('%s: %d (+ %d) out of %d tasks completed. (Time to last Change: %s)\n', self.jobName, nCompleted, diff_completed, self.nTasks, dtstr);
                
                if self.nTasks - nCompleted < 10
                    disp(['Waiting for tasks: ' num2str(tasksNotCompleted)]);
                end
                
                if nErrors > 0
                    disp(['Errors in tasks: ' num2str(tasksWithErrors)]);
                end
                
                if nCompleted + nErrors == self.nTasks break; end
                pause(pause_duration);
            end
            
            out = true;
            if ~nErrors
                disp([self.jobName ': All tasks completed successfully...']);
            else
                disp([self.jobName ': ATTENTION: Waiting ended without completing all tasks!']);
                disp([num2str(nErrors) ' errors detected!']);
                out = false;
            end
            self.summarizeReports();
        end
        
        % -----------------------------------------------------------------
        function [nCompleted, tasksNotCompleted, nErrors, tasksWithErrors] = summarySnapshot(self, bNoCmdLineOutput)
            nCompleted = 0;
            nErrors = 0;
            tasksNotCompleted  = [];
            tasksWithErrors = [];
            if nargin == 1
                myDisp = @(x) disp(x);
            elseif bNoCmdLineOutput
                myDisp = @(x) 0;
            end
            for ii = 1:self.nTasks
                try
                    rep = hdsort.filewrapper.util.BinaryFileMatrix(self.files.report{ii}, [1 2], 'writable', false);
                catch
                    myDisp([self.jobName ': Binary file for task ' num2str( self.taskIDs(ii) ) ' threw exception!']);
                    rep = [0 1];
                end
                
                if rep(1,2) > 0
                    myDisp([self.jobName ': Error in task ' num2str( self.taskIDs(ii) ) ]);
                    nErrors = nErrors + 1;
                    tasksWithErrors = [tasksWithErrors self.taskIDs(ii)];
                end
                if rep(1,1) > 0
                    nCompleted = nCompleted + 1;
                else
                    tasksNotCompleted = [tasksNotCompleted, self.taskIDs(ii)];
                end
            end
            str = sprintf('%s: %d out of %d tasks completed. \n', self.jobName, nCompleted, self.nTasks);
            myDisp(str);
        end
        
        % -----------------------------------------------------------------
        function [bFinished, nErrors, nUnfinished, nCompleted, nTotal] = isTaskFinished(self, bNoCmdLineOutput)
            % Checks the status of all tasks
            nTotal = self.nTasks;
            bFinished = true;
            if nargin == 1
                bNoCmdLineOutput = false;
            end
            
            [nCompleted, tasksNotCompleted, nErrors, tasksWithErrors] = self.summarySnapshot(bNoCmdLineOutput);
            
            nUnfinished = nTotal - nErrors - nCompleted;
            
            if nUnfinished > 0
                bFinished = false;
            end
        end
        
        % -----------------------------------------------------------------
        function summarizeReports(self, cmdLineOutput)
            % This function reads all reports files and creates a summary of
            % them.
            if nargin == 1
                cmdLineOutput = false;
            end
            
            %if  ~cmdLineOutput
            %    return
            %end
            
            %% Check for most recent job id:
            try
                currentJobId = fullfile(self.folders.report, 'currentJobID.txt')
                fid = fopen(currentJobId);
                job_id_str = fread(fid, '*char')'; fclose(fid);
            catch
                job_id_str = 'undetermined_job';
            end
            
            %% Create a summary file for this job:
            self.files.summary = fullfile(self.folders.main, [job_id_str '_summary.txt'] );
            self.files.summaryMAT = fullfile(self.folders.main, [job_id_str '_summary.mat'] );
            summary = struct;
            
            if cmdLineOutput
                myDisp = @(varargin) fprintf(varargin{:});
            else
                fid = fopen(self.files.summary, 'a+');
                myDisp = @(fid, varargin) fprintf(fid, varargin{:});
            end
            
            myDisp(['Job: ' job_id_str '\n'])
            
            summary.tasksError = [];
            summary.nCompleted = 0;
            summary.completedTasks = zeros(1, self.nTasks);
            summary.times = {};
            for ii = 1:self.nTasks
                rep = hdsort.filewrapper.util.BinaryFileMatrix(self.files.report{ii}, [1 2], 'writable', false);
                if rep(1,2) > 0
                    summary.tasksError = [summary.tasksError self.taskIDs(ii)];
                end
                if rep(1,1) > 0
                    summary.nCompleted = summary.nCompleted + 1;
                    summary.completedTasks(ii) = 1;
                end
                
                %% Extract the time of this task:
                try
                    startTimeFile = fullfile(self.folders.report, ['startTime'  job_id_str '_'  num2str(self.taskIDs(ii)) '.mat']);
                    load(startTimeFile);
                    
                    endTimeFile = fullfile(self.folders.report, ['endTime'  job_id_str '_' num2str(self.taskIDs(ii)) '.mat']);
                    load(endTimeFile);
                    
                    summary.times{ii, 1} = startTime;
                    summary.times{ii, 2} = endTime;
                    summary.times{ii, 3} = etime(endTime, startTime);
                catch
                    summary.times{ii, 1} = clock;
                    summary.times{ii, 2} = clock;
                    summary.times{ii, 3} = 0;
                end
            end
            summary.totalDuration = self.findTotalDuration( summary.times );
            myDisp( ['Total Duration (in seconds): ' num2str(summary.totalDuration) '\n']);
            
            myDisp( ['Number of tasks completed: ' num2str(summary.nCompleted) '\n']);
            myDisp( ['Number of errors:          ' num2str( length(summary.tasksError) ) '\n']);
            
            if ~isempty(summary.tasksError)
                myDisp( ['Errors in tasks: ' num2str( summary.tasksError ) '\n']);
            end
            if any(summary.completedTasks==0)
                myDisp( ['Not completed tasks: ' num2str( find(summary.completedTasks==0) ) '\n']);
            end
            if ~cmdLineOutput
                fclose(fid);
                save(self.files.summaryMAT, 'summary');
            end
            
        end
        
        % -----------------------------------------------------------------
        function totalDuration = findTotalDuration(self, times )
            minTime = times{1,1};
            maxTime = times{1,2};
            for ii = 1:size(times, 1)
                startT = times{ii,1};
                stopT = times{ii,2};
                
                if etime(startT, minTime) < 0
                    minTime = startT;
                end
                if etime(stopT, maxTime) > 0
                    maxTime = stopT;
                end
            end
            
            totalDuration = etime(maxTime, minTime);
        end
        
        % -----------------------------------------------------------------
        function str = constructConsoleJobStartCommand(self)
            str = sprintf('qsub -q %s.q %s_job.sh', self.P.queue, self.jobName);
        end
        
        % -----------------------------------------------------------------
        function filename = constructJobSH(self)
            
            str = '#!/bin/bash\n';
            
            if strcmp(self.gridType, 'QSUB')
                str = [str '#$ -V\n'];
                str = [str '#$ -cwd\n'];
                
                %if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                    %warning('Sorting not started from a linux machine might cause problems!')
                    logFolder = hdsort.grid.GridJob.convertToLinux(self.folders.log);
                    str = [str '#$ -o ' logFolder '\n'];
                    str = [str '#$ -e ' logFolder '\n'];
                %else
                %    str = [str '#$ -o ' self.folders.log '\n'];
                %    str = [str '#$ -e ' self.folders.log '\n'];
                %end
                
                str = [str '#$ -q ' self.P.queue '.q\n'];
                if strcmp(self.P.queue, 'regular')
                    % do nothing
                elseif strcmp(self.P.queue, 'immediate')
                    str = [str '#$ -l immediate\n'];
                elseif strcmp(self.P.queue, 'background')
                    str = [str '#$ -l background\n'];
                else
                    error('Unknown Queue!');
                end
                str = [str sprintf('#$ -t %d-%d\n', self.startIndex, self.endIndex)];
                %if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                    taskFile = hdsort.grid.GridJob.convertToLinux(self.files.tasks)
                    str = [str sprintf('matlab -nodisplay -r "hdsort.grid.GridJob.runQSubTask(''%s\''); exit();"', taskFile)];
                %else
                %    str = [str sprintf('matlab -nodisplay -r "hdsort.grid.GridJob.runQSubTask(''%s\''); exit();"', self.files.tasks)];
                %end
                
            elseif strcmp(self.gridType, 'BSUB')
                str = [str sprintf('#BSUB -W %d:%d\n', self.runtime.hours, self.runtime.minutes)]; %# wall-clock time (hrs:mins)
                str = [str '#BSUB -n ' num2str(self.nWorkersPerNode) '\n']; % # number of tasks in job
                str = [str '#BSUB -J ' self.jobName sprintf('[%d-%d]', self.startIndex, self.endIndex) '\n' ]; % # name and task id range in brackets
                str = [str '#BSUB -o ' self.folders.log '/o.%%J.%%I\n']; % # output file name in which %J is replaced by the job ID
                str = [str '#BSUB -e ' self.folders.log '/e.%%J.%%I\n']; % # error file name in which %J is replaced by the job ID
                
                if ~isempty(self.memoryusage)
                    str = [str sprintf('#BSUB -R "rusage[mem=MMMM]"', self.memoryusage)]; % # memort in MB MMMM
                end
                
                str = [str 'cd ~/trunk/matlab\n'];
                if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                    error('You cannot create a bsub on a windows or mac machine yet!')
                else
                    str = [str sprintf('matlab -nodisplay -r "hdsort.grid.GridJob.runBSubTask(''%s\''); exit();"', self.files.tasks) '\n'];
                end
            else
                error(['GridType ' self.gridType ' unknown!'])
            end
            
            filename = fullfile(self.folders.main, [self.jobName '_job.sh']);
            fh = fopen(filename, 'w+');
            fprintf(fh, str );
            fclose(fh);
        end
        
        
        
        % -----------------------------------------------------------------
        function copyDataFilesToScratch(self)
            error('Not supported anymore!')
            
            disp('Copy data file(s)...');
            for i = 1:length(self.startlocation.files.data)
                self.scratchlocation.files.data{i} = self.copyFile(self.startlocation.files.data{i}, self.scratchlocation.folders.data);
                disp([ num2str(i) ' of ' num2str(length(self.startlocation.files.data)) ' copied.']);
            end
            
        end
        
        
        % -----------------------------------------------------------------
        % Copy one file to destination folder and return new file name:
        function fileName = copyFile(self, fileName, destFolder, forceCopy)
            if nargin < 4 forceCopy = false; end
            flags = '-q';%  '-vP'
            
            if forceCopy
                system(['rsync ' flags ' ' fileName ' ' destFolder]);
            else
                system(['rsync ' flags ' --ignore-existing ' fileName ' ' destFolder]);
            end
            
            [pathstr,name,ext] = fileparts(fileName);
            fileName = fullfile(destFolder, [name ext]);
        end
        
        % -----------------------------------------------------------------
        function folderName = copyFolder(self, folderName, destFolder, forceCopy)
            if nargin < 4 forceCopy = false; end
            
            if forceCopy
                system(['rsync -vP -r' folderName ' ' destFolder]);
            else
                system(['rsync -vP -r --ignore-existing ' folderName ' ' destFolder]);
            end
            
            folderName = fullfile(destFolder, folderName);
        end
        
        % -----------------------------------------------------------------
        function folderName = copyFolderContent(self, folderName, destFolder, forceCopy)
            if nargin < 4 forceCopy = false; end
            
            if forceCopy
                system(['rsync -vP -r' folderName '/* ' destFolder]);
            else
                system(['rsync -vP -r --ignore-existing ' folderName '/* ' destFolder]);
            end
            
            folderName = fullfile(destFolder);
        end
        
        % -----------------------------------------------------------------
        function N = nTasks(self)
            N = self.endIndex - self.startIndex + 1;
        end
        
        
        % -----------------------------------------------------------------
        %              FOR DEBUGGING:
        % -----------------------------------------------------------------
        function runTaskLocally(self, task_id)
            if 0
                setenv('SGE_TASK_ID', num2str(task_id))
                
                setenv('JOB_ID', '0')
                jobFolder = fullfile(self.folders.report, '0')
                try
                    rmdir(jobFolder, 's');
                catch
                end
                mkdir(jobFolder);
                hdsort.grid.GridJob.runTask(self.files.tasks, false);
            else
                task_id_str = num2str(task_id);
                job_id_str = '0';
                jobFolder = fullfile(self.folders.report, job_id_str)
                debugFlag = false;
                try
                    rmdir(jobFolder, 's');
                catch
                end
                mkdir(jobFolder);
                hdsort.grid.GridJob.runTask(self.files.tasks, task_id_str, job_id_str, debugFlag);
            end
        end
        
        
        % -----------------------------------------------------------------
        function prepareTest(self)
            self.startIndex = 3;
            self.endIndex = 12;
            self.files.tasks = fullfile( self.gridConfig.home, self.jobName, '/taskFile');
            
            var1 = 'var1';
            var2 = 'var2';
            num1 = 12345;
            nTasks = self.numTasks();
            outputPath = fullfile( self.gridConfig.home, self.jobName, 'results');
            
            for task_id = self.startIndex:self.endIndex
                save( [self.files.tasks num2str(task_id)], 'outputPath', 'var1', 'var2', 'num1', 'task_id', 'nTasks');
            end
            
            test = 1;
            self.constructJobSH(test);
        end
        
        % -----------------------------------------------------------------
        function createAutoSubmitToken(self)
            token_file = fullfile(self.gridConfig.tokenFilesFolder, ['start_' self.jobName '.mat']);
            
            shFile = fullfile(self.folders.main, [self.jobName '_job.sh']);
            %if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                shFile = hdsort.grid.GridJob.convertToLinux(shFile);
            %end
            
            save(token_file, 'shFile');
            disp('Token file created.')
            while true
                pause(2);
                if ~exist(token_file, 'file')
                    disp('Token file used and moved.')
                    break;
                end
            end
        end
        
    end
    
    methods (Abstract, Static)
        run(taskFile)
    end
    
    methods (Static)
        
        % -----------------------------------------------------------------
        function linuxPath = convertToLinux(path)
            % This function transforms absolute file names from a Windows
            % or Mac system to a linux system. The necessary settings are
            % must be set in +hdsort/+grid/config.m.
            % The path is reduced to a string called
            % 'gridConfig.commonFolderName' and then rebuilt based on
            % 'gridConfig.linuxSortingPath'.
            if iscell(path)
                linuxPath = {};
                for ii = 1:numel(path)
                    linuxPath{ii} = hdsort.grid.GridJob.convertToLinux(path{ii});
                end
                return
            end
            linuxPath = path;
            
            if ~isempty(strfind(computer, 'WIN')) || ~isempty(strfind(computer, 'MACI64'))
                gridConfig = hdsort.grid.config();
                
                p = regexp(path, filesep, 'split');
                while ~strcmp(p{2}, gridConfig.commonFolderName)
                    p = {p{2:end}};
                end
                
                linuxPath = gridConfig.linuxSortingPath;
                while ~isempty(p)
                    linuxPath =  [linuxPath '/' p{1}];
                    p = {p{2:end}};
                end
            end
        end
        
        % -----------------------------------------------------------------
        function localPath = convertToLocal(path)
            if iscell(path)
                localPath = {};
                for ii = 1:numel(path)
                    localPath{ii} = hdsort.grid.GridJob.convertToLocal(path{ii});
                end
                return
            end
            
            gridConfig = hdsort.grid.config();
            
            p = regexp(path, filesep, 'split');
            while ~strcmp(p{2}, 'Mea1k')
                p = {p{2:end}};
            end
            
            localPath = gridConfig.localSortingPath;
            while ~isempty(p)
                localPath =  [localPath '/' p{1}];
                p = {p{2:end}};
            end
        end
        
        % -----------------------------------------------------------------
        function out = runQSubTask(taskName, debugFlag)
            if nargin < 2 debugFlag = false; end
            task_id_str = getenv('SGE_TASK_ID')
            job_id_str = getenv('JOB_ID')
            taskName
            out = hdsort.grid.GridJob.runTask(taskName, task_id_str, job_id_str, debugFlag)
        end
        
        % -----------------------------------------------------------------
        function out = runBSubTask(taskName, debugFlag)
            if nargin < 2 debugFlag = false; end
            task_id_str = getenv('LSB_JOBINDEX');
            job_id_str = getenv('LSB_JOBID');
            out = hdsort.grid.GridJob.runTask(taskName, task_id_str, job_id_str, debugFlag);
        end
        
        % -----------------------------------------------------------------
        function out = runTask(taskName, task_id_str, job_id_str, debugFlag)
            %Static function of GridJob called by each node on the grid.
            % OUT = hdsort.grid.GridJob.runTask(TASKNAME, DEBUGFLAG)
            %
            % OUT is boolean indicating success or failure of task.
            disp(['Starting hdsort.grid.GridJob.runTask(' taskName ', ' task_id_str ', ' job_id_str ') ...'])
            
            if nargin < 4
                debugFlag = false;
            end
            
            %% Find task-ID and load its taskFile:
            assert( ~isempty(task_id_str), 'No task_id_str found!');
            taskFile = fullfile([taskName task_id_str '.mat']);
            f = load(taskFile);
            
            %% Check if the same task in the same job already exists:
            processFolder = fullfile( f.taskParameters.reportFolder, job_id_str)
            [dir_exists,mess,messid] = mkdir(processFolder);
            uniqueProcessFile = fullfile( processFolder, ['task' task_id_str 'exists.dat'])
            
            if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                uniqueProcessFile = hdsort.grid.GridJob.convertToLocal(uniqueProcessFile);
                f.taskParameters.reportFolder = hdsort.grid.GridJob.convertToLocal(f.taskParameters.reportFolder);
            end
            
            if exist(uniqueProcessFile, 'file') == 2
                error(['Task ' task_id_str ' in job ' job_id_str ' already exists!'])
            end
            fileID = fopen(uniqueProcessFile,'w+'); fwrite(fileID, zeros(1), 'int8'); fclose(fileID);
            disp('uniqueProcessFile created.')
            
            %% Create a file that stores the current job_id:
            fileID = fopen( fullfile(f.taskParameters.reportFolder, 'currentJobID.txt'),'w+'); fwrite(fileID, job_id_str); fclose(fileID);
            clear fileID
            disp('currentJobID.txt file created.')
            
            %% Create a file that stores the start time vector:
            startTime = clock;
            save(  fullfile(f.taskParameters.reportFolder, ['startTime' job_id_str '_' task_id_str '.mat']), 'startTime' );
            disp(['startTime' task_id_str '.mat file created.'])
            
            hdsort.grid.(f.taskType).run(taskFile, debugFlag)
            out = true;
            
            %% Create a file that stores the current time vector:
            endTime = clock;
            save(  fullfile(f.taskParameters.reportFolder, ['endTime'  job_id_str '_'  task_id_str '.mat']), 'endTime' );
            disp(['endTime' task_id_str '.mat file created.'])
            
            disp('hdsort.grid.GridJob.runTask() finished.')
        end
        
        
        function submit_daemon(gridType)
            if nargin < 1
                gridType = 'QSUB';
            end
            
            gridConfig = hdsort.grid.config();
            
            cd(gridConfig.tokenFilesFolder)
            while 1
                fnames = dir(fullfile(gridConfig.tokenFilesFolder, 'start_*.mat'));
                if ~isempty(fnames)
                    disp('Found token, processing...')
                    for i=1:length(fnames)
                        token_file = fullfile(gridConfig.tokenFilesFolder, fnames(i).name);
                        t = load(token_file);
                        submit_token_file = fullfile(gridConfig.tokenFilesFolder, 'submitted', fnames(i).name);
                        disp(token_file)
                        
                        disp(['Sorting shFile: ' t.shFile])
                        
                        if strcmp(gridType, 'BSUB')
                            submit_str = sprintf('bsub < %s', t.shFile);
                        elseif strcmp(gridType, 'QSUB')
                            submit_str = sprintf('qsub %s', t.shFile);
                        end
                        move_str   = sprintf('mv %s %s', token_file, submit_token_file);
                        
                        cd('~')
                        hdsort.util.logToFile(gridConfig.log_file, submit_str)
                        [status, result] = system(submit_str);
                        hdsort.util.logToFile(gridConfig.log_file, result)
                        cd(gridConfig.tokenFilesFolder)
                        
                        if status == 0
                            disp('Submit successful')
                            disp(result)
                        else
                            disp('Submit failed')
                        end
                        
                        hdsort.util.logToFile(gridConfig.log_file, move_str)
                        [status, result] = system(move_str);
                        hdsort.util.logToFile(gridConfig.log_file, result)
                        pause(2.5)
                        disp('Done processing. Waiting...')
                    end
                end
                pause(10)
            end
            disp('Done.')
        end
    end
end
