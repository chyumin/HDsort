%%
pd = pdefs();

rootFolder = pd.eulerRoot

dj = grid.DummyJob('dummyjob04', rootFolder, 'onEuler', true)
%%
dj.setTaskParameters();
dj.prepareTasks();
%%
dj.createAutoSubmitToken();
dj.waitForTasksToFinish(5);



%
% #!/bin/bash
% #
% #BSUB -P test_project         # project code
% #BSUB -J dummy_job_name       # job name
% #BSUB -W 00:10                # wall-clock time (hrs:mins)
% #BSUB -n 3                    # number of tasks in job
% #BSUB -q normal               # queue
% #BSUB -e errors.%J.hybrid     # error file name in which %J is replaced by the job ID
% #BSUB -o output.%J.hybrid     # output file name in which %J is replaced by the job ID
% #BSUB -R "rusage[mem=MMMM]"   # memort in MB MMMM
%
% matlab -nojvm -singleCompThread -nodisplay -r "grid.GridJob.runBSubTask('~/data/dummyjob01/taskFile/taskFile'); exit();"
% 
% #!/bin/bash
% #
% #BSUB -W 00:02                # wall-clock time (hrs:mins)
% #BSUB -n 1                    # number of tasks in job
% #BSUB -J dummyjob01[1-3]      # name and task id range in brackets
% #BSUB -P dummytestproject     # project code
% #BSUB -e /cluster/home/rolandd/data/dummyjob01/log/e.%J.%I     # error file name in which %J is replaced by the job ID 
% #BSUB -o /cluster/home/rolandd/data/dummyjob01/log/o.%J.%I     # output file name in which %J is replaced by the job ID
% # #BSUB -R "rusage[mem=MMMM]"   # memort in MB MMMM
% 
% 
% matlab -nojvm -singleCompThread -nodisplay -r "grid.GridJob.runBSubTask('~/data/dummyjob01/taskFiles/taskFile'); exit();"


