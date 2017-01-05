function def = pathDefinitions()
    bIsOffline = true;
    bIsOffline = false; 
%     this_path = pwd;
%     offlineFile = fullfile(this_path, 'offline_status_indicator.txt');
%     if exist(offlineFile, 'file')
%         bIsOffline = true;
%     end
    def.hdsort.pathDefinitionsPath = pwd;
    def.callerPath = [];
    
    user_config = hdsort.global_user_config();
    user_name = user_config.user_name;
    
%     if strcmp(user_name, 'frankef');
%         disp('Hey, sunshine!');
%     elseif strcmp(user_name, 'rolandd')
%         disp('M&Ms !!!!!!!!!');
%     end
    
    if ~isempty(strfind(computer, 'WIN'))
        def.sortingOutPath = 'C:\localData\sortingOut\';
        def.networkTempShare = 'S:\group\hierlemann\Temp\FelixFranke\';
        def.localData = 'C:\localData\';
        def.localDataOnNetwork = fullfile(def.networkTempShare, 'LocalData');
        def.serverData = 'S:\group\hierlemann\recordings\collaborations\';
        def.serverDataRoska = 'S:\group\hierlemann\recordings\HiDens\Roska';  
        def.svnDataRoska = '';
        def.ravaDropBox = 'C:\Users\frankef\Dropbox\Paper NoiseCorrelations\';
        def.ravaSimulations = fullfile(def.localDataOnNetwork, 'RavaSimulations');
        def.serverSortingOut = 'S:\group\hierlemann\recordings\SpikeSortingOut\';
        def.micheleDSPaperDataPath = 'S:\group\hierlemann\recordings\HiDens\Roska\Fiscella_2014\backup_results';
        def.analysedRoland = 'S:\group\hierlemann\AnalyzedData\Mea1k\rolandd';
        
        def.mea1kData = 'S:\group\hierlemann\recordings\Mea1k';
        def.mea1kIntermediate = 'S:\group\hierlemann_Mea1k\intermediate_data\Mea1k\';
        def.tokenFiles = [def.mea1kIntermediate user_name '\tokens\'];
        
        if bIsOffline
            def.ravaSimulations = 'C:\LocalData\RavaSimulations\';
            def.localDataOnNetwork = def.localData;
            def.micheleDSPaperDataPath = 'C:\LocalData\Michele\backup_results';
        end
        
    elseif ~isempty(strfind(computer, 'MACI64'))
        def.hima01 = fullfile('/Volumes', 'hierlemann-1'); %% Normal storage
        def.hima02 = fullfile('/Volumes', 'hierlemann'); %% Hima large storage (not backuped!)
        
        def.mea1kRoot = fullfile(def.hima02, 'recordings', 'Mea1k');
        def.mea1kIntermediate = fullfile(def.hima02, 'intermediate_data', 'Mea1k');
        def.mea1kAnalyzed = fullfile(def.hima01, 'AnalyzedData', 'Mea1k');
        def.Mea1kShared = fullfile(def.hima02, 'recordings', 'Mea1k', 'shared');
        
        def.networkTempShare = fullfile('/Volumes', 'hierlemann-1', 'Temp', 'FelixFranke');
        
        def.localData = '/Users/rolandd/tmp/';
        def.serverData = fullfile(def.hima01, 'recordings','collaborations');        
        def.localDataOnNetwork = fullfile(def.networkTempShare, 'LocalData');
        
        def.tokenFiles = fullfile(def.mea1kIntermediate, user_name, 'tokens');
         
    else
        def.hima01 = fullfile('/net', 'bs-filesvr01', 'export', 'group', 'hierlemann');
        def.hima02 = fullfile('/net', 'bs-filesvr02', 'export', 'group', 'hierlemann');
        
        def.mea1kRoot = fullfile(def.hima02, 'recordings/Mea1k/');
        def.mea1kIntermediate = fullfile(def.hima02, 'intermediate_data/Mea1k/');
        def.mea1kData = fullfile(def.hima01, 'recordings/Mea1k/');
        def.mea1kAnalyzed = fullfile(def.hima01, 'AnalyzedData', 'Mea1k')
        def.Mea1kShared = fullfile(def.hima01, 'recordings', 'Mea1k', 'shared'); 
        
        def.tokenFiles = fullfile(def.mea1kIntermediate, user_name, 'tokens');
        def.serverData = fullfile(def.hima01, 'recordings/collaborations/');        
        
        def.eulerRoot = fullfile('/cluster','scratch', user_name);
        
        def.sortingOutPath = fullfile(def.hima01, 'recordings/HiDens/SpikeSorting/');
        def.serverSortingOut = fullfile(def.hima01, 'recordings/SpikeSortingOut/');
        def.serverDataRoska = fullfile(def.hima01, 'recordings/HiDens/Roska/');
        
        %% Felix:
        def.networkTempShare = fullfile('/links/groups/hierlemann/Temp/FelixFranke/');
        def.localData = fullfile('/net/bs-filesvr01/export/group/hierlemann/Temp/FelixFranke/LocalData/');
	def.localDataOnNetwork = fullfile(def.networkTempShare, 'LocalData');            
    
        %% Michele:
        def.micheleDSPaperDataPath = fullfile('/links/groups/hima/recordings/HiDens/Roska/Fiscella_2014/backup_results');
        def.svnDataRoska = fullfile('/home/frankef/bel.svn/hima_internal/cmosmea_recordings/trunk/Roska');        
        def.ravaDropBox = fullfile(def.localData, 'RavaDropboxLinux');
        def.ravaSimulations = fullfile(def.localData, 'RavaSimulations');
        
    end
    
    def.linuxSortingPath = '/net/bs-filesvr02/export/group/hierlemann';
        
