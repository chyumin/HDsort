function sortOneLEG(sortingName, wfsFile, covFile, legFolder, sortingParameters, groupidx, preprocessedFileList)
% INPUT:
% wfsFile:              Cut waveforms of all detected spikes (in the format
%                       that can be opened with hdsort.waveforms.WaveFormFile
%
% covFile:              File containing the noise covariance (.060cov.mat)
%
% legFolder:            Folder that contains the output files of this LEG
%                       (former: taskparameters.groupFolder)
% sortingName:          Name of this sorting (former: taskparameters.name)
%
% sortingParameters:    All parameters
%
% groupidx:             LEG.groupidx
%
% preprocessedFileList: List of preprocessed files for full template
%                       estimation

WFS = hdsort.file.WaveFormFile(wfsFile);

% Load noise covariance:
N = load(covFile);
%input.noise = N.noise;

%% Run sort function:
disp('Start sorting function...')
[S, P_] = hdsort.simplifiedSortScript(WFS, N.noise, legFolder, sortingName, sortingParameters);
disp('Sorting function finished.')

assert( all(S.clusteringMatched.ts > 0), '!!!')

%% Save gdf:
gdf = [S.clusteringMerged.ids, S.clusteringMatched.ts, S.clusteringMatched.amps, groupidx(S.clusteringMatched.ampsChIdx)];

gdfFile = fullfile(legFolder, [sortingName '.gdf.mat']);
save(gdfFile, 'gdf')

%% Template estimation:
templateFile    = fullfile(legFolder, [sortingName '_templates.mat']);
try
    load(templateFile);
    disp('Templates already estimated.')
catch
    disp('Starting template estimation...')
    DS              = hdsort.file.CMOSMEA(preprocessedFileList);
    MES             = DS.MultiElectrode;
    elGroupIndices  = groupidx;
    cutleft         = P_.templateEstimation.cutLeft;
    Tf              = P_.templateEstimation.Tf;
    maxSpikes       = P_.templateEstimation.maxSpikes;
    
    t_templateEstimation = tic;
    [wfs, nSourceSpikesPerTemplateAndChannel] = hdsort.waveforms.templateEstimation(...
        DS, gdf, 'Tf', Tf, 'cutLeft', cutleft, 'maxSpikes', maxSpikes);
    time = toc(t_templateEstimation);
    
    save(templateFile, 'wfs', 'cutleft', 'MES', 'elGroupIndices', 'nSourceSpikesPerTemplateAndChannel', 'time');
    disp('Template estimation finished.')
end

end