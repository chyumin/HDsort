%%
LEGs = {[1 2 3], [3 4 5], [5], [6 7 8]};

%%
X1 = randn(210000,10);
X1(10000:10003,3) = -10;
X1(10100:10103,3) = -10;
X2 = randn(220000,10);
X2(2000:2003,3) = -10;
X2(2100:2003,3) = -10;
DS1 = hdsort.file.DataMatrix(X1, 'samplesPerSecond', 20000);
DS2 = hdsort.file.DataMatrix(X2, 'samplesPerSecond', 20000);

%%
outFolder = 'test01';
jobName = 'J180214Bconf2';
DSlist = DS1;
fnlist = {'DS1.1'};
prep = hdsort.Preprocessor(DSlist, LEGs, jobName, ...
    'prefilter', 1, 'chunkSize', 10000, 'minChunkSize', 1000,...
    'saveRawH5FileNameList', fnlist, 'debug', 1, 'parfor', 1, ...
    'forceFileDeletionIfExists', 1);

preprocessorFile = fullfile(outFolder, [jobName '_preprocessor.mat'])

R = prep.preprocess(outFolder, preprocessorFile);

%%
F = hdsort.file.CMOSMEA(R.preprocessedFiles)

%%
figure; histogram(DS1(:, 1), 30)

%%
figure; histogram(F(:, 1), 30)
