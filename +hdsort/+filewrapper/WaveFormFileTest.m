if ismac
    wfsFile = '/Users/rolandd/tmp/testwfsfile.h5'
else
    wfsFile = '/home/rolandd/tmp/testwfsfile.h5'
end

cutLeft = 4
Tf = 10
nC = 3

%%
delete(wfsFile);
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'cutLeft', cutLeft, 'Tf', Tf, 'nC', nC);
x = W(:,:);
y = W(:,:,:);
xgdf = W.getGdf();
clear W
assert(isempty(x), '!')
assert(isempty(y), '!')
assert(isempty(xgdf), '!')


%%
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'writable', false);
x = W(:,:);
y = W(:,:,:);
clear W
assert(isempty(x), '!')
assert(isempty(y), '!')

%%
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'writable', true);
x = W(:,:);
y = W(:,:,:);
assert(isempty(x), '!')
assert(isempty(y), '!')

%% Make sure that the initialized WaveFormFile is empty:
wGDF = W.getGdf();
vWFS = W(:,:);
tWFS = W(:,:, :);

assert( isempty(wGDF), '!')
assert( isempty(vWFS), '!')
assert( isempty(tWFS), '!')

%% Also after reopening:
clear W
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'writable', true);
wGDF = W.getGdf();
vWFS = W(:,:);
tWFS = W(:,:, :);

assert( isempty(wGDF), '!')
assert( isempty(vWFS), '!')
assert( isempty(tWFS), '!')


%% Add some waveforms:
N = 10;
gdf = [repmat(1, 1, N), (1:N); repmat(2, 1, N), (1:N)+20]';
wfs = 100*randn(N*2, Tf*nC); % waveforms should be double and somewhere between ~350 and 350
W.addWaveforms(wfs, gdf);

WFS = wfs;
GDF = gdf;

assert(isequal(size(wfs), size(W)), '!')


%% Add more waveforms:
clear W
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'writable', true);
% N=1
gdf = [repmat(1, 1, N), (1:N)+30; repmat(2, 1, N), (1:N)+40]';
wfs = 100*randn(N*2, Tf*nC);
W.addWaveforms(wfs, gdf);

WFS = [WFS; wfs];
GDF = [GDF; gdf];

tWFS = hdsort.waveforms.v2t(WFS, nC);

clear W;

%%
W = hdsort.filewrapper.WaveFormFile(wfsFile);
Wgdf = W.getGdf();
assert(isequal( Wgdf(:, 1:2), GDF), '!')

reducedGDF = W.getGdf([1:4, 7]);
assert(isequal( reducedGDF(:, 1:2), GDF([1:4, 7],1:2)), '!')

%% Get full waveforms back:
tW = W(:,:, :);
assert(isequal(tW, tWFS), '!')

%% Get part of the waveforms back (1):
wfidx = [1, 13, 31];
tW = W(:,:,wfidx);
assert(isequal(tW, tWFS(:, :,wfidx)), '!')

%% Get part of the waveforms back (2):
chidx = [1, 3];
tW = W(:,chidx,:);
assert(isequal(tW, tWFS(:,chidx,:)), '!')

%% Get part of the waveforms back (3):
tfidx = [4:8];
tW = W(tfidx,:,:);
assert(isequal(tW, tWFS(tfidx,:,:)), '!')

%% Get part of the waveforms back (4):
tW = W(tfidx,chidx,wfidx);
assert(isequal(tW, tWFS(tfidx,chidx,wfidx)), '!')

%% Get part of the waveforms back (5):
tfchidx = 1:15;
vW = W(wfidx, tfchidx);
assert(isequal(vW, WFS(wfidx, tfchidx)), '!')

clear W


%% === Append data to the file ===
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'writable', true);

N = 14;
gdf = [repmat(1, 1, N), (1:N)+130; repmat(2, 1, N), (1:N)+140]';
wfs = 100*randn(N*2, Tf*nC);
W.addWaveforms(wfs, gdf);

WFS = [WFS; wfs];
tWFS = hdsort.waveforms.v2t(WFS, nC);

clear W;
clear wfs;

%% Check size function:
W = hdsort.filewrapper.WaveFormFile(wfsFile)
s = size(W)
clear W

%% === Append data to the file ===
W = hdsort.filewrapper.WaveFormFile(wfsFile, 'writable', true);

N = 1;
gdf = [repmat(1, 1, N)   (1:N)+130]';
wfs = 100*randn(N*1, Tf*nC);
W.addWaveforms(wfs);

WFS = [WFS; wfs];
tWFS = hdsort.waveforms.v2t(WFS, nC);

clear W;
clear wfs;
%% Get part of the waveforms back (6):
wfidx = [1, 14, 31, 67];
W = hdsort.filewrapper.WaveFormFile(wfsFile);
tW = W(tfidx,chidx,wfidx);
assert(isequal(tW, tWFS(tfidx,chidx,wfidx)), '!')


%% Create new file and write directly into it:
wfsFile2 = '/Users/rolandd/tmp/testwfsfile2.h5'
delete(wfsFile2);
W = hdsort.filewrapper.WaveFormFile(wfsFile2, 'cutLeft', cutLeft, 'Tf', Tf, 'nC', nC);

N = 10;
gdf = [repmat(1, 1, N), (1:N); repmat(2, 1, N), (1:N)+20]';
wfs = 100*randn(Tf, nC, N*2); % waveforms should be double and somewhere between ~350 and 350
W.addWaveforms(wfs, gdf);

tWFS = wfs;
GDF = gdf;

tW = W(:,:, :);
assert(isequal(tW, tWFS), '!')
Wgdf = W.getGdf();
assert(isequal( Wgdf(:, 1:2), GDF), '!')

clear W
W2 = hdsort.filewrapper.WaveFormFile(wfsFile2, 'writable', false);
tW2 = W2(:,:, :);
assert(isequal(tW2, tWFS), '!')
clear W2

W3 = hdsort.filewrapper.WaveFormFile(wfsFile2, 'writable', true);
W3.addWaveforms(wfs, gdf);
tWFS3 = cat(3, tWFS, wfs);
tW3 = W3(:,:, :);
assert(isequal(tW3, tWFS3), '!')

%% Try it in a parfor:
cutLeft = 4
Tf = 10
nC = 3
WF = {}; wfsFile = {};
for ii = 1:10
    if ismac
        wfsFile{ii} = ['/Users/rolandd/tmp/testwfs/file' num2str(ii) '.h5']
    else
        wfsFile{ii} = ['/home/rolandd/tmp/testwfs/file' num2str(ii) '.h5']
    end
    delete(wfsFile{ii});
    WF{ii} = hdsort.filewrapper.WaveFormFile(wfsFile{ii}, 'cutLeft', cutLeft, 'Tf', Tf, 'nC', nC+ii);
end
parfor ii = 1:10
    % This should not work!!!!! However, on MAC it does!
    N = 10;
    gdf = [repmat(1, 1, N), (1:N); repmat(2, 1, N), (1:N)+20]';
    wfs = 100*randn(N*2, Tf*(nC+ii));
    WF{ii}.addWaveforms(wfs, gdf);
end 
clear WF
for ii = 1:10
    WF{ii} = hdsort.filewrapper.WaveFormFile(wfsFile{ii});
    size(WF{ii})
end
%%
clear WF
parfor ii = 1:10
    % This should not work!!!!! However, on MAC it does!
    N = 10;
    gdf = [repmat(1, 1, N), (1:N); repmat(2, 1, N), (1:N)+20]';
    wfs = 100*randn(N*2, Tf*(nC+ii));
    WF = hdsort.filewrapper.WaveFormFile(wfsFile{ii}, 'writable', true);
    WF.addWaveforms(wfs, gdf);
    WF = [];
end 



