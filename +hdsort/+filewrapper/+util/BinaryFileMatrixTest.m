X = magic(10);
X = X(:,1:5);

pd = pdefs();
ffile = fullfile(pd.localData, 'binaryFileMatrixTest.bin');
dims = size(X);
Y = hdsort.filewrapper.util.binaryFileMatrix(ffile, dims, 'writable', true, 'precision', 'single');
Y(:,:) = X;

K = hdsort.filewrapper.util.binaryFileMatrix(ffile, dims, 'writable', false, 'precision', 'single');

X
K(:,:)
delete(ffile)