X = magic(10);
X = X(:,1:5);

pd = pdefs();
ffile = fullfile(pd.localData, 'binaryFileMatrixTest.bin');
dims = size(X);
Y = mysort.ds.binaryFileMatrix(ffile, dims, 'writable', true, 'precision', 'single');
Y(:,:) = X;

K = mysort.ds.binaryFileMatrix(ffile, dims, 'writable', false, 'precision', 'single');

X
K(:,:)
delete(ffile)