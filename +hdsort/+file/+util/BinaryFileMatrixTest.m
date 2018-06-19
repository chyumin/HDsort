X = magic(10);
X = X(:,1:5);

pd = pdefs();
ffile = fullfile(pd.localData, 'BinaryFileMatrixTest.bin');
dims = size(X);
Y = hdsort.file.util.BinaryFileMatrix(ffile, dims, 'writable', true, 'precision', 'single');
Y(:,:) = X;

K = hdsort.file.util.BinaryFileMatrix(ffile, dims, 'writable', false, 'precision', 'single');

X
K(:,:)
delete(ffile)