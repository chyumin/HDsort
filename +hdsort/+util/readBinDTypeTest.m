% Just write random stuff in a binary file

% C = {'hallo', pi, 'nonsense', 2*pi, 1};
% PRECS = {'uchar', 'float32', 'uchar', 'float64', 'uint16'};
% REPET = [5, 1, 8, 1, 1];
% fh = fopen('file.bin','w');
% for i=1:10%0000*2*100
%    for k=1:length(REPET)
%        fwrite(fh, C{1,k}, PRECS{k});
%    end
% end
% fclose(fh);



% Define the DType and read the binary file

DType = {'name', '5*uint8=>char',   5
        'pi',   'float32=>single', 4
        '',     'skip',            8
        'pipi', 'float64=>double', 8
        'one',  'uint16=>int',     2};
tic
Y = hdsort.util.readBinDtype('C:\\Data\\file.bin', DType);
toc