function rawFile = download_testfile(rawFile)
if exist(rawFile, 'file') == 0
    userInput = input(['Do you want to download test file to ' mainFolder ' (~2GB)? Y/N:'],'s');
    assert(strcmp(userInput, 'Y'), 'User denied');
    mkdir(mainFolder);
    url = 'https://polybox.ethz.ch/index.php/s/4JqS6b0q2VGZKs6/download';
    websave(rawFile,url)
else
    disp(['Test file ' rawFile ' already downloaded.'])
end
end