dpath = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/frankef/configs_7Jul2015_CNO_michele/SortOut/Sorting1/Config1/groups/';

filefilter = 'J*_proj_1*.bin';
fhandle = @hdsort.util.renameFilesHelper;

hdsort.util.runFunctionOnFilesRecursively(dpath, filefilter, fhandle)

