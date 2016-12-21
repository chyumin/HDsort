function c = config()

c.submit_host = 'bs-submit01';
c.home = fullfile('/net/bs-filesvr02/export/group/hierlemann/');
%c.logFolder = fullfile('/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/');
%c.home = fullfile('/links/grid/scratch/hierlemann_sortings');

if ismac
    pd = pdefs();
    %c.home = fullfile('/Users', getenv('USER'), 'tmp' )
    c.home = fullfile('/Volumes', 'hierlemann');
    %c.logFolder = fullfile('/Volumes', 'hierlemann', 'intermediate_data', 'Mea1k');
end

end


