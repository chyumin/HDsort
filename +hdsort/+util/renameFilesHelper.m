function renameFilesHelper(x)
    % x is a full qualifier for a file
    [a,b,c] = fileparts(x);
    
    idx = strfind(b, '_proj_1');
    if ~isempty(idx)
        b(idx:idx+6) = [];
        xn = fullfile(a, [b c]);
        do_str = sprintf('!mv %s %s', x, xn);
        fprintf('%s\n', do_str);
        eval(do_str);
    end
end