function tagFilesWithHeader(dpath, filetype, templatefile, template_content)
    % recursively walks through dpath and adds the content of templatefile 
    % to every file with the ending filetype.
    
    if nargin < 4
        % We have to load the templatefile once
        fh = fopen(templatefile);
        template_content =  [];
        tline = fgets(fh);
        while ischar(tline)
            template_content = [template_content tline];
            tline = fgets(fh);
        end
        fclose(fh);
    end
    
    fprintf('Processing directory: %s\n', dpath);
    fnames = dir(dpath);
    for i=1:length(fnames)
        myFile = fnames(i);
        if myFile.name(1) == '.'
            continue
        end
           
        if myFile.isdir
            hdsort.util.tagFilesWithHeader([dpath filesep myFile.name], filetype, templatefile, template_content);
        elseif strcmp(myFile.name(end-length(filetype)+1:end),filetype)
            fprintf('    %s\n', myFile.name);
            tagfile([dpath filesep myFile.name], template_content);
        end
    end
    
    function tagfile(file, header_str)
        header_end = '%@_________________________________________________________________________';
        fh = fopen(file, 'r');
        tline = fgets(fh);
        still_in_header = 1;
        content = [];
        while ischar(tline) 
            if  isempty(tline) || tline(1) ~= '%'
                still_in_header = 0;
            end
            if still_in_header == 0
                content = [content tline];
            else
                k = strfind(tline, header_end);
                if ~isempty(k)
                    still_in_header=0;
                end
            end
            tline = fgets(fh);
        end
        fclose(fh);
        
        content = [header_str sprintf('\r\n') content];
        
        fh = fopen(file, 'w');
        fwrite(fh, content);
        fclose(fh);
    end   
end