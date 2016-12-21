function [dS iseq] = compareStructs(struct1, struct2)

fields = sort(fieldnames(struct1));
fields2 = sort(fieldnames(struct2));

differentFields= {};
differentFieldValues = {};

if any(~cellfun(@strcmp, fields, fields2))
    iseq = false;
    dS = struct();
    return;
end
assert( ~any(~cellfun(@strcmp, fields, fields2)), 'Structures that are not of the same type can not be compared!')

dS1 = struct();
dS2 = struct();

n = 1;
for f_ = fields'
    fn = f_{1};
    
    iseq = false;
    if numel([{struct1.(fn)}]) == 1
        if isnumeric(struct1.(fn))
            iseq = isequal(struct1.(fn), struct2.(fn));
        elseif ischar(struct1.(fn))
            iseq = strcmp(struct1.(fn), struct2.(fn));
        elseif iscell(struct1.(fn)) && isnumeric(struct1.(fn){1})
            iseq =  ~any(~cellfun(@isequal, struct1.(fn), struct2.(fn) ));
        elseif iscell(struct1.(fn)) && ischar(struct1.(fn){1})
            iseq = numel(struct1.(fn)) == numel(struct2.(fn));
            if iseq
                iseq =  ~any(~cellfun(@strcmp, struct1.(fn), struct2.(fn) ));
            end
        elseif isempty(struct1.(fn))
            iseq = isempty(struct2.(fn));
        %elseif strcmp( class(struct1.(fn)), class(struct2.(fn)))
        %    [dS_ iseq] = hdsort.util.compareStructs(struct1.(fn), struct2.(fn));
        end
    end
    
    if ~iseq
        
        %differentFields{n} = fn;
        %differentFieldValues{n, 1} = struct1.(fn);
        %differentFieldValues{n, 2} = struct2.(fn);
        if numel([{struct1.(fn)}]) == 1
            dS1.(fn) = struct1.(fn);
            dS2.(fn) = struct2.(fn);
        else
            dS1.(fn) = [{struct1.(fn)}];
            dS2.(fn) = [{struct2.(fn)}];
        end
         
        n = n + 1;
    end
    
end

dS = struct(dS1);
dS(2) = dS2;
iseq = n == 1;


