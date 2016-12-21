function javaaddpath_silent(p)
    % adds p to the java path only if it is not already on the pathlist.
    % This would cause an annoying warning.
    A = javaclasspath();
    if ~any(strcmp(A, p))
        javaaddpath(p);
    end
    