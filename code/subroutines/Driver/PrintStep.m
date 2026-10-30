function PrintStep(verbose,msg)
%PRINTSTEP Print a progress message followed by a separator line.
%
%   PrintStep(verbose, msg) prints msg if verbose is true.

if verbose
    disp(msg)
    disp(repmat('-',1,75))
end

end
