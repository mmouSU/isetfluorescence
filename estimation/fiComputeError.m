function [ avgRmse, stdRmse ] = fiComputeError( est, ref, type )

p = inputParser;
p.addRequired('est',@(x) (isnumeric(x) || iscell(x)));
p.addRequired('ref',@(x) (isnumeric(x) || iscell(x)));
p.addOptional('type','default',@ischar);

p.parse(est,ref,type);
inputs = p.Results;


switch inputs.type

    case 'normalized'
        % Normalize each measurement by the maximal value.
        if iscell(est)
            rmse = cellfun(@(x,y) sqrt(mean((x(:)/max(x(:)) - y(:)/max(y(:))).^2)),est,ref);
        else
            rmse = sqrt(mean((est*diag(1./max(est)) - ref*diag(1./max(ref))).^2));
        end

    otherwise
        % Make a comparison as is.
        if iscell(est)
            rmse = cellfun(@(x,y) sqrt(mean((x(:) - y(:)).^2)),est,ref);
        else
            rmse = sqrt(mean((est - ref).^2));
        end
end

avgRmse = mean(rmse);
stdRmse = std(rmse);


end

