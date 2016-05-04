function [ avgRmse, stdRmse ] = fiComputeError( est, ref, type )

% [ avgRmse, stdRmse ] = fiComputeError( est, ref, type )
%
% Compute the mean and standard deviation of the RMSE error between the
% estimated and ground truth reference. If inputs are matrices then the
% RMSE is computed for every column, and the mean and standard deviation
% are estimated across columns. If inputs are cell arrays, then the RMSE is
% computed for every cell array entry, and means and standard deviation
% coputed across cell array entries.
%
% Inputs (required):
%    est - estimated quantity. A 2D (m x s) array or a cell array.
%    ref - reference quantity. A 2D (m x s) array or a cell array.
%    type - a string indicating if the quantities are normalized (i.e.
%      divided by their maximal values) before RMSE is computed. Allowed
%      values are {'','default','normalized'}.
%
% Outputs:
%    avgRMSE - the average RMSE across the s samples.
%    stdRMSE - the RMSE standard deviation across the s samples.
%
% Copytight, Henryk Blasinski 2016

p = inputParser;
p.addRequired('est',@(x) (isnumeric(x) || iscell(x)));
p.addRequired('ref',@(x) (isnumeric(x) || iscell(x)));
p.addOptional('type','default',@(x) validatestring(x,{'','default','normalized'});

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

