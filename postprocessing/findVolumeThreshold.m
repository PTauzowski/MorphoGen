function t = findVolumeThreshold(x, targetVF)
% FINDVOLUMETHRESHOLD  Density threshold that preserves a target volume fraction.
%
%   t = findVolumeThreshold(x, targetVF)
%
%   Returns the largest t such that mean(x >= t) >= targetVF.
%   If targetVF > mean(x > 0), returns the minimum positive value in x.
%
%   Inputs:
%     x        [n x 1]  density field (any range)
%     targetVF scalar   target fraction of elements to keep, in (0, 1]
%
%   Output:
%     t  scalar threshold

    xs = sort(x(:), 'descend');
    n  = numel(xs);
    k  = min(n, max(1, ceil(targetVF * n)));
    t  = xs(k);
end
