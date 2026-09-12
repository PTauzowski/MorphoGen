function [candIters, bestIter] = selectOptimalIteration(volHistory, VolFrac, window)
% SELECTOPTIMALITERATION  Find candidate iterations near a target volume fraction.
%
%   [candIters, bestIter] = selectOptimalIteration(volHistory, VolFrac, window)
%
%   Searches the ESO volume-fraction history for iterations whose volume
%   fraction lies within [VolFrac - window, VolFrac + window].  Returns the
%   indices of all such iterations and a single "best" recommendation.
%
%   The best iteration is the last one whose volume fraction is at or below
%   VolFrac + eps — i.e. the iteration where the ESO first achieves the
%   target.  If no iteration is that close, the one with minimum |VF - VolFrac|
%   is returned.
%
%   Inputs:
%     volHistory  [nIter x 1]  volume fraction at each ESO iteration
%     VolFrac     scalar       target volume fraction
%     window      scalar       half-width of search window (default: 0.05)
%
%   Outputs:
%     candIters   [k x 1]  iteration indices within the window (ascending)
%     bestIter    scalar   recommended single iteration

    if nargin < 3 || isempty(window), window = 0.05; end

    volHistory = volHistory(:);
    nIter      = numel(volHistory);

    % Candidate iterations within [VolFrac - window, VolFrac + window]
    inWindow  = abs(volHistory - VolFrac) <= window;
    candIters = find(inWindow);

    if isempty(candIters)
        % Fallback: closest iteration to target
        [~, bestIter] = min(abs(volHistory - VolFrac));
        candIters = bestIter;
        return;
    end

    % Among candidates, prefer the last iteration at or below VolFrac
    atOrBelow = candIters(volHistory(candIters) <= VolFrac + 1e-6);
    if ~isempty(atOrBelow)
        bestIter = atOrBelow(end);
    else
        % All candidates are above target — take the one closest from above
        [~, minIdx] = min(volHistory(candIters) - VolFrac);
        bestIter = candIters(minIdx);
    end
end
