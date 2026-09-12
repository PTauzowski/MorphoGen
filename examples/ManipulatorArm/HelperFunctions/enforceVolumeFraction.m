function x = enforceVolumeFraction(x, VolFrac, xmin, xmax)
    x = min(max(x(:), xmin), xmax);
    target = VolFrac * numel(x);
    target = min(max(target, sum(xmin)), sum(xmax));

    lo = min(xmin - x);
    hi = max(xmax - x);
    for k = 1:80
        shift = 0.5 * (lo + hi);
        candidate = min(max(x + shift, xmin), xmax);
        if sum(candidate) < target
            lo = shift;
        else
            hi = shift;
        end
    end
    x = min(max(x + 0.5 * (lo + hi), xmin), xmax);
end
