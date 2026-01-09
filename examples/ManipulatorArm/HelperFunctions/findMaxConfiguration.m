function [max_value, config, segment] = findMaxConfiguration(values)
    [max_value, i_lin] = max(abs(values(:)));
    [config, segment] = ind2sub(size(values), i_lin);
end