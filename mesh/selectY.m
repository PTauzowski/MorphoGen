function selector = selectY(y0)
    selector = Selector(@(x) x(:,2) - y0);
end

