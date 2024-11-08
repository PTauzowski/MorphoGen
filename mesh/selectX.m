function selector = selectX(x0)
    selector = Selector(@(x) x(:,1) - x0);
end

