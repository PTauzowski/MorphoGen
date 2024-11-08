function selector = selectZ(z0)
    selector = Selector(@(x) x(:,3) - z0);
end

