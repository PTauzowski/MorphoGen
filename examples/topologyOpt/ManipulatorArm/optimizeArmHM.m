function xopt = optimizeArmHM(E,nu,segmentLength,R,r,res, alpha, lb, ub, x0, ShapeFn)

    of = @(x) -maxArmHM(E,nu,segmentLength,R,r,res, alpha, [0 x], ShapeFn);

    %xopt = rs(of, size(x0,2), 500000, [], [], [], [], lb, ub, []);
    %xopt = ga(of, size(x0,2), [], [], [], [], lb, ub, []);
    xopt = fmincon(of,x0,[],[],[],[],lb,ub,[]);

end

