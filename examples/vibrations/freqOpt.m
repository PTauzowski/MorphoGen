L=8; H=1; nelx=240; nely=30; volFrac=0.5;

%topFreqOptimization_MMA_bound( L,H,nelx,nely,volfrac,penal,rmin,maxiter,supportType,J)

[omega_c, xc] = topFreqOptimization_MMA_bound(L,H, nelx, nely, volFrac, 3.0, 1.5, 120, "CC", 3);

% [omega_c, xc] = topFreqOptimization_MMA(L,H,nelx,nely,volFrac,3.0,1.5,120,"CC");    % (c)
% [omega_a, xa] = topFreqOptimization_MMA(L,H,nelx,nely,volFrac,3.0,1.5,120,"HHmid"); % (a)
% [omega_b, xb] = topFreqOptimization_MMA(L,H,nelx,nely,volFrac,3.0,1.5,120,"CH");    % (b)

