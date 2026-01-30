function eps_val = continuation_eps(iter, ep_old, rank, mid_num, ini)
%CONTINUATION_EPS  Appendix-style epsilon schedule used in MMC codes.
%   Mirrors the 2D refactored schedule; returns scalar epsilon.
    eps_val = ini + ep_old - (1/ep_old + exp(-rank*(iter - mid_num)))^(-1);
end

function H = heaviside_smooth(phi, alpha, epsilon)
%HEAVISIDE_SMOOTH Smoothed Heaviside used for density interpolation.
    H = 3*(1-alpha)/4 * (phi/epsilon - phi.^3/(3*(epsilon)^3)) + (1+alpha)/2;
    H(phi >  epsilon) = 1;
    H(phi < -epsilon) = alpha;
end

function [eigVec, eigVal] = eigs_wrapper(K, M, dof, n_modes)
%EIGS_WRAPPER Thin wrapper to keep eigen solve identical across 2D/3D.
    [eigVec, eigVal] = eigs(K(dof, dof), M(dof, dof), n_modes, 'sm');
end

