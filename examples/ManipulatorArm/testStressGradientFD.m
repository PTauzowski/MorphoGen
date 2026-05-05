% testStressGradientFD
% Finite-difference validation of the adjoint stress gradient in
% solveSIMPVolumeStressMMA.  Builds a coarse single-segment solid model,
% applies a bending load, then checks dS/dx via centred FD.

clear; close all; clc;
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(42, 'twister');

%% ---- Build a coarse solid model ----------------------------------------
E  = 200e9;
nu = 0.3;
R  = 0.05;
r  = 0.03;
h  = 0.06;
alpha_deg = 5;
res_th    = 1;   % very coarse: fast FD

mdl      = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
analysis = mdl.analysis;

lc.My = E * mdl.Iy / R;
mdl.applyLoadCase(lc);
analysis.prepareRHSVectors();

nElems = analysis.getTotalElemsNumber();
fprintf('Mesh: %d elements, %d DOFs\n', nElems, analysis.getTaskDim());

%% ---- Problem parameters ------------------------------------------------
penal          = 3.0;
stressPNorm    = 8;
q              = 0.5;
nClusters      = 3;
xMinVal        = 0.1;

x    = xMinVal + (1 - xMinVal) * rand(nElems, 1);
xmin = xMinVal * ones(nElems, 1);
xmax = ones(nElems, 1);

%% ---- Compute adjoint gradient ------------------------------------------
fprintf('\nComputing adjoint gradient...\n');
[S0, dS0] = evalStress(analysis, x, penal, stressPNorm, q, nClusters);
fprintf('S0 = %s\n', mat2str(S0', 4));

%% ---- Centred finite-difference check -----------------------------------
h_fd = 1e-5;
candidateIds = find(xmax - xmin > 4*h_fd);
nTest = min(30, numel(candidateIds));
testIds = candidateIds(randperm(numel(candidateIds), nTest));

fprintf('\nFD check: %d elements, h=%.1e\n', nTest, h_fd);
fprintf('%-8s  %-8s  %14s  %14s  %10s\n', 'elemId', 'cluster', 'analytic', 'FD', 'relErr');
fprintf('%s\n', repmat('-', 1, 62));

allRelErrors = nan(nTest, nClusters);
for i = 1:nTest
    e    = testIds(i);
    hUse = min([h_fd, 0.4*(xmax(e)-x(e)), 0.4*(x(e)-xmin(e))]);
    if hUse <= 0, continue; end

    xp = x; xp(e) = x(e) + hUse;
    xm = x; xm(e) = x(e) - hUse;
    [Sp, ~] = evalStress(analysis, xp, penal, stressPNorm, q, nClusters);
    [Sm, ~] = evalStress(analysis, xm, penal, stressPNorm, q, nClusters);
    fdGrad  = (Sp - Sm) / (2*hUse);

    for c = 1:nClusters
        an  = dS0(e, c);
        fd  = fdGrad(c);
        rel = abs(an - fd) / max([abs(an), abs(fd), eps]);
        allRelErrors(i, c) = rel;
        if abs(an) > 1e-14 || abs(fd) > 1e-14
            fprintf('e=%5d  c=%d  an=% .4e  fd=% .4e  relErr=%.2e\n', e, c, an, fd, rel);
        end
    end
end

maxRel = max(allRelErrors(:), [], 'omitnan');
fprintf('\nMax relative error: %.3e\n', maxRel);
if maxRel < 0.01
    fprintf('PASS\n');
else
    fprintf('FAIL\n');
end

%% ======== Local helper functions ========================================

function [S, dS] = evalStress(analysis, x, penal, pNorm, q, nC)
    nElems = analysis.getTotalElemsNumber();
    x      = x(:);
    xPenal = x .^ penal;
    analysis.solveWeighted(xPenal);
    analysis.computeElementResults(xPenal);

    sigma   = gpHuberMises(analysis, nElems);
    relaxed = (max(x, eps) .^ q) .* sigma;
    target  = ones(nC, 1);

    clusters  = stressClusters(relaxed, nC);
    S         = zeros(nC, 1);
    dSdRelAll = zeros(nElems, nC);

    for c = 1:nC
        eIds  = clusters{c};
        ratio = relaxed(eIds) / target(c);
        mPow  = mean(ratio .^ pNorm);
        S(c)  = target(c) * mPow ^ (1/pNorm);
        if mPow <= eps, continue; end
        dSdr  = target(c) * (1/pNorm) * mPow^(1/pNorm-1) ...
            * pNorm * ratio.^(pNorm-1) / numel(eIds);
        dSdRelAll(eIds, c) = dSdr / target(c);
    end

    dS_direct = dSdRelAll .* ((penal + q) * relaxed ./ max(x, eps));

    P_adj  = adjLoad(analysis, x, penal, q, dSdRelAll, sigma);
    lam    = analysis.solveAdjointWithLoad(xPenal, P_adj);
    dS_ind = adjIndirect(analysis, lam, x, penal);

    dS = dS_direct + dS_ind;
end

function sigma = gpHuberMises(analysis, nElems)
    sigma      = zeros(nElems, 1);
    ei = analysis.getElemIndices();
    for i = 1:numel(analysis.felems)
        fe = analysis.felems{i};
        if ~isfield(fe.results,'gp') || ~isfield(fe.results.gp,'stress'), continue; end
        s  = fe.results.gp.stress;
        s1=s(:,:,1); s2=s(:,:,2); s3=s(:,:,3);
        s4=s(:,:,4); s5=s(:,:,5); s6=s(:,:,6);
        hm = sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2)+3*(s4.^2+s5.^2+s6.^2));
        sigma(ei{i}) = mean(hm, 2);
    end
end

function cl = stressClusters(s, nC)
    [~, ord] = sort(s(:), 'descend');
    n  = numel(ord);
    nC = min(max(1,nC), n);
    cl = cell(nC,1);
    edges = round(linspace(0, n, nC+1));
    for c = 1:nC
        ids = ord((edges(c)+1):edges(c+1));
        if isempty(ids), ids = ord(1); end
        cl{c} = ids;
    end
end

function P_adj = adjLoad(analysis, x, penal, q, dSdRelAll, ~)
    % Exact adjoint load: sigma_e = (1/nip)*sum_ip sHM(s_ip), so
    % d(sigma_e)/du = (1/nip)*sum_ip v_ip^T * x^p * D * B_ip (per-IP v).
    nTD  = analysis.getTaskDim();
    nC   = size(dSdRelAll, 2);
    ndN  = size(analysis.ndofs, 2);
    P_adj = zeros(nTD, nC);
    ei   = analysis.getElemIndices();
    nodes = analysis.mesh.nodes;
    for fi = 1:numel(analysis.felems)
        fe = analysis.felems{fi};
        if ~isa(fe, 'SolidElasticElem'), continue; end
        eIds = ei{fi};
        nEI  = numel(eIds);
        nnpE = size(fe.elems, 2);
        dimE = ndN * nnpE;
        s    = fe.results.gp.stress;   % (nEI, nip, 6)
        D    = fe.mat.D;
        intg = fe.sf.createIntegrator();
        dN   = fe.sf.computeGradient(intg.points);
        dNtr = permute(dN,[2,1,3]);
        nip  = size(intg.points, 1);
        fAll = zeros(dimE, nEI);
        cSz  = fe.assemblyChunkSize(dimE);
        for first = 1:cSz:nEI
            chunk = first:min(first+cSz-1, nEI);
            nc = numel(chunk);
            [Jinv, ~] = fe.jacobianInversePages(nodes, chunk, dNtr);
            B = fe.strainBPages(Jinv, dNtr);   % (6, dimE, nc, nip)
            for ip = 1:nip
                % Per-IP stress gradient direction v_ip for each element in chunk
                sip = reshape(s(chunk, ip, :), nc, 6);
                s1=sip(:,1); s2=sip(:,2); s3=sip(:,3);
                s4=sip(:,4); s5=sip(:,5); s6=sip(:,6);
                sHMip = sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2)+3*(s4.^2+s5.^2+s6.^2));
                vip = [(2*s1-s2-s3),(2*s2-s1-s3),(2*s3-s1-s2),6*s4,6*s5,6*s6] ...
                      ./ (2*max(sHMip, eps));   % (nc, 6)
                DvipP = reshape(D * vip', 6, 1, nc);   % (6, 1, nc)
                Bip = B(:,:,:,ip);   % (6, dimE, nc)
                BtDv = reshape(squeeze(pagemtimes(Bip,'transpose',DvipP,'none')), dimE, nc);
                fAll(:,chunk) = fAll(:,chunk) + BtDv / nip;
            end
        end
        nIPD = ceil((1:dimE)' / ndN);
        dIPD = mod((0:dimE-1)', ndN) + 1;
        aGD  = (fe.elems(:, nIPD') - 1)*ndN + dIPD';
        xpq  = x(eIds).^(penal+q);
        for c = 1:nC
            wc = dSdRelAll(eIds,c) .* xpq;
            if all(wc==0), continue; end
            wf = fAll .* wc';
            P_adj(:,c) = P_adj(:,c) + accumarray(aGD(:), reshape(wf', [], 1), [nTD,1]);
        end
    end
end

function dSi = adjIndirect(analysis, lam, x, penal)
    nE  = numel(x);
    nC  = size(lam, 2);
    dSi = zeros(nE, nC);
    ei  = analysis.getElemIndices();
    nodes = analysis.mesh.nodes;
    uN  = analysis.qnodal;
    sfn = 'computeStifnessMatrix';
    if analysis.isConst, sfn = 'computeStifnessMatrixConst'; end
    for fi = 1:numel(analysis.felems)
        fe = analysis.felems{fi};
        if ~isa(fe,'SolidElasticElem'), continue; end
        eIds = ei{fi};
        nEI  = numel(eIds);
        nnpE = size(fe.elems,2);
        dim  = nnpE * numel(fe.ndofs);
        K0   = reshape(fe.(sfn)(nodes, ones(nEI,1)), dim, dim, nEI);
        uE   = fe.createElemSolutionVectors(uN);
        K0u  = pagemtimes(K0, reshape(uE, dim, 1, nEI));
        pref = -penal * x(eIds).^(penal-1);
        for c = 1:nC
            lN  = analysis.fromFEMVector(lam(:,c));
            lE  = fe.createElemSolutionVectors(lN);
            bil = squeeze(pagemtimes(reshape(lE,1,dim,nEI), K0u));
            dSi(eIds,c) = dSi(eIds,c) + pref .* bil(:);
        end
    end
end
