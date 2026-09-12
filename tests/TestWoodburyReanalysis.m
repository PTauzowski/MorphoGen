classdef TestWoodburyReanalysis < matlab.unittest.TestCase
    properties (Constant)
        ElemDofs = 24
        KeVecLen = 24 * 24
    end

    methods (TestClassSetup)
        function addMathPath(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(projectRoot, 'math')));
        end
    end

    methods (TestMethodSetup)
        function setDeterministicRng(~)
            rng(1, 'twister');
        end
    end

    methods (Test)
        function testDimensionsAndDiagonalC(testCase)
            nFree = 200;
            nFixed = 40;
            nChanged = 20;
            eigTol = 1.0e-14;

            [supports, freedofs, ~, freeMap] = createSupportAndMaps(nFree, nFixed);
            edofs_list = sampleElementDofs(freedofs, nChanged, testCase.ElemDofs);
            Ke0_vec_list = randomSPDKeVectors(nChanged);

            w_old = linspace(0.30, 0.78, nChanged)';
            dw = 0.03 * sin((1:nChanged)' * 0.7);
            w_new = w_old + dw;

            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, eigTol);

            testCase.verifySize(U, [nFree, size(U, 2)]);
            testCase.verifyEqual(size(C, 1), size(U, 2));
            testCase.verifyEqual(size(C, 2), size(U, 2));
            testCase.verifyEqual(nnz(C - diag(diag(C))), 0);
            testCase.verifyTrue(all(diag(C) == 1 | diag(C) == -1));
            testCase.verifyEqual(info.rank, size(U, 2));
            testCase.verifyEqual(info.rank, sum(info.keptPerElem));

            solver = LinearEquationsSystem(1, 1, supports);
            freeMapFromSolver = solver.createFreeMap();
            testCase.verifyEqual(freeMapFromSolver, freeMap);

            [Ue, Ce, infoe] = lowRankFromElementUpdates_vecKe_array( ...
                zeros(0, testCase.ElemDofs), zeros(0, testCase.KeVecLen), ...
                zeros(0, 1), zeros(0, 1), freeMap, nFree, eigTol);
            testCase.verifySize(Ue, [nFree, 0]);
            testCase.verifySize(Ce, [0, 0]);
            testCase.verifyEqual(infoe.rank, 0);
            testCase.verifyEqual(infoe.nChanged, 0);
            testCase.verifySize(infoe.keptPerElem, [0, 1]);
        end

        function testWoodburyMatchesDirectSmall(testCase)
            nFree = 500;
            nFixed = 60;
            nChanged = 4;

            [Kbase, decBase, b] = buildSPDBaseSystem(nFree, 0.02, 10.0);
            [supports, freedofs, fixeddofs, freeMap] = createSupportAndMaps(nFree, nFixed);

            edofs_list = sampleElementDofs(freedofs, nChanged, testCase.ElemDofs);
            Ke0_vec_list = randomSPDKeVectors(nChanged);

            w_old = [0.30; 0.40; 0.50; 0.60];
            w_new = w_old + [0.06; -0.04; 0.03; -0.02];

            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-14);
            testCase.verifyGreaterThan(info.rank, 0);

            dK = U * C * U';
            Knew = Kbase + dK;

            x_direct = Knew \ b;
            x_helper = woodburySolve(decBase, b, U, C);
            relerrHelper = relativeError(x_direct, x_helper);
            testCase.verifyLessThan(relerrHelper, 1.0e-9);

            solver = LinearEquationsSystem(1, 1, supports);
            x_class_free = solver.solveWoodburyFree(decBase, b, U, C);
            relerrClassFree = relativeError(x_direct, x_class_free);
            testCase.verifyLessThan(relerrClassFree, 1.0e-9);

            Pfull = zeros(nFree + nFixed, 1);
            Pfull(freedofs) = b;
            q_full = solver.solveWoodbury(decBase, Pfull, U, C);
            x_class_full = q_full(freedofs);
            relerrClassFull = relativeError(x_direct, x_class_full);
            testCase.verifyLessThan(relerrClassFull, 1.0e-9);
            testCase.verifyEqual(norm(q_full(fixeddofs)), 0.0, 'AbsTol', 1.0e-14);
        end

        function testNoUpdatesDwZero(testCase)
            nFree = 320;
            nFixed = 30;
            nChanged = 5;

            [Kbase, decBase, b] = buildSPDBaseSystem(nFree, 0.02, 8.0);
            [supports, freedofs, ~, freeMap] = createSupportAndMaps(nFree, nFixed);

            edofs_list = sampleElementDofs(freedofs, nChanged, testCase.ElemDofs);
            Ke0_vec_list = randomSPDKeVectors(nChanged);

            w_old = linspace(0.2, 0.8, nChanged)';
            w_new = w_old;

            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-14);

            testCase.verifySize(U, [nFree, 0]);
            testCase.verifySize(C, [0, 0]);
            testCase.verifyEqual(info.rank, 0);

            x_base = decBase \ b;
            x_helper = woodburySolve(decBase, b, U, C);
            testCase.verifyLessThan(relativeError(x_base, x_helper), 1.0e-12);

            solver = LinearEquationsSystem(1, 1, supports);
            x_class = solver.solveWoodburyFree(decBase, b, U, C);
            testCase.verifyLessThan(relativeError(x_base, x_class), 1.0e-12);
        end

        function testAllElementDofsConstrained(testCase)
            nFree = 180;
            nFixed = 100;
            nChanged = 4;

            [~, decBase, b] = buildSPDBaseSystem(nFree, 0.02, 8.0);
            [supports, ~, fixeddofs, freeMap] = createSupportAndMaps(nFree, nFixed);

            edofs_list = sampleElementDofs(fixeddofs, nChanged, testCase.ElemDofs);
            Ke0_vec_list = randomSPDKeVectors(nChanged);
            w_old = 0.4 * ones(nChanged, 1);
            w_new = w_old + [0.10; -0.05; 0.03; -0.02];

            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-14);

            testCase.verifySize(U, [nFree, 0]);
            testCase.verifySize(C, [0, 0]);
            testCase.verifyEqual(info.rank, 0);

            x_base = decBase \ b;
            x_helper = woodburySolve(decBase, b, U, C);
            testCase.verifyLessThan(relativeError(x_base, x_helper), 1.0e-12);

            solver = LinearEquationsSystem(1, 1, supports);
            x_class = solver.solveWoodburyFree(decBase, b, U, C);
            testCase.verifyLessThan(relativeError(x_base, x_class), 1.0e-12);
        end

        function testNegativeUpdatesStillWorkSPD(testCase)
            nFree = 300;
            nFixed = 50;
            nChanged = 20;

            [Kbase, decBase, b] = buildSPDBaseSystem(nFree, 0.03, 120.0);
            [~, freedofs, ~, freeMap] = createSupportAndMaps(nFree, nFixed);
            edofs_list = sampleElementDofs(freedofs, nChanged, testCase.ElemDofs);
            Ke0_vec_list = randomSPDKeVectors(nChanged);

            w_old = 0.6 * ones(nChanged, 1);
            w_new = w_old - linspace(1.0e-4, 5.0e-4, nChanged)';

            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-14);
            testCase.verifyGreaterThanOrEqual(info.rank, 1);
            testCase.verifyTrue(any(diag(C) == -1));

            Knew = full(Kbase + U * C * U');
            Knew = 0.5 * (Knew + Knew');
            [~, cholFlag] = chol(Knew);
            testCase.verifyEqual(cholFlag, 0);

            x_direct = Knew \ b;
            x_helper = woodburySolve(decBase, b, U, C);
            testCase.verifyLessThan(relativeError(x_direct, x_helper), 1.0e-9);
        end

        function testRankTruncationEffect(testCase)
            nFree = 400;
            nFixed = 40;
            nChanged = 20;

            [Kbase, decBase, b] = buildSPDBaseSystem(nFree, 0.015, 2.0);
            [~, freedofs, ~, freeMap] = createSupportAndMaps(nFree, nFixed);
            edofs_list = sampleElementDofs(freedofs, nChanged, testCase.ElemDofs);
            Ke0_vec_list = spectrumControlledKeVectors(nChanged, logspace(0, -10, 24));

            w_old = zeros(nChanged, 1);
            w_new = linspace(3.0, 5.0, nChanged)';

            [Ufull, Cfull, infoFull] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-14);
            [Utr, Ctr, infoTr] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-2);

            testCase.verifyGreaterThan(infoFull.rank, infoTr.rank);

            KnewFull = Kbase + Ufull * Cfull * Ufull';
            x_direct = KnewFull \ b;
            x_approx = woodburySolve(decBase, b, Utr, Ctr);

            relerr = relativeError(x_direct, x_approx);
            testCase.verifyGreaterThan(relerr, 1.0e-6);
        end

        function testPerformanceSanityOptional(testCase)
            nFree = 2000;
            nChanged = 2;

            [Kbase, decBase, b] = buildSPDBaseSystem(nFree, 0.003, 6.0);
            freeMap = (1:nFree)';
            edofs_list = sampleElementDofs((1:nFree)', nChanged, testCase.ElemDofs);
            Ke0_vec_list = randomSPDKeVectors(nChanged);

            w_old = [0.4; 0.6];
            w_new = w_old + [0.02; -0.015];

            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, 1.0e-14);

            dK = U * C * U';
            Knew = Kbase + dK;

            tDirect = tic;
            x_direct = Knew \ b;
            timeDirect = toc(tDirect);

            tWoodbury = tic;
            x_woodbury = woodburySolve(decBase, b, U, C);
            timeWoodbury = toc(tWoodbury);

            relerr = relativeError(x_direct, x_woodbury);
            fprintf('[performanceSanity] nFree=%d rank=%d direct=%.4fs woodbury=%.4fs relerr=%.3e\n', ...
                nFree, info.rank, timeDirect, timeWoodbury, relerr);

            testCase.verifyLessThan(relerr, 1.0e-8);
        end

        function testStabilityNoInvAndInformativeErrors(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            linearSystemSrc = fileread(fullfile(projectRoot, 'math', 'LinearEquationsSystem.m'));
            lowRankSrc = fileread(fullfile(projectRoot, 'math', 'lowRankFromElementUpdates_vecKe_array.m'));

            testCase.verifyFalse(contains(linearSystemSrc, 'inv('));
            testCase.verifyFalse(contains(lowRankSrc, 'inv('));

            solver = LinearEquationsSystem(1, 1, false(20, 1));
            decBase = decomposition(speye(20), 'chol');
            b = randn(20, 1);
            U = randn(20, 3);
            Cbad = eye(2);

            try
                solver.solveWoodburyFree(decBase, b, U, Cbad);
                testCase.verifyFail('Expected an error for mismatched C size.');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'C must be an (m x m) matrix'), ...
                    sprintf('Unexpected error message: %s', ME.message));
            end

            try
                solver.solveWoodburyFree(randn(20, 1), zeros(20, 0), zeros(0, 0));
                testCase.verifyFail('Expected an error when base decomposition is not set.');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'Base matrix not set'), ...
                    sprintf('Unexpected error message: %s', ME.message));
            end
        end
    end
end

function [supports, freedofs, fixeddofs, freeMap] = createSupportAndMaps(nFree, nFixed)
    nGlobal = nFree + nFixed;
    dofs = (1:nGlobal)';
    shuffled = dofs(randperm(nGlobal));

    freedofs = sort(shuffled(1:nFree));
    fixeddofs = setdiff(dofs, freedofs, 'stable');

    supports = false(nGlobal, 1);
    supports(fixeddofs) = true;

    freeMap = zeros(nGlobal, 1);
    freeMap(freedofs) = (1:nFree)';
end

function [Kbase, decBase, b] = buildSPDBaseSystem(nFree, density, alpha)
    A = sprandn(nFree, nFree, density);
    Kbase = A' * A + alpha * speye(nFree);
    Kbase = 0.5 * (Kbase + Kbase');
    decBase = decomposition(Kbase, 'chol');
    b = randn(nFree, 1);
end

function edofs_list = sampleElementDofs(pool, nChanged, elemDofs)
    assert(numel(pool) >= elemDofs, 'Not enough DOFs in the sampling pool.');
    edofs_list = zeros(nChanged, elemDofs);
    for i = 1:nChanged
        pick = randperm(numel(pool), elemDofs);
        edofs_list(i, :) = pool(pick);
    end
end

function Ke0_vec_list = randomSPDKeVectors(nChanged)
    Ke0_vec_list = zeros(nChanged, 24 * 24);
    for i = 1:nChanged
        B = randn(24, 24);
        Ke = B' * B + 1.0e-3 * eye(24);
        Ke0_vec_list(i, :) = Ke(:)';
    end
end

function Ke0_vec_list = spectrumControlledKeVectors(nChanged, spectrum)
    Ke0_vec_list = zeros(nChanged, 24 * 24);
    for i = 1:nChanged
        [Q, ~] = qr(randn(24, 24), 0);
        Ke = Q * diag(spectrum) * Q';
        Ke = 0.5 * (Ke + Ke');
        Ke0_vec_list(i, :) = Ke(:)';
    end
end

function x = woodburySolve(decBase, b, U, C)
    y = decBase \ b;
    if isempty(U)
        x = y;
        return;
    end
    AU = decBase \ U;
    cdiag = diag(C);
    Sinv = diag(1 ./ cdiag) + (U' * AU);
    x = y - AU * (Sinv \ (U' * y));
end

function relerr = relativeError(xRef, xTest)
    relerr = norm(xRef - xTest) / max(1.0, norm(xRef));
end
