classdef TestJ2PlaneStressReturnMapping < matlab.unittest.TestCase
    % J2 (von Mises) plane stress-projected return mapping.
    %
    % Source of truth:
    %   de Souza Neto, Peric & Owen, "Computational Methods for Plasticity:
    %   Theory and Applications", Wiley 2008 — Section 9.4, Boxes 9.3-9.6.
    %
    % The implementation under test is
    %   elements/PlaneStressElastoPlasticElem.returnMapping / NR4RetMap
    %   materials/PlaneStressMaterial.tangentD          (consistent tangent)
    %
    % Box 9.3 defines the yield function in the squared form used by the code:
    %
    %       Phi = 1/2 * s' * P * s  -  1/3 * sy^2       P = 1/3 [ 2 -1 0
    %                                                            -1  2 0
    %                                                             0  0 6 ]
    %
    % Box 9.4 gives the update; Box 9.5 the Newton loop on Phi(dgamma) = 0.
    %
    % These tests exercise returnMapping directly rather than through an
    % analysis, so they check the constitutive update in isolation.

    properties (Constant)
        E   = 210000;
        Nu  = 0.3;
        Sy  = 250;
        Tol = 1e-10;     % relative to Sy^2 — the return mapping is exact
    end

    properties (Access = private)
        fe
        mat
    end

    methods (TestClassSetup)
        function addRepoToPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture({ ...
                    fullfile(root, 'analysis'), fullfile(root, 'elements'), ...
                    fullfile(root, 'mesh'),     fullfile(root, 'math'), ...
                    fullfile(root, 'materials')}));
        end
    end

    methods (TestMethodSetup)
        function buildElement(testCase)
            % A single element is enough; returnMapping is point-wise.
            [~, fe] = evalc( ...
                'PlaneStressElastoPlasticElem(ShapeFunctionL4, [1 2 3 4])');
            [~, m]  = evalc('PlaneStressMaterial(''j2'')');
            evalc(sprintf('m.setElastoPlasticIzo(%g, %g, %g)', ...
                          testCase.E, testCase.Nu, testCase.Sy));
            fe.setMaterial(m);
            testCase.fe  = fe;
            testCase.mat = m;
        end
    end

    methods (Test)

        function equibiaxialTensionYieldsAtSy(testCase)
            % Equibiaxial is one of the two states a fixed strain increment
            % maps onto itself: de = [e; e; 0] gives s_trial = [t; t; 0], and
            % A(dgamma) preserves that form, so the returned state is
            % [s; s; 0]. Box 9.3 then gives s'*P*s = (2/3)s^2, hence
            % Phi = 0  <=>  s = Sy.
            %
            % NOTE: uniaxial tension is NOT such a state. A strain increment
            % producing a uniaxial trial stress returns with syy ~= 0 unless
            % nu = 0.5, because A11 == A22 requires 6G(1-nu) == E. Testing
            % uniaxial stress properly needs a driver that iterates on eyy to
            % hold syy = 0; that is deliberately not done here.
            e  = 3 * testCase.Sy * (1 - testCase.Nu) / testCase.E;
            de = [e; e; 0];

            [s, ~, ~, dg] = testCase.fe.returnMapping( ...
                zeros(3,1), zeros(3,1), de, 1.0);

            testCase.verifyGreaterThan(dg, 0, 'Expected a plastic step.');
            testCase.verifyEqual(s(1), testCase.Sy, 'RelTol', 1e-8, ...
                'Equibiaxial stress must return to Sy.');
            testCase.verifyEqual(s(2), testCase.Sy, 'RelTol', 1e-8, ...
                'Equibiaxial state must stay equibiaxial.');
            testCase.verifyEqual(s(3), 0, 'AbsTol', testCase.Sy*1e-8);
        end

        function equivalentPlasticStrainIsScalarNotBroadcast(testCase)
            % Box 9.4 carries TWO distinct quantities: the plastic strain
            % tensor eps^p and the scalar accumulated equivalent plastic
            % strain ebar^p = ebar^p_n + dgamma*sqrt(2/3*xi).
            %
            % returnMapping stores its third output in results.gp.pstrain, a
            % 3-component array, but computes it with the SCALAR formula — so
            % MATLAB broadcasts one value into all three slots. The result is
            % neither the plastic strain tensor (whose components differ, with
            % eps^p_33 = -(eps^p_11 + eps^p_22)) nor a clean scalar.
            e = 3 * testCase.Sy * (1 - testCase.Nu) / testCase.E;
            [~, ~, ep, ~] = testCase.fe.returnMapping( ...
                zeros(3,1), zeros(3,1), [e; e; 0], 1.0);

            testCase.verifyFalse(ep(1) == ep(2) && ep(2) == ep(3), ...
                ['pstrain holds one scalar broadcast into three components, ' ...
                 'conflating the plastic strain tensor with the equivalent ' ...
                 'plastic strain (Box 9.4).']);
        end

        function pureShearYieldsAtSyOverSqrt3(testCase)
            % Engineering shear strain gamma: trial stress [0; 0; G*gamma].
            % von Mises in pure shear yields at Sy/sqrt(3).
            G     = testCase.E / (2 * (1 + testCase.Nu));
            gamma = 4 * testCase.Sy / (sqrt(3) * G);
            de    = [0; 0; gamma];

            [s, ~, ~, dg] = testCase.fe.returnMapping( ...
                zeros(3,1), zeros(3,1), de, 1.0);

            testCase.verifyGreaterThan(dg, 0, 'Expected a plastic step.');
            testCase.verifyEqual(s(3), testCase.Sy/sqrt(3), 'RelTol', 1e-8, ...
                'Shear stress must return to Sy/sqrt(3).');
            testCase.verifyEqual(s(1), 0, 'AbsTol', testCase.Sy*1e-8);
            testCase.verifyEqual(s(2), 0, 'AbsTol', testCase.Sy*1e-8);
        end

        function returnedStressLiesOnYieldSurface(testCase)
            % Box 9.5 terminates on |Phi| <= tol. Whatever the tolerance, the
            % returned stress must satisfy the consistency condition Phi = 0.
            % Several load paths, all well past yield.
            paths = { ...
                [ 4e-3;  0;      0    ], ...
                [ 3e-3; -1e-3;   0    ], ...
                [ 2e-3;  2e-3;   0    ], ...   % equibiaxial
                [ 0;     0;      6e-3 ], ...   % pure shear
                [ 3e-3;  1e-3;   2e-3 ]};      % mixed

            for k = 1:numel(paths)
                [s, ~, ~, dg] = testCase.fe.returnMapping( ...
                    zeros(3,1), zeros(3,1), paths{k}, 1.0);
                testCase.assertGreaterThan(dg, 0, ...
                    sprintf('Path %d was expected to yield.', k));

                phi = testCase.yieldFunction(s);
                testCase.verifyEqual(phi, 0, ...
                    'AbsTol', testCase.Sy^2 * testCase.Tol, ...
                    sprintf(['Path %d: returned stress is off the yield ' ...
                             'surface. Phi = %g (Sy^2/3 = %g).'], ...
                            k, phi, testCase.Sy^2/3));
            end
        end

        function elasticStepPreservesPlasticStrain(testCase)
            % Box 9.4 step (ii): on an elastic step set (.)_{n+1} := (.)_trial.
            % Accumulated plastic strain must therefore survive unchanged.
            % Yield first, then apply a small elastic unloading increment.
            e   = 3 * testCase.Sy / testCase.E;
            [~, e1, ep1, ~] = testCase.fe.returnMapping( ...
                zeros(3,1), zeros(3,1), [e; -testCase.Nu*e; 0], 1.0);
            testCase.assertGreaterThan(norm(ep1), 0, ...
                'Setup failed: no plastic strain was produced.');

            % A small unloading step stays inside the yield surface.
            deUnload = -0.01 * [e; -testCase.Nu*e; 0];
            [~, ~, ep2, dg2] = testCase.fe.returnMapping(e1, ep1, deUnload, 1.0);

            testCase.verifyEqual(dg2, 0, ...
                'Unloading step should be elastic.');
            testCase.verifyEqual(ep2, ep1, 'AbsTol', 1e-14, ...
                ['Plastic strain must be carried through an elastic step ' ...
                 '(Box 9.4 (ii)), not reset.']);
        end

        function consistentTangentMatchesNumericalDerivative(testCase)
            % Box 9.6: the consistent (algorithmic) tangent is the exact
            % derivative of the updated stress with respect to the strain
            % increment. A wrong tangent still converges to the right answer
            % but ruins Newton convergence, so nothing else catches it.
            de0 = [3e-3; 1e-3; 1e-3];
            [s0, ~, ~, dg0] = testCase.fe.returnMapping( ...
                zeros(3,1), zeros(3,1), de0, 1.0);
            testCase.assertGreaterThan(dg0, 0, 'Expected a plastic step.');

            Dt = testCase.mat.tangentD(s0, dg0);

            h  = 1e-7;
            Dnum = zeros(3);
            for j = 1:3
                dep = de0;  dep(j) = dep(j) + h;
                dem = de0;  dem(j) = dem(j) - h;
                sp = testCase.fe.returnMapping(zeros(3,1), zeros(3,1), dep, 1.0);
                sm = testCase.fe.returnMapping(zeros(3,1), zeros(3,1), dem, 1.0);
                Dnum(:, j) = (sp - sm) / (2*h);
            end

            scale = norm(Dnum, 'fro');
            testCase.verifyEqual(Dt, Dnum, 'AbsTol', 1e-4 * scale, ...
                'Consistent tangent does not match the numerical derivative.');
        end

    end

    methods (Access = private)

        function phi = yieldFunction(testCase, s)
            % Box 9.3, squared form.
            P   = (1/3) * [2 -1 0; -1 2 0; 0 0 6];
            phi = 0.5 * (s' * P * s) - (1/3) * testCase.Sy^2;
        end

    end

end
