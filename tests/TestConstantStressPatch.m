classdef TestConstantStressPatch < matlab.unittest.TestCase
    % Constant-stress patch test — validates inter-element load application.
    %
    % A rectangular beam carries a uniform pressure on one edge (2D) or face
    % (3D) and is shear-supported on the opposite edge/face: restrained normal
    % to it, free to slide, with one node pinned in the remaining directions to
    % remove rigid-body motion.
    %
    % The exact solution is a uniform uniaxial stress state. At EVERY Gauss
    % point of EVERY element the stress in the loading direction must equal the
    % applied traction, and all other components must vanish.
    %
    % This is the only test in the suite that exercises the edge/face load
    % integral, which is where a wrong Jacobian or a lumped-instead-of-
    % consistent nodal force produces plausible-looking wrong answers.
    %
    % Tolerance is machine precision, not engineering tolerance: any element
    % able to represent constant strain reproduces this exactly, on any mesh at
    % any resolution. A result that is merely close is a failure.
    %
    % Assertions are made on Gauss-point values rather than nodal ones. Nodal
    % results are extrapolated and averaged between elements — exactly the
    % smoothing that would mask a bad load integral.

    properties (Constant)
        Length   = 3;        % beam length parameter l
        Res      = 4;        % deliberately coarse: the test is exact at any resolution
        E        = 210000;
        Nu       = 0.3;
        Traction = -100;     % applied traction along x
        RelTol   = 1e-10;    % relative to |Traction|
    end

    properties (TestParameter)
        % 2D: quadrilateral orders plus the triangular variant.
        planeElement = {'Q4', 'Q9', 'Q16', 'T3'};
        % 3D: hexahedral orders.
        solidElement = {'H8', 'H27'};
    end

    methods (TestClassSetup)

        function addRepoToPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture({ ...
                    fullfile(root, 'analysis'), ...
                    fullfile(root, 'elements'), ...
                    fullfile(root, 'mesh'), ...
                    fullfile(root, 'math'), ...
                    fullfile(root, 'materials'), ...
                    fullfile(root, 'design'), ...
                    fullfile(root, 'examples', 'models')}));
        end

        function suppressFigures(testCase)
            previous = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
            testCase.addTeardown(@() set(0, 'DefaultFigureVisible', previous));
        end

    end

    methods (Test)

        function planeStressConstantAtEveryGaussPoint(testCase, planeElement)
            model = testCase.buildPlaneModel(planeElement);
            evalc('model.solveWeighted()');

            % PlaneStressElem leaves gp.stress as (component, elem, ip) —
            % its permute call is commented out. See note in extractComponent.
            testCase.assertStressState( ...
                model.fe.results.gp.stress, 'componentFirst', ...
                [testCase.Traction, 0, 0], ...
                ["sxx" "syy" "sxy"], ...
                planeElement);
        end

        function solidStressConstantAtEveryGaussPoint(testCase, solidElement)
            model = testCase.buildSolidModel(solidElement);
            evalc('model.solveWeighted()');

            % SolidElasticElem permutes gp.stress to (elem, ip, component).
            testCase.assertStressState( ...
                model.fe.results.gp.stress, 'componentLast', ...
                [testCase.Traction, 0, 0, 0, 0, 0], ...
                ["sxx" "syy" "szz" "sxy" "syz" "sxz"], ...
                solidElement);
        end

    end

    methods (Access = private)

        function model = buildPlaneModel(testCase, elementName)
            pressure = [testCase.Traction 0];
            switch elementName
                case 'T3'
                    % Triangular meshes use their own fixture.
                    ctor = @() ConstPlaneStressModelTriangular( ...
                        ShapeFunctionT3, testCase.Length, testCase.Res, ...
                        testCase.E, testCase.Nu, pressure);
                otherwise
                    sf = feval(sprintf('ShapeFunction%s', elementName));
                    ctor = @() ConstPlaneStressModel( ...
                        sf, testCase.Length, testCase.Res, ...
                        testCase.E, testCase.Nu, pressure);
            end
            % Fixtures print problem info on construction; keep test output clean.
            [~, model] = evalc('ctor()');
        end

        function model = buildSolidModel(testCase, elementName)
            sf = feval(sprintf('ShapeFunction%s', elementName));
            % The solid fixture takes a scalar and applies [-pressure 0 0].
            pressure = abs(testCase.Traction);
            ctor = @() ConstStressSolidModel( ...
                sf, testCase.Length, testCase.Res, ...
                testCase.E, testCase.Nu, pressure);
            [~, model] = evalc('ctor()');
        end

        function assertStressState(testCase, gpStress, layout, expected, names, label)
            % NOTE: the two element families store gp.stress with opposite
            % index order. PlaneStressElem has its permute commented out and
            % indexes (component, elem, ip); SolidElasticElem permutes to
            % (elem, ip, component). This is identical on develop, Vibrations
            % and CAS_Arm, so it is pre-existing rather than a merge artefact —
            % but it means `gp.stress(1,:,:)` means "sxx everywhere" for plane
            % elements and "element 1" for solids. Recorded, not fixed (Rule 1).

            testCase.assertNotEmpty(gpStress, ...
                'No Gauss-point stresses were produced.');

            ncomp = numel(expected);
            switch layout
                case 'componentFirst'
                    testCase.assertEqual(size(gpStress, 1), ncomp, ...
                        'Unexpected component count for (component, elem, ip).');
                case 'componentLast'
                    testCase.assertEqual(size(gpStress, 3), ncomp, ...
                        'Unexpected component count for (elem, ip, component).');
                otherwise
                    testCase.assertFail(sprintf('Unknown layout ''%s''.', layout));
            end

            tol = abs(testCase.Traction) * testCase.RelTol;

            for c = 1:ncomp
                switch layout
                    case 'componentFirst'
                        actual = reshape(gpStress(c, :, :), 1, []);
                    case 'componentLast'
                        actual = reshape(gpStress(:, :, c), 1, []);
                end
                worst = max(abs(actual - expected(c)));
                testCase.verifyEqual(actual, ...
                    repmat(expected(c), 1, numel(actual)), ...
                    'AbsTol', tol, ...
                    sprintf(['%s: %s should be %g at every Gauss point ' ...
                             '(worst deviation %g over %d points).'], ...
                            label, names(c), expected(c), worst, numel(actual)));
            end
        end

    end

end
