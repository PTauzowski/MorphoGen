classdef TestMultiGroupAssembly < matlab.unittest.TestCase
%TESTMULTIGROUPASSEMBLY Assembly bookkeeping across several element groups.
%
%   An FEAnalysis may hold more than one finite-element group - several
%   materials on one mesh, or frames and solids sharing one mesh, which is the
%   case develop's DOF refactor exists to support. Call sites write the group
%   cell array both ways: CorbelModelMultiMat uses a row, {fe1 fe2 fe3}, and
%   PlaneStressCracked uses a row of two.
%
%   Every assertion here passes trivially with a single group, which is why
%   this path went untested and why the smoke examples never caught it.

    properties (Constant)
        NX = 4      % elements along x
        NY = 2      % elements along y
    end

    properties
        mesh
        fe1
        fe2
        totalElems
    end

    methods (TestClassSetup)
        function addRepoToPath(testCase)
            here = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(here, 'IncludingSubfolders', true));
        end
    end

    methods (TestMethodSetup)
        function buildTwoGroupsOnOneMesh(testCase)
            sf = ShapeFunctionL4;
            testCase.mesh = Mesh();
            elems = testCase.mesh.addRectMesh2D( ...
                0, 0, testCase.NX, testCase.NY, testCase.NX, testCase.NY, sf.pattern);

            testCase.totalElems = size(elems, 1);
            split = floor(testCase.totalElems / 2);

            mat = PlaneStressMaterial('mat');
            mat.setElasticIzo(1.0E5, 0.3);

            testCase.fe1 = PlaneStressElem(sf, elems(1:split, :));
            testCase.fe1.props.h = 1;
            testCase.fe1.setMaterial(mat);

            testCase.fe2 = PlaneStressElem(sf, elems(split+1:end, :));
            testCase.fe2.props.h = 1;
            testCase.fe2.setMaterial(mat);
        end
    end

    methods (Test)

        function totalElemsCountsEveryGroup(testCase)
            % getTotalElemsNumber sizes the density vector for the whole
            % design/ topology-optimisation stack. Counting one group leaves
            % the optimiser working on a fraction of the mesh.
            for orientation = ["row", "column"]
                analysis = testCase.buildAnalysis(orientation);
                testCase.verifyEqual(analysis.getTotalElemsNumber(), testCase.totalElems, ...
                    sprintf('getTotalElemsNumber must count every group (%s cell array)', orientation));
            end
        end

        function elemIndicesPartitionEveryGroup(testCase)
            % getElemIndices slices the density vector per group - x(ei{k}).
            % It must return one index range per group, together covering
            % every element exactly once.
            for orientation = ["row", "column"]
                analysis = testCase.buildAnalysis(orientation);
                ei = analysis.getElemIndices();

                testCase.verifyNumElements(ei, 2, ...
                    sprintf('one index range per element group (%s cell array)', orientation));
                testCase.verifyEqual(sort([ei{:}]), 1:testCase.totalElems, ...
                    sprintf('ranges must partition the elements (%s cell array)', orientation));
            end
        end

        function weightedAssemblyCoversEveryGroup(testCase)
            % The values array from the weighted aggregation must line up with
            % the index arrays from globalMatrixIndices - sparse(I,J,K) consumes
            % all three together. A dropped group makes them different lengths.
            for orientation = ["row", "column"]
                analysis = testCase.buildWeightedAnalysis(orientation);
                [I, J, ~, ~] = analysis.globalMatrixIndices();
                x = ones(testCase.totalElems, 1);

                K = analysis.globalMatrixAggregationWeighted('computeStifnessMatrix', x);

                testCase.verifyNumElements(K, numel(I), ...
                    sprintf('weighted values must match the index arrays (%s cell array)', orientation));
                testCase.verifyNumElements(J, numel(I), 'I and J must agree');
            end
        end

        function unweightedAssemblyCoversEveryGroup(testCase)
            for orientation = ["row", "column"]
                analysis = testCase.buildAnalysis(orientation);
                [I, ~, ~, ~] = analysis.globalMatrixIndices();

                K = analysis.globalMatrixAggregation('computeStifnessMatrix');

                testCase.verifyNumElements(K, numel(I), ...
                    sprintf('unweighted values must match the index arrays (%s cell array)', orientation));
            end
        end

        function orientationDoesNotChangeResults(testCase)
            % A row and a column cell array describe the same model, so every
            % assembly quantity must be identical between them.
            byRow    = testCase.buildAnalysis("row");
            byColumn = testCase.buildAnalysis("column");

            testCase.verifyEqual(byRow.getTotalElemsNumber(), byColumn.getTotalElemsNumber());
            testCase.verifyEqual(byRow.getElemIndices(), byColumn.getElemIndices());
            testCase.verifyEqual( ...
                byRow.globalMatrixAggregation('computeStifnessMatrix'), ...
                byColumn.globalMatrixAggregation('computeStifnessMatrix'));
        end

        function dofNamesUnionEveryGroup(testCase)
            % obj.ndofs seeds from felems{1} and unions the rest. With groups
            % of one element family the union is of identical name sets, so
            % this asserts the loop runs rather than that it changes anything.
            analysis = testCase.buildAnalysis("row");
            testCase.verifyEqual(sort(analysis.ndofs), sort(testCase.fe1.eDofs), ...
                'one element family: the union is that family''s DOF names');
        end

    end

    methods (Access = private)

        function analysis = buildAnalysis(testCase, orientation)
            analysis = LinearElasticity(testCase.groupCell(orientation), testCase.mesh);
        end

        function analysis = buildWeightedAnalysis(testCase, orientation)
            analysis = LinearElasticityWeighted(testCase.groupCell(orientation), testCase.mesh, false);
        end

        function groups = groupCell(testCase, orientation)
            if orientation == "row"
                groups = { testCase.fe1  testCase.fe2 };     % 1x2
            else
                groups = { testCase.fe1; testCase.fe2 };     % 2x1
            end
        end

    end
end
