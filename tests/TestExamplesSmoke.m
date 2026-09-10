classdef TestExamplesSmoke < matlab.unittest.TestCase
    % Smoke harness — does each representative example still RUN?
    %
    % This answers a different question from the unit tests: it does not check
    % that results are correct, only that no code path throws. A silently wrong
    % stiffness matrix still converges to something plausible and passes here,
    % which is why TestConstantStressPatch and the numerical baselines from
    % `main` exist alongside it.
    %
    % Examples are skipped (not failed) when absent, so the same harness can be
    % copied into a detached checkout of any branch to capture its baseline.
    % Heavy topology-optimisation runs are deliberately excluded — this suite is
    % meant to be run after every conflict resolution.

    properties (TestParameter)
        example = { ...
            'elasticity/planeProblems/ConstStressTest.m', ...
            'elasticity/solidProblems/ConstStressSolidTest.m', ...
            'elasticity/planeProblems/CantileverTest.m', ...
            'elasticity/planeProblems/LameProblemTest.m', ...
            'elasticity/solidProblems/CantileverSolidTest.m', ...
            'elasticity/planeProblems/InclinedSupportTest.m'};
    end

    properties (Access = private)
        RepoRoot
    end

    methods (TestClassSetup)

        function addRepoToPath(testCase)
            testCase.RepoRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(testCase.RepoRoot, ...
                    'IncludingSubfolders', true));
        end

        function suppressFigures(testCase)
            previous = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
            testCase.addTeardown(@() set(0, 'DefaultFigureVisible', previous));
            testCase.addTeardown(@close, 'all');
        end

    end

    methods (Test)

        function exampleRunsWithoutError(testCase, example)
            scriptPath = fullfile(testCase.RepoRoot, 'examples', example);
            testCase.assumeTrue(isfile(scriptPath), ...
                sprintf('%s is not present on this branch.', example));

            % Only an exception is a failure — these examples warn legitimately.
            % evalc keeps their printed output out of the test log.
            try
                evalc(sprintf('run(''%s'')', scriptPath));
            catch err
                where = '';
                if ~isempty(err.stack)
                    where = sprintf(' at %s line %d', ...
                        err.stack(1).name, err.stack(1).line);
                end
                testCase.verifyFail(sprintf('%s threw: %s%s', ...
                    example, err.message, where));
            end
        end

    end

end
