function results = runTests(varargin)
%RUNTESTS Run the MorphoGen test suite.
%
%   runTests                 % everything under tests/
%   runTests('patch')        % constant-stress patch test only
%   runTests('smoke')        % example smoke harness only
%
%   results = runTests(...)  % returns the matlab.unittest.TestResult array
%
% Exits with a non-zero status when run non-interactively (matlab -batch),
% so it can gate a merge step or, later, CI.

    here = fileparts(mfilename('fullpath'));

    switch lower(strjoin(varargin, ''))
        case {'', 'all'}
            suite = matlab.unittest.TestSuite.fromFolder(here);
        case 'patch'
            suite = matlab.unittest.TestSuite.fromFile( ...
                fullfile(here, 'TestConstantStressPatch.m'));
        case 'smoke'
            suite = matlab.unittest.TestSuite.fromFile( ...
                fullfile(here, 'TestExamplesSmoke.m'));
        otherwise
            error('runTests:unknownSelection', ...
                'Unknown selection ''%s''. Use ''all'', ''patch'' or ''smoke''.', ...
                strjoin(varargin, ' '));
    end

    runner = matlab.unittest.TestRunner.withTextOutput();
    results = runner.run(suite);

    summary = table(results);
    fprintf('\n%s\n', repmat('-', 1, 60));
    fprintf('passed %d   failed %d   incomplete %d   skipped %d   (%.1f s)\n', ...
        nnz([results.Passed]), nnz([results.Failed]), ...
        nnz([results.Incomplete]), nnz(~[results.Passed] & ~[results.Failed] ...
                                       & ~[results.Incomplete]), ...
        sum([results.Duration]));
    fprintf('%s\n', repmat('-', 1, 60));

    if nargout > 0
        return;   % caller wants the results; never exit out from under them
    end

    disp(summary);

    % Non-zero exit for batch use: `matlab -batch "runTests"` gates a step.
    if ~usejava('desktop') && any([results.Failed])
        exit(1);
    end
    clear results;
end
