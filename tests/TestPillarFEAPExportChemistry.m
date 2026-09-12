classdef TestPillarFEAPExportChemistry < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            dirs = {'math', 'mesh', 'analysis', 'elements', ...
                    fullfile('examples', 'Pillar')};
            for d = dirs
                testCase.applyFixture( ...
                    matlab.unittest.fixtures.PathFixture(fullfile(root, d{1})));
            end
        end
    end

    methods (Test)
        function testEDISUsesLocalChemistryChannels(testCase)
            p = struct();
            p.sf = ShapeFunctionH27();

            p.top_R = 5;
            p.pillar_inclination_deg = 5;
            p.pillar_layers = [2 4];
            p.pillar_res = [1 1];
            p.ground_layers = [3 5];
            p.ground_res = [1 2];

            p.int_th = 0.2;
            p.pillar_chem = [0.30 0.60];
            p.ground_chem = [0 1];
            p.z_offset = -p.pillar_layers(1);

            p.depression_width = NaN;
            p.depression_depth = NaN;
            p.depression_r_min = NaN;
            p.tile_size = NaN;

            p.res_cyl = 2;
            p.res_ring = 2;
            p.res_tile = 2;
            p.res_bank = 1;

            model = PillarModel(p);

            outFile = [tempname '.i'];
            cleanupObj = onCleanup(@() deleteIfPresent(outFile)); %#ok<NASGU>
            model.FEAP_Export(outFile);

            edis = readEDIS(outFile);
            z_local = edis.z - p.z_offset;
            tol = 1e-6;

            bottomPillar = z_local > (p.int_th/2 + tol) & ...
                           z_local < (p.pillar_layers(1) - p.int_th/2 - tol);
            testCase.assertTrue(any(bottomPillar), ...
                'No exported z-level found inside the bottom pillar layer.');
            testCase.verifyTrue(any(edis.z(bottomPillar) < 0), ...
                'Regression setup is invalid: bottom pillar layer did not shift below exported z=0.');
            testCase.verifyLessThanOrEqual( ...
                max(abs(edis.chem1(bottomPillar) - p.pillar_chem(1))), tol);
            testCase.verifyLessThanOrEqual(max(abs(edis.chem2(bottomPillar))), tol);

            margin = 0.5;
            topGround = z_local < -(p.int_th/2 + margin) & ...
                        z_local > -(p.ground_layers(end) - p.int_th/2 - margin);
            testCase.assertTrue(any(topGround), ...
                'No exported z-level found inside the top ground layer.');
            testCase.verifyLessThanOrEqual(max(abs(edis.chem1(topGround))), tol);
            testCase.verifyLessThanOrEqual( ...
                max(abs(edis.chem2(topGround) - p.ground_chem(end))), tol);

            seamMid = abs(z_local) <= tol;
            testCase.assertTrue(any(seamMid), ...
                'No exported z-level found at the centered pillar-ground interface.');
            testCase.verifyLessThanOrEqual( ...
                max(abs(edis.chem1(seamMid) - 0.5 * p.pillar_chem(1))), tol);
            testCase.verifyLessThanOrEqual( ...
                max(abs(edis.chem2(seamMid) - 0.5 * p.ground_chem(end))), tol);
        end
    end
end

function edis = readEDIS(filename)
    lines = regexp(fileread(filename), '\r\n|\n|\r', 'split');
    startIdx = find(strcmp(strtrim(lines), 'EDIS'), 1, 'first');
    if isempty(startIdx)
        error('TestPillarFEAPExportChemistry:MissingEDIS', ...
            'Could not find EDIS block in %s.', filename);
    end

    z = [];
    chem1 = [];
    chem2 = [];
    for idx = startIdx + 1:numel(lines)
        line = strtrim(lines{idx});
        if isempty(line)
            break;
        end

        vals = sscanf(line, '%f');
        if numel(vals) < 7
            continue;
        end

        z(end+1,1) = vals(2); %#ok<AGROW>
        chem1(end+1,1) = vals(6);
        chem2(end+1,1) = vals(7);
    end

    edis = struct('z', z, 'chem1', chem1, 'chem2', chem2);
end

function deleteIfPresent(filename)
    if exist(filename, 'file')
        delete(filename);
    end
end
