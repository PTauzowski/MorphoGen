%% ========================================================================
%  MMC Topology Optimization for Natural Frequency Maximization
%  ========================================================================
%
%  Description:
%    Maximizes the first natural frequency of a 2D structure using the
%    Moving Morphable Components (MMC) method. The structure has fixed
%    supports at both ends and a concentrated mass at the center.
%
%  Reference:
%    Huang, X. et al. - Topology optimization using MMC method
%    (Appendix F/G of the original paper)
%
%  Method Overview:
%    1. Components are represented as super-ellipsoids with 6 design vars:
%       [x0, y0, L, t1, t2, theta] = [center_x, center_y, half-length,
%                                      half-width1, half-width2, angle]
%    2. Component shapes are aggregated using KS (Kreisselmeier-Steinhauser)
%    3. Density field derived via smoothed Heaviside function
%    4. Design updated using MMA (Method of Moving Asymptotes)
%
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. PROBLEM PARAMETERS
% =========================================================================

% --- Geometry ---
domainWidth  = 4;       % [m] Width of design domain
domainHeight = 1;       % [m] Height of design domain
thickness    = 0.01;    % [m] Out-of-plane thickness

% --- Mesh ---
numElemX = 400;         % Number of elements in x-direction
numElemY = 100;         % Number of elements in y-direction

% --- Material Properties ---
youngsModulus = 2e11;   % [Pa] Young's modulus (steel)
poissonsRatio = 0.3;    % [-] Poisson's ratio
density       = 7800;   % [kg/m^3] Material density
concentratedMass = 1e5; % [kg] Concentrated mass at center

% --- MMC Component Parameters ---
componentSpacingX = 0.25;   % [m] Initial x-spacing between component centers
componentSpacingY = 0.25;   % [m] Initial y-spacing between component centers
initialHalfLength = 0.4;    % [m] Initial half-length of components
initialHalfWidth1 = 0.04;   % [m] Initial half-width (end 1)
initialHalfWidth2 = 0.04;   % [m] Initial half-width (end 2)
initialAngle      = pi/4;   % [rad] Initial inclination angle

% --- Optimization Parameters ---
volumeFraction     = 0.3;   % [-] Target volume fraction (upper bound)
superEllipsoidPower = 6;    % [-] Power for super-ellipsoid shape
ksParameter        = 100;   % [-] KS aggregation parameter

% --- Convergence Settings ---
maxIterations    = 500;     % Maximum optimization iterations
convergenceTol   = 1e-4;    % Convergence tolerance (relative change)
moveLimit        = 0.1;     % Maximum move per iteration (relative)

% --- Epsilon Schedule (Heaviside sharpness) ---
eps_initial = 0.02;         % Initial epsilon value
eps_final   = 0.3;          % Final epsilon value
eps_midIter = 30;           % Iteration at sigmoid midpoint
eps_rate    = 0.5;          % Sigmoid rate parameter

% --- Eigenvalue Analysis ---
targetMode      = 1;        % Which eigenmode to optimize (1 = fundamental)
numModesCompute = 6;        % Number of modes to compute

% --- Numerical Parameters ---
minDensity        = 1e-6;   % Minimum density (void)
sensitivityDigits = 5;      % Significant digits for sensitivity truncation
objScaleFactor    = 100;    % Initial objective scaling factor

%% ========================================================================
%  2. FINITE ELEMENT MESH SETUP
% =========================================================================

% Helper function for sparse matrix assembly
fsparse = @(i,j,s,sz) sparse(i,j,s,sz(1),sz(2));

% Mesh dimensions
numElements = numElemX * numElemY;
numNodes    = (numElemX + 1) * (numElemY + 1);
numDofs     = 2 * numNodes;

% Element dimensions
elemWidth  = domainWidth / numElemX;
elemHeight = domainHeight / numElemY;
minElemSize = min(elemWidth, elemHeight);

% Element stiffness and mass matrices
Ke_tril = computeElementStiffness(youngsModulus, poissonsRatio, ...
                                   elemWidth, elemHeight, thickness);
KE = triLowerToFull(Ke_tril, 8);

Me_tril = computeElementMass(density, elemWidth, elemHeight, thickness);
ME = triLowerToFull(Me_tril, 8);

% Node numbering matrix
nodeMatrix = int32(reshape(1:numNodes, 1+numElemY, 1+numElemX));

% DOF connectivity matrix (8 DOFs per element: 2 per node x 4 nodes)
edofVec = reshape(2*nodeMatrix(1:end-1, 1:end-1) - 1, numElements, 1);
edofMat = edofVec + int32([0 1 2*numElemY+[2 3 4 5] 2 3]);

% Node IDs for each element (for density interpolation)
elementNodeIDs = edofMat(:, 2:2:8) ./ 2;

% Assembly index arrays (lower triangular, for symmetric matrices)
[indexI, indexII] = deal([]);
for j = 1:8
    indexI  = cat(2, indexI, j:8);
    indexII = cat(2, indexII, repmat(j, 1, 8-j+1));
end
[iK, jK] = deal(edofMat(:, indexI)', edofMat(:, indexII)');
assemblyIndex = sort([iK(:), jK(:)], 2, 'descend');
clear iK jK;

% Nodal coordinates (centered at origin)
[nodeX, nodeY] = meshgrid(elemWidth  * (-numElemX/2 : numElemX/2), ...
                          elemHeight * (-numElemY/2 : numElemY/2));
gridCoords.x = nodeX(:);
gridCoords.y = nodeY(:);

% Volume weight for each node (for constraint computation)
volumeWeight = sparse(double(elementNodeIDs(:)), 1, 1/4);

%% ========================================================================
%  3. BOUNDARY CONDITIONS & LOADS
% =========================================================================

% Fixed nodes: left and right edges (clamped-clamped beam)
fixedNodes = union(1:numElemY+1, ...
                   numElemX*(numElemY+1)+1 : (numElemX+1)*(numElemY+1));
fixedDofs  = [2*fixedNodes-1, 2*fixedNodes];
freeDofs   = setdiff(1:numDofs, fixedDofs);

% Elements at fixed boundaries (for connectivity check)
fixedElements = union(1:numElemY, (numElemX-1)*numElemY+1 : numElemY*numElemX);

% Concentrated mass location: center of domain
massNodeID  = (numElemY+1) * numElemX/2 + numElemY/2 + 1;
massDofs    = 2*massNodeID-1 : 2*massNodeID;

% Elements near mass location (for connectivity check)
massElements = [(numElemX-1)*numElemY/2, ...
                (numElemX-1)*numElemY/2+1, ...
                (numElemX+1)*numElemY/2, ...
                (numElemX+1)*numElemY/2+1];

%% ========================================================================
%  4. MMC COMPONENT INITIALIZATION
% =========================================================================

% Generate initial component center coordinates
centerX = componentSpacingX : 2*componentSpacingX : domainWidth;
centerY = componentSpacingY : 2*componentSpacingY : domainHeight;
numCompX = length(centerX);
numCompY = length(centerY);

% Duplicate coordinates for crossing components at each location
centerX = kron(centerX, ones(1, 2*numCompY));
centerY = repmat(kron(centerY, ones(1, 2)), 1, numCompX);

numComponents = length(centerX);

% Initialize design variables for all components
% Each component: [x0, y0, L, t1, t2, theta]
halfLengths = repmat(initialHalfLength, 1, numComponents);
halfWidths1 = repmat(initialHalfWidth1, 1, numComponents);
halfWidths2 = repmat(initialHalfWidth2, 1, numComponents);
angles      = repmat([initialAngle, -initialAngle], 1, numComponents/2);

% Assemble design variable vector
% Shift centers to be relative to domain center
designVars = [centerX - elemWidth*numElemX/2; ...
              centerY - elemHeight*numElemY/2; ...
              halfLengths; halfWidths1; halfWidths2; angles];

numDesignVars    = numel(designVars);
varsPerComponent = numDesignVars / numComponents;

% Active components and design variables (may be reduced during optimization)
activeComponents = 1:numComponents;
activeDesignVars = 1:numDesignVars;

% Non-design region (empty for this problem)
nonDesignTDF = [];

% Initialize component TDF matrix
componentTDF = zeros(numNodes, numComponents);

%% ========================================================================
%  5. MMA OPTIMIZER SETUP
% =========================================================================

% MMA parameters
numConstraints = 1;
mma_c = 1000 * ones(numConstraints, 1);
mma_d = zeros(numConstraints, 1);
mma_a0 = 1;
mma_a = zeros(numConstraints, 1);

% Design variable history
xval  = designVars(:);
xold1 = xval;
xold2 = xval;

% Design variable bounds (per component)
xmin_single = [-domainWidth/2; -domainHeight/2; minElemSize; ...
               minElemSize; minElemSize; -pi];
xmax_single = [domainWidth/2; domainHeight/2; ...
               sqrt(domainWidth^2 + domainHeight^2)/2; ...
               0.5*min(domainWidth, domainHeight) * [1; 1]; pi];

xmin = repmat(xmin_single, numComponents, 1);
xmax = repmat(xmax_single, numComponents, 1);

% MMA asymptote bounds
low = xmin;
upp = xmax;

% Best design tracking (for recovery from disconnection)
xval_best  = xval;
bestFrequency = 0;

%% ========================================================================
%  6. OPTIMIZATION LOOP
% =========================================================================

% History storage
objectiveHistory   = zeros(1, maxIterations);
constraintHistory  = zeros(1, maxIterations);
objRelativeChange  = 1.0;
iter = 1;

fprintf('\n=== Starting Optimization ===\n');
fprintf('Target volume fraction: %.2f\n', volumeFraction);
fprintf('Number of components: %d\n', numComponents);
fprintf('Number of design variables: %d\n\n', numDesignVars);

while objRelativeChange > convergenceTol && iter <= maxIterations

    % --- Compute adaptive epsilon (Heaviside sharpness) ---
    epsilon = eps_initial + eps_final - ...
              (1/eps_final + exp(-eps_rate*(iter - eps_midIter)))^(-1);

    % =====================================================================
    %  STEP 1: Compute Component TDFs and Derivatives
    % =====================================================================
    componentTDFderiv = sparse(numNodes, numDesignVars);

    for iComp = activeComponents
        [componentTDF, componentTDFderiv, xval, activeComponents, activeDesignVars] = ...
            computeComponentTDF(componentTDF, componentTDFderiv, xval, iComp, ...
                               gridCoords, superEllipsoidPower, varsPerComponent, ...
                               epsilon, activeComponents, activeDesignVars, minElemSize);
    end

    % Aggregate TDFs using KS function
    activeTDF = [componentTDF(:, activeComponents), nonDesignTDF];
    expTDF    = exp(ksParameter * activeTDF);
    globalTDF = max(-1e3, log(sum(expTDF, 2)) / ksParameter);

    % Global TDF sensitivity w.r.t. component TDFs
    activeTDFderiv = componentTDFderiv(:, activeDesignVars);
    dPhimax_dPhi   = expTDF(:, 1:length(activeComponents)) ./ (sum(expTDF, 2) + eps);
    dPhimax_dPhi   = kron(dPhimax_dPhi, ones(1, varsPerComponent));
    globalTDFderiv = dPhimax_dPhi .* activeTDFderiv;

    % =====================================================================
    %  STEP 2: Visualization
    % =====================================================================
    figure(1); clf;
    contourf(reshape(nodeX, [numElemY+1, numElemX+1]), ...
             reshape(nodeY, [numElemY+1, numElemX+1]), ...
             reshape(globalTDF, [numElemY+1, numElemX+1]), [0, 0], 'LineColor', 'none');
    colormap([1 1 1; 0 0 0]);  % White = void, Black = material
    axis equal;
    axis([-domainWidth/2, domainWidth/2, -domainHeight/2, domainHeight/2]);
    title(sprintf('Iteration %d | Frequency: %.1f', iter, bestFrequency));
    drawnow;

    % =====================================================================
    %  STEP 3: Finite Element Analysis
    % =====================================================================

    % Compute element densities from global TDF
    heavisideH    = smoothHeaviside(globalTDF, minDensity, epsilon);
    elementDensity = sum(heavisideH(elementNodeIDs), 2) / 4;

    % Initialize eigenvector
    modeShape = zeros(numDofs, 1);

    % Check structural connectivity (mass must be connected to supports)
    [hasLoadPath, connectedLabel] = checkStructuralConnectivity(...
        elementDensity, numElemY, numElemX, minDensity, fixedElements, massElements);

    if hasLoadPath
        isConnected = true;

        % Find elements in connected region
        bwLabels = reshape(bwlabel(reshape(elementDensity, numElemY, numElemX) > minDensity, 4), ...
                          numElemX*numElemY, 1);
        connectedElements = find(bwLabels == connectedLabel);

        % Use only connected elements for FEA
        connectedDensity = zeros(size(elementDensity));
        connectedDensity(connectedElements) = elementDensity(connectedElements);

        edofMatConnected = edofMat(connectedElements, :);
        freeDofConnected = setdiff(edofMatConnected(:), fixedDofs);

        % Build assembly index for connected region
        [iK1, jK1] = deal(edofMatConnected(:, indexI)', edofMatConnected(:, indexII)');
        assemblyIndexConnected = sort([iK1(:), jK1(:)], 2, 'descend');
        clear iK1 jK1;

        % Assemble stiffness and mass matrices
        [K, M] = assembleGlobalMatrices(Ke_tril, Me_tril, connectedDensity, ...
                    connectedElements, assemblyIndexConnected, numDofs, ...
                    massDofs, concentratedMass, fsparse);

        % Solve eigenvalue problem
        [eigenvectors, eigenvalues] = eigs(K(freeDofConnected, freeDofConnected), ...
                                           M(freeDofConnected, freeDofConnected), ...
                                           numModesCompute, 'sm');
        frequency = eigenvalues(targetMode, targetMode);
        modeShape(freeDofConnected) = eigenvectors(:, targetMode);
        modeShape = modeShape / sqrt(modeShape' * M * modeShape);

        % Track best design
        if frequency > bestFrequency
            bestFrequency = frequency;
            xval_best = xval;
        end
    else
        % No load path - penalize disconnected design
        isConnected = false;
        fprintf('  WARNING: Structure disconnected - applying penalty\n');

        % Use full domain for FEA
        [K, M] = assembleGlobalMatrices(Ke_tril, Me_tril, elementDensity, ...
                    1:numElements, assemblyIndex, numDofs, ...
                    massDofs, concentratedMass, fsparse);

        [eigenvectors, eigenvalues] = eigs(K(freeDofs, freeDofs), ...
                                           M(freeDofs, freeDofs), ...
                                           numModesCompute, 'sm');
        frequency = eigenvalues(targetMode, targetMode);
        modeShape(freeDofs) = eigenvectors(:, targetMode);
        modeShape = modeShape / sqrt(modeShape' * M * modeShape);

        % Heavy penalty for disconnection
        if bestFrequency > 0
            frequency = 0.1 * bestFrequency;
        end
        freeDofConnected = freeDofs;  % For sensitivity calculation
    end

    % Store objective and constraint
    objectiveHistory(iter) = frequency;
    volumeConstraint = sum(elementDensity) * elemWidth * elemHeight / ...
                       (domainWidth * domainHeight) - volumeFraction;
    constraintHistory(iter) = volumeConstraint + volumeFraction;

    % =====================================================================
    %  STEP 4: Sensitivity Analysis
    % =====================================================================

    % Derivative of Heaviside function
    dH_dPhi = 3*(1 - minDensity) / (4*epsilon) * (1 - globalTDF.^2 / epsilon^2);
    dH_dPhi(abs(globalTDF) > epsilon) = 0;

    % Adjoint method for eigenvalue sensitivity
    [adjointVec, adjointScalar] = computeAdjointSensitivity(K, M, freeDofConnected, ...
                                                             modeShape, numDofs, frequency);

    % Element sensitivity
    elemSensitivity = sum((adjointVec(edofMat) * (KE - frequency*ME)) .* modeShape(edofMat) - ...
                          0.5 * (adjointScalar * modeShape(edofMat) * ME) .* modeShape(edofMat), 2);

    % Convert to nodal sensitivity
    elemSensNode = elemSensitivity * ones(1, 4) / 4;
    nodalSensitivity = sparse(double(elementNodeIDs(:)), 1, elemSensNode(:));

    % Chain rule: objective sensitivity w.r.t. design variables
    df0dx = zeros(1, numDesignVars);
    df0dx(activeDesignVars) = (nodalSensitivity .* dH_dPhi)' * globalTDFderiv;

    % Constraint sensitivity w.r.t. design variables
    dfdx = zeros(1, numDesignVars);
    dfdx(activeDesignVars) = (volumeWeight .* dH_dPhi)' * globalTDFderiv * ...
                              elemWidth * elemHeight / (domainWidth * domainHeight);

    % Sensitivity truncation (numerical stability)
    digits = sensitivityDigits - floor(log10([max(abs(df0dx(:)))+eps, ...
                                               max(abs(dfdx(:)))+eps]));

    % Objective scaling
    scale = max(abs(df0dx));
    if scale < eps || ~isfinite(scale)
        scale = objScaleFactor;
    else
        objScaleFactor = scale;
    end

    f0val = -frequency / scale;
    df0dx = round(df0dx * 10^digits(1)) / 10^digits(1) / scale;
    dfdx  = round(dfdx * 10^digits(2)) / 10^digits(2);

    % =====================================================================
    %  STEP 5: Design Update (MMA)
    % =====================================================================

    % Apply move limits
    xmin_move = max(xmin, xval - moveLimit*(xmax - xmin));
    xmax_move = min(xmax, xval + moveLimit*(xmax - xmin));

    % Call MMA optimizer
    [xmma, ~,~,~,~,~,~,~,~, low, upp] = mmasub(numConstraints, numDesignVars, iter, ...
        xval, xmin_move, xmax_move, xold1, xold2, f0val, df0dx(:), ...
        volumeConstraint, dfdx, low, upp, mma_a0, mma_a, mma_c, mma_d);

    % Update design variable history
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;

    % =====================================================================
    %  STEP 6: Convergence Check
    % =====================================================================

    if iter >= 5 && volumeConstraint/volumeFraction < 1e-4
        recentObj = objectiveHistory(iter-4:iter);
        objRelativeChange = abs(max(abs(recentObj - mean(recentObj))) / mean(recentObj));
    end

    % Print iteration info
    fprintf('It.: %4d  |  Freq: %8.2f  |  Vol: %+7.4f  |  Change: %.4f\n', ...
            iter, frequency, volumeConstraint, objRelativeChange);

    iter = iter + 1;
end

fprintf('\n=== Optimization Complete ===\n');
fprintf('Final frequency: %.2f\n', bestFrequency);
fprintf('Total iterations: %d\n', iter-1);

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function [TDF, TDFderiv, xval, activeComp, activeVars] = ...
    computeComponentTDF(TDF, TDFderiv, xval, compIdx, grid, p, varsPerComp, ...
                        epsilon, activeComp, activeVars, minSize)
% COMPUTECOMPONENTTDF Compute TDF and derivatives for a single MMC component
%
%   The component is a super-ellipsoid defined by:
%     phi(x,y) = 1 - ((x'/L)^p + (y'/l)^p)^(1/p)
%   where (x',y') are local coordinates and l varies linearly with x'.

    % Extract component parameters
    varIdx = (compIdx-1)*varsPerComp + 1 : compIdx*varsPerComp;
    params = xval(varIdx);

    x0 = params(1);     % Center x
    y0 = params(2);     % Center y
    L  = params(3) + eps;  % Half-length (avoid division by zero)
    t1 = params(4);     % Half-width at one end
    t2 = params(5);     % Half-width at other end
    theta = params(6);  % Rotation angle

    % Rotation to local coordinates
    st = sin(theta);
    ct = cos(theta);
    x_local = ct*(grid.x - x0) + st*(grid.y - y0) + eps;
    y_local = -st*(grid.x - x0) + ct*(grid.y - y0) + eps;

    % Variable half-width along component length
    halfWidth = (t1 + t2)/2 + (t2 - t1)/2/L * x_local + eps;

    % Super-ellipsoid implicit function
    superEllipsoid = abs(x_local).^p / L^p + abs(y_local).^p ./ halfWidth.^p;
    TDF(:, compIdx) = 1 - superEllipsoid.^(1/p);

    % Check if component should be deleted (too small or outside domain)
    if (t1/minSize < 1.01 && t2/minSize < 1.01) || min(abs(TDF(:, compIdx))) >= epsilon
        TDF(:, compIdx) = -1e3;
        xval(varIdx(4:5)) = 0;
        activeComp = setdiff(activeComp, compIdx);
        activeVars = setdiff(activeVars, varIdx);
        return;
    end

    % --- Compute derivatives ---
    % d(x_local)/d(params)
    dx_local = [-ct + 0*x_local, -st + 0*x_local, 0*x_local, 0*x_local, 0*x_local, y_local];
    dy_local = [st + 0*y_local, -ct + 0*y_local, 0*y_local, 0*y_local, 0*y_local, -x_local];

    % d(halfWidth)/d(params)
    dw_dx = (t2 - t1) / (2*L);
    dw_dL = -(t2 - t1)/2 * x_local / L^2;
    dw_dt1 = 1/2 - x_local/(2*L);
    dw_dt2 = 1/2 + x_local/(2*L);
    dHalfWidth = [0*halfWidth, 0*halfWidth, dw_dL, dw_dt1, dw_dt2, 0*halfWidth] + ...
                 repmat(dw_dx, 1, varsPerComp) .* dx_local;

    % d(phi)/d(x_local), d(phi)/d(y_local), d(phi)/d(halfWidth)
    dPhi_dx = -(superEllipsoid).^(1/p - 1) .* (x_local/L).^(p-1) / L;
    dPhi_dy = -(superEllipsoid).^(1/p - 1) .* (y_local./halfWidth).^(p-1) ./ halfWidth;
    dPhi_dw = (superEllipsoid).^(1/p - 1) .* (y_local./halfWidth).^p ./ halfWidth;

    % d(phi)/d(L) - special case
    dPhi_dL = zeros(size(dx_local));
    dPhi_dL(:, 3) = (superEllipsoid).^(1/p - 1) .* (x_local/L).^p / L;

    % Chain rule for all parameters
    TDFderiv(:, varIdx) = repmat(dPhi_dx, 1, varsPerComp) .* dx_local + ...
                          repmat(dPhi_dy, 1, varsPerComp) .* dy_local + ...
                          dPhi_dL + ...
                          repmat(dPhi_dw, 1, varsPerComp) .* dHalfWidth;
end


function [adjoint1, adjoint2] = computeAdjointSensitivity(K, M, freeDofs, psi, nDof, lambda)
% COMPUTEADJOINTSENSITIVITY Solve adjoint system for eigenvalue sensitivity
%
%   Solves the augmented system:
%     [K - lambda*M,  -M*psi] [a1]   [0]
%     [-psi'*M,          0  ] [a2] = [1]

    L11 = K(freeDofs, freeDofs) - lambda * M(freeDofs, freeDofs);
    L12 = -M(freeDofs, freeDofs) * psi(freeDofs);
    L21 = L12';

    augmentedMatrix = [L11, L12; L21, 0];
    rhs = [zeros(length(freeDofs), 1); 1];

    solution = augmentedMatrix \ rhs;

    adjoint1 = zeros(nDof, 1);
    adjoint1(freeDofs) = solution(1:end-1);
    adjoint2 = solution(end);
end


function H = smoothHeaviside(phi, alpha, epsilon)
% SMOOTHHEAVISIDE Smoothed Heaviside function for density interpolation
%
%   H = alpha + (1-alpha) * smooth_step(phi/epsilon)
%   where smooth_step is a cubic polynomial in the transition region.

    H = 3*(1 - alpha)/4 * (phi/epsilon - phi.^3/(3*epsilon^3)) + (1 + alpha)/2;
    H(phi > epsilon)  = 1;
    H(phi < -epsilon) = alpha;
end


function [hasPath, connectedLabel] = checkStructuralConnectivity(density, nely, nelx, ...
                                                                  minDen, fixEle, massEle)
% CHECKSTRUCTURALCONNECTIVITY Check if mass is connected to supports
%
%   Uses connected component labeling to verify load path exists.

    binaryDensity = reshape(density, nely, nelx) > minDen;
    labeledRegions = reshape(bwlabel(binaryDensity, 4), nelx*nely, 1);

    fixedLabels = unique(labeledRegions(fixEle));
    massLabels  = nonzeros(unique(labeledRegions(massEle)));

    hasPath = false;
    connectedLabel = [];

    if ~isempty(massLabels)
        for label = massLabels'
            if ismember(label, fixedLabels)
                hasPath = true;
                connectedLabel = label;
                return;
            end
        end
    end
end


function [K, M] = assembleGlobalMatrices(Ke, Me, density, elements, assemblyIdx, ...
                                          nDof, massDofs, massValue, fsparse)
% ASSEMBLEGLOBALMATRICES Assemble global stiffness and mass matrices
%
%   Uses vectorized assembly for efficiency.

    % Stiffness matrix
    sK = reshape(Ke(:) * density(elements)', length(Ke)*length(elements), 1);
    K = fsparse(assemblyIdx(:,1), assemblyIdx(:,2), sK, [nDof, nDof]);
    K = K + K' - diag(diag(K));
    K = K + fsparse(1:nDof, 1:nDof, eps*ones(1,nDof), [nDof, nDof]);  % Regularization

    % Mass matrix
    sM = reshape(Me(:) * density(elements)', length(Me)*length(elements), 1);
    M = fsparse(assemblyIdx(:,1), assemblyIdx(:,2), sM, [nDof, nDof]);
    M = M + M' - diag(diag(M));
    M(massDofs, massDofs) = M(massDofs, massDofs) + massValue * eye(length(massDofs));
    M = M + fsparse(1:nDof, 1:nDof, eps*ones(1,nDof), [nDof, nDof]);  % Regularization
end


function fullMatrix = triLowerToFull(triLower, n)
% TRILOWERTOFULL Convert lower triangular vector to full symmetric matrix

    fullMatrix = zeros(n, n);
    fullMatrix(tril(ones(n)) == 1) = triLower';
    fullMatrix = fullMatrix + fullMatrix' - diag(diag(fullMatrix));
end


function Ke = computeElementStiffness(E, nu, a, b, h)
% COMPUTEELEMENTSTIFFNESS 4-node quadrilateral element stiffness matrix
%
%   Returns lower triangular part (36 unique entries for 8x8 symmetric matrix)

    k1 = [-1/6/a/b*(nu*a^2-2*b^2-a^2), 1/8*nu+1/8, ...
          -1/12/a/b*(nu*a^2+4*b^2-a^2), 3/8*nu-1/8, ...
           1/12/a/b*(nu*a^2-2*b^2-a^2), -1/8*nu-1/8, ...
           1/6/a/b*(nu*a^2+b^2-a^2), -3/8*nu+1/8];

    k2 = [-1/6/a/b*(nu*b^2-2*a^2-b^2), 1/8*nu+1/8, ...
          -1/12/a/b*(nu*b^2+4*a^2-b^2), 3/8*nu-1/8, ...
           1/12/a/b*(nu*b^2-2*a^2-b^2), -1/8*nu-1/8, ...
           1/6/a/b*(nu*b^2+a^2-b^2), -3/8*nu+1/8];

    Ke = E*h/(1-nu^2) * ...
        [k1(1); k1(2); k1(3); k1(4); k1(5); k1(6); ...
         k1(7); k1(8); k2(1); k2(8); k2(7); k2(6); ...
         k2(5); k2(4); k2(3); k1(1); k1(6); k1(7); ...
         k1(4); k1(5); k1(2); k2(1); k2(8); k2(3); ...
         k2(2); k2(5); k1(1); k1(2); k1(3); k1(4); ...
         k2(1); k2(8); k2(7); k1(1); k1(6); k2(1)];
end


function Me = computeElementMass(rho, EW, EH, h)
% COMPUTEELEMENTMASS 4-node quadrilateral consistent mass matrix
%
%   Returns lower triangular part (36 unique entries for 8x8 symmetric matrix)

    A = EW * EH;
    Me = (rho*A*h/36) * ...
        [4; 0; 2; 0; 1; 0; ...
         2; 0; 4; 0; 2; 0; ...
         1; 0; 2; 4; 0; 2; ...
         0; 1; 0; 4; 0; 2; ...
         0; 1; 4; 0; 2; 0; ...
         4; 0; 2; 4; 0; 4];
end
