# MorphoGen

Computational MorphoGenesis System

## Authors
[Piotr Tauzowski](mailto:ptauzow@ippt.pan.pl)

## Introduction
The __MorphoGen__ system is used to solve tasks in the field of structural topology optimization. The system is written in MATLAB using an object-oriented programming paradigm using a layered architecture allows easy development of the system. The system also includes a finite element and reliability layer. Structural analysis is an essential element of structural topological optimization; moreover, the system can only be used for FEM structural analysis and allows for easy expansion in scope. new FEM analyzes and new types of finite elements. MorphoGen directory consists of several subdirectories with source code described below:

* __examples__ - here users find ready-to-use example scripts showcasing the implementation and use of the software.
* design layer implements algorithms for topology optimization such
as derivative-free approach called stress intensity driven topology optimization, as well as sensitivity-based approach known as Method of
Moving Asymptotes.
* __analysis__ layer contains algorithms dedicated to solving specific problems for e.g.: algorithms for topology optimization such
as derivative-free approach called stress intensity driven topology optimization, as well as sensitivity-based approach known as Method of
Moving Asymptotes. This directory contains also finite element algorithms for example: linear elastic analysis with sensitivity, elastoplastic analy-
sis, etc. Additionally, reliability algorithms such as FORM (First-Order
Reliability Method), AMV (Advanced Mean Value), and Monte Carlo
simulation are stored here.
* __element__ stores the finite element files, including planar and
spatial element definitions.
* __material__ layer contains material definitions specific to each finite element type, e.g. isotropic materials are defined for plane stress/strain
and solid elements.
* __mesh__ contains nodes and elements organized into mesh class
designed for finite element mesh storage and generation.
* __math__ provides various mathematical utilities, including Gauss
numerical integration and linear equation solver class. They facilitate
data manipulation of global matrices and right-hand vectors. Addi-
tionally, shape function classes compatible with finite elements are included.
These directories provide a comprehensive framework to explore and utilize the software’s features and capabilities.
5


## Features
The basic functions of the system include the following:  
* Topology optimization in based on stress intensity gradientless approach.
* Topology optimization gradient based algorithm with application of the Method of Moving Asymptotes (MMA) implementented by Krister Svanberg in MATLAB.
* Various constraint: volume, compliance, reliability.
* Finite Element Analysis for plane stress and solid elements.
* Reliability analysis: Monte Carlo, FORM, HMV methods.

## Installation
After cloning the repository, the 'MorphoGen' directory will appear. After adding it to the path along with its subdirectories, sysyem is ready for use. The user can enter the 'examples' subdirectory and run sample task files.

## Configuration
The system is currently written and tested in MATLAB R2023a. So it is recommended to use this version. A regular adding MorphoGen directory to the path along with its subdirectories is needed for be able to use the system.

## Usage
We will present here a file with the definition of an example topological optimization of a cantilever using two methods:

```matlab
clear;
close all;

% Cantilever topology optimization elastic task.

% Resolution of shortest (vertical) edge
res = 50;

% height of the cantilever
h = 1;

% Aspect ratio length/height
aspect=2;

% Filtering radius
Rfilter = 2*h/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

% Type of shape function to be used (here: four node Langrange)
sfL4 = ShapeFunctionL4;

% Creating FE mesh object
mesh = Mesh();

% Generating rectangular mesh ( aspect*h x h )
mesh.addRectMesh2D(0, 0, aspect*h, h, aspect*res, res, sfL4.pattern);

% Create plane stress finite element object
fe=PlaneStressElem( sfL4, mesh.elems );

% Create isotropic material object
material = PlaneStressMaterial('mat1');
fe.props.h=1;
material.setElasticIzo(1, 0.3);

% Assigning material to finite element
fe.setMaterial( material );

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
analysis = LinearElasticityWeighted( fe, mesh, true );

% Creating node selector object to select fixed edge (left)
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.001 ) );

% Fixing structure according to above defined node selector object
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );

% Creating load vector with one node loaded at the middle of right edge
analysis.loadClosestNode([aspect*h, h/2 ], ["ux" "uy"], [0 -1] );

% Defining an object for stress intensity based topology optimization with volume constraint.
tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );

% Executing the above mentioned algorithm 
[objF, xopt]  = topOpt.solve();
toc

figure;
tic

% Defining an object for MMA based compliance minimization topology optimization algorithm with volume constraint.
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, true);

% Executing the above mentioned algorithm 
[objF, xopt]  = topOpt.solve();
toc



