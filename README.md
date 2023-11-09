# MorphoGen

Computational MorphoGenesis System

## Authors
[Piotr Tauzowski](mailto:ptauzow@ippt.pan.pl)

## Introduction
__MorphoGen__ system is used to solve tasks in the field of structural topology optimization. The system is written in MATLAB using an object-oriented programming (OOP) paradigm and a layered architecture, both of which allow efficient development of the system. The system also includes two main layers: the first one is responsible for finite element analysis (FEA) and the second one for reliability assesment. FEA is an essential component of structural topology optimization. Due to OOP implenetation new types of finite elements and analyzes can be easily added. Main directory of the system consists of several subdirectories with source code described below:

* __examples__ - here users find ready-to-use example scripts showcasing the implementation and use of the software.
* __design__ - layer implements algorithms for topology optimization such as derivative-free approach called stress intensity driven topology optimization, as well as sensitivity-based approach known as Method of
Moving Asymptotes.
* __analysis__ layer contains algorithms dedicated to solving specific problems for e.g.: algorithms for topology optimization such as derivative-free approach called stress intensity driven topology optimization, as well as sensitivity-based approach known as Method of Moving Asymptotes. This directory contains also finite element algorithms for example: linear elastic analysis with sensitivity, elastoplastic analy-
sis, etc. Additionally, reliability algorithms such as FORM (First-Order Reliability Method), AMV (Advanced Mean Value), and Monte Carlo simulation are stored here.
* __element__ stores the finite element files, including planar and spatial element definitions.
* __material__ layer contains material definitions specific to each finite element type, e.g. isotropic materials are defined for plane stress/strain and solid elements.
* __mesh__ contains nodes and elements organized into mesh class designed for finite element mesh storage and generation.
* __math__ provides various mathematical utilities, including Gauss numerical integration and linear equation solver class. They facilitate
data manipulation of global matrices and right-hand vectors. Additionally, shape function classes compatible with finite elements are included.
These directories provide a comprehensive framework to explore and utilize the software’s features and capabilities.

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
% Cantilever topology optimization elastic task
clear; close all;

% Initialization 
res = 50; % resolution of shortest (vertical) edge
h = 1; % height of the cantilever
aspect = 2; % aspect ratio length/height
Rfilter = 2*h/res; % filtering radius
cutTreshold = 0.005; % paramter selecting the removal intensity threshold
penal = 3; % penalty factor

% MESH layer (creating FE mesh object)
mesh = Mesh();
mesh.addRectMesh2D(0, 0, aspect*h, h, aspect*res, res, sfL4.pattern); % generating rectangular mesh ( aspect times h by h )

% MATERIAL layer (creating isotropic material object)
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
fe.setMaterial( material ); % assigning material to finite element

% ELEMENT layer (creating node selector object to select fixed edge (left))
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.001 ) );
sfL4 = ShapeFunctionL4; % type of shape function to be used (here: four node Langrange)
fe = PlaneStressElem( sfL4, mesh.elems ); % create plane stress finite element object

% ANALYSIS layer (creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density)
analysis = LinearElasticityWeighted( fe, mesh, true );
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] ); % fixing structure according to above defined node selector object
analysis.loadClosestNode([aspect*h, h/2 ], ["ux" "uy"], [0 -1] ); % creating load vector with one node loaded at the middle of right edge

% DESIGN layer (defining an object for stress intensity based topology optimization with volume constraint)
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve(); % executing the above mentioned algorithm 
```
