function u_solid = applyRigidBodyKinematics(solid_nodes, frame_node, q_frame6)
% APPLYRIGIDBODYKINEMATICS  Compute solid-node displacements from a 6-DOF
% frame node using linearised rigid-body kinematics.
%
%   u(p) = u(c) + theta(c) x (p - c)
%
% This implements the Euler-Bernoulli cross-section hypothesis at the segment
% connection level: every node on the cross-section undergoes the rigid-body
% motion defined by the frame-node translation and rotation.  The formula is
% exact for rigid bodies and is the leading-order approximation for slender
% elastic members (small displacement / small rotation).
%
% INPUTS
%   solid_nodes  - N x 3  coordinates of solid nodes on the cross-section
%   frame_node   - 1 x 3  position of the frame node (cross-section centroid)
%   q_frame6     - 1 x 6  [ux, uy, uz, theta_x, theta_y, theta_z] at the frame node
%
% OUTPUT
%   u_solid      - N x 3  prescribed displacements [ux, uy, uz] for each node

u_c   = q_frame6(1:3);      % translational DOFs
theta = q_frame6(4:6);      % rotational DOFs (small-rotation vector)

N  = size(solid_nodes, 1);
dp = solid_nodes - frame_node;            % N x 3 relative position vectors

% cross(theta, dp)  row-wise: theta x (p_i - c)
u_solid = repmat(u_c, N, 1) + cross(repmat(theta, N, 1), dp, 2);
end
