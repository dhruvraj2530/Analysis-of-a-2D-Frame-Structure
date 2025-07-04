







% Initialize workspace
clear; clc;

% Define scale factor for plotting
applied_scale_factor = 2000;

% Specify restrained degrees of freedom (DOFs)
restrained_dofs = [205, 206, 260, 301, 302, 311, 323, 338, 344, 346, 347, 386, 388, 389];


% Read joint coordinates and member connectivity from Excel
Node_Coordinates = xlsread("DRSR.xlsx", "Joint Coordinates1");
Member_Connectivity = xlsread("DRSR.xlsx", "Member Connectivity1");

% Get the number of joints and members
[num_joints, ~] = size(Node_Coordinates);
[num_members, ~] = size(Member_Connectivity);






% Initialize global stiffness matrix and force vectors
K_global = zeros(3 * num_joints);
fixed_end_forces = zeros(3 * num_joints, 1);
nodal_forces = zeros(3 * num_joints, 1);
settlements = zeros(3 * num_joints, 1);

% Assemble global stiffness matrix and load vectors
for member = 1:num_members
    % Extract member properties
    area = Member_Connectivity(member, 4);
    elasticity = Member_Connectivity(member, 5);
    inertia = Member_Connectivity(member, 6);
    load_w = Member_Connectivity(member, 7);
    joint_case = Member_Connectivity(member, 8);
    node_i = Member_Connectivity(member, 2);
    node_j = Member_Connectivity(member, 3);
    
    % Extract settlements
    settlement_x_i = Node_Coordinates(node_i, 4);
    settlement_y_i = Node_Coordinates(node_i, 5);
    settlement_x_j = Node_Coordinates(node_j, 4);
    settlement_y_j = Node_Coordinates(node_j, 5);
    
    % Extract concentrated forces and moments
    conc_force_x_i = Node_Coordinates(node_i, 6);
    conc_force_y_i = Node_Coordinates(node_i, 7);
    conc_moment_z_i = Node_Coordinates(node_i, 8);
    conc_force_x_j = Node_Coordinates(node_j, 6);
    conc_force_y_j = Node_Coordinates(node_j, 7);
    conc_moment_z_j = Node_Coordinates(node_j, 8);
    
    % Coordinates of nodes
    x_i = Node_Coordinates(node_i, 2);
    y_i = Node_Coordinates(node_i, 3);
    x_j = Node_Coordinates(node_j, 2);
    y_j = Node_Coordinates(node_j, 3);
    
    % Calculate member length and orientation
    member_length = sqrt((x_j - x_i)^2 + (y_j - y_i)^2);
    cos_theta = (x_j - x_i) / member_length;
    sin_theta = (y_j - y_i) / member_length;
    
    % Determine local stiffness matrix and fixed-end forces based on joint case
    switch joint_case
        case 1
            kloc = ((inertia * elasticity) / member_length^3) * ...
                [(area * member_length^2) / inertia, 0, 0, -(area * member_length^2) / inertia, 0, 0;
                 0, 12, 6 * member_length, 0, -12, 6 * member_length;
                 0, 6 * member_length, 4 * member_length^2, 0, -6 * member_length, 2 * member_length^2;
                 -(area * member_length^2) / inertia, 0, 0, (area * member_length^2) / inertia, 0, 0;
                 0, -12, -6 * member_length, 0, 12, -6 * member_length;
                 0, 6 * member_length, 2 * member_length^2, 0, -6 * member_length, 4 * member_length^2];
            Q = [0; (load_w * member_length) / 2; (load_w * member_length^2) / 12; 0; (load_w * member_length) / 2; -(load_w * member_length^2) / 12];
        case 2
            kloc = ((inertia * elasticity) / member_length^3) * ...
                [(area * member_length^2) / inertia, 0, 0, -(area * member_length^2) / inertia, 0, 0;
                 0, 3, 0, 0, -3, 3 * member_length;
                 0, 0, 0, 0, 0, 0;
                 -(area * member_length^2) / inertia, 0, 0, (area * member_length^2) / inertia, 0, 0;
                 0, -3, 0, 0, 3, -3 * member_length;
                 0, 3 * member_length, 0, 0, -3 * member_length, 3 * member_length^2];
            Q = [0; (load_w * member_length) / 2 - (1.5 / member_length) * (load_w * member_length^2) / 12; 0; 0; (load_w * member_length) / 2 + (1.5 / member_length) * (load_w * member_length^2) / 12; -(load_w * member_length^2) / 12 - 0.5 * (load_w * member_length^2) / 12];
        case 3
            kloc = ((inertia * elasticity) / member_length^3) * ...
                [(area * member_length^2) / inertia, 0, 0, -(area * member_length^2) / inertia, 0, 0;
                 0, 3, 3 * member_length, 0, -3, 0;
                 0, 3 * member_length, 3 * member_length^2, 0, -3 * member_length, 0;
                 -(area * member_length^2) / inertia, 0, 0, (area * member_length^2) / inertia, 0, 0;
                 0, -3, -3 * member_length, 0, 3, 0;
                 0, 0, 0, 0, 0, 0];
            Q = [0; (load_w * member_length) / 2 - (1.5 / member_length) * (-load_w * member_length^2) / 12; (load_w * member_length^2) / 12 - 0.5 * (-load_w * member_length^2) / 12; 0; (load_w * member_length) / 2 + (1.5 / member_length) * (-load_w * member_length^2) / 12; 0];
        case 4
            kloc = (elasticity * area) / member_length * ...
                [1, 0, 0, -1, 0, 0;
                 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0;
                 -1, 0, 0, 1, 0, 0;
                 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0];
            Q = [0; (load_w * member_length) / 2; 0; 0; (load_w * member_length) / 2; 0];
        otherwise
            error('Invalid joint case specified.');
    end
    
    % Transformation matrix
    T = [cos_theta, sin_theta, 0, 0, 0, 0;
        -sin_theta, cos_theta, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0;
         0, 0, 0, cos_theta, sin_theta, 0;
         0, 0, 0, -sin_theta, cos_theta, 0;
         0, 0, 0, 0, 0, 1];
    
    % Global stiffness matrix for the member
    k_global_member = T' * kloc * T;
    
    % Global DOF indices
    global_dof = [3 * (node_i - 1) + 1, 3 * (node_i - 1) + 2, 3 * (node_i - 1) + 3, 3 * (node_j - 1) + 1, 3 * (node_j - 1) + 2, 3 * (node_j - 1) + 3];
    
    % Assemble into global stiffness matrix
    K_global(global_dof, global_dof) = K_global(global_dof, global_dof) + k_global_member;
    
    % Account for settlements
    delta = zeros(6, 1);
    delta(1) = settlement_x_i;
    delta(2) = settlement_y_i;
    delta(4) = settlement_x_j;
    delta(5) = settlement_y_j;
    settlements(global_dof([1, 2, 4, 5])) = [settlement_x_i; settlement_y_i; settlement_x_j; settlement_y_j];
    z = k_global_member * delta;
    fixed_end = T' * Q + z;
    
    % Assemble fixed-end forces
    fixed_end_forces(global_dof) = fixed_end_forces(global_dof) + fixed_end;
    
    % Assemble nodal forces
    nodal_forces(global_dof) = nodal_forces(global_dof) + [conc_force_x_i; conc_force_y_i; conc_moment_z_i; conc_force_x_j; conc_force_y_j; conc_moment_z_j];
end

% Apply boundary conditions by removing restrained DOFs
K_reduced = K_global;
nodal_forces_reduced = nodal_forces;
fixed_end_forces_reduced = fixed_end_forces;

for idx = 1:length(restrained_dofs)
    dof = restrained_dofs(idx) - (idx - 1);
    K_reduced(dof, :) = [];
    K_reduced(:, dof) = [];
    nodal_forces_reduced(dof) = [];
    fixed_end_forces_reduced(dof) = [];
end

% Solve for displacements
displacements_reduced = K_reduced \ (nodal_forces_reduced - fixed_end_forces_reduced);

% Reconstruct full displacement vector including restrained DOFs
displacements_full = zeros(3 * num_joints, 1);
free_dofs = setdiff(1:3 * num_joints, restrained_dofs);
displacements_full(free_dofs) = displacements_reduced;
displacements_full = displacements_full + settlements;

% Write member forces to Excel
headers = ["MEMBER_NUMBER", "F_x_NODEi", "F_y_NODEi", "M_z_NODEi", "F_x_NODEj", "F_y_NODEj", "M_z_NODEj", "Length", "Uniformly_distributed_Load"];
writematrix(headers, 'Output DRSR.xlsx', 'Range', 'A1');

for member = 1:num_members
    % Extract member properties
    area = Member_Connectivity(member, 4);
    elasticity = Member_Connectivity(member, 5);
    inertia = Member_Connectivity(member, 6);
    load_w = Member_Connectivity(member, 7);
    joint_case = Member_Connectivity(member, 8);
    node_i = Member_Connectivity(member, 2);
    node_j = Member_Connectivity(member, 3);
    
    % Coordinates of nodes
    x_i = Node_Coordinates(node_i, 2);
    y_i = Node_Coordinates(node_i, 3);
    x_j = Node_Coordinates(node_j, 2);
    y_j = Node_Coordinates(node_j, 3);
    
    % Calculate member length and orientation
    member_length = sqrt((x_j - x_i)^2 + (y_j - y_i)^2);
    cos_theta = (x_j - x_i) / member_length;
    sin_theta = (y_j - y_i) / member_length;
    
    % Determine local stiffness matrix and fixed-end forces based on joint case
    switch joint_case
        case 1
            kloc = ((inertia * elasticity) / member_length^3) * ...
                [(area * member_length^2) / inertia, 0, 0, -(area * member_length^2) / inertia, 0, 0;
                 0, 12, 6 * member_length, 0, -12, 6 * member_length;
                 0, 6 * member_length, 4 * member_length^2, 0, -6 * member_length, 2 * member_length^2;
                 -(area * member_length^2) / inertia, 0, 0, (area * member_length^2) / inertia, 0, 0;
                 0, -12, -6 * member_length, 0, 12, -6 * member_length;
                 0, 6 * member_length, 2 * member_length^2, 0, -6 * member_length, 4 * member_length^2];
            Q = [0; (load_w * member_length) / 2; (load_w * member_length^2) / 12; 0; (load_w * member_length) / 2; -(load_w * member_length^2) / 12];
        case 2
            kloc = ((inertia * elasticity) / member_length^3) * ...
                [(area * member_length^2) / inertia, 0, 0, -(area * member_length^2) / inertia, 0, 0;
                 0, 3, 0, 0, -3, 3 * member_length;
                 0, 0, 0, 0, 0, 0;
                 -(area * member_length^2) / inertia, 0, 0, (area * member_length^2) / inertia, 0, 0;
                 0, -3, 0, 0, 3, -3 * member_length;
                 0, 3 * member_length, 0, 0, -3 * member_length, 3 * member_length^2];
            Q = [0; (load_w * member_length) / 2 - (1.5 / member_length) * (load_w * member_length^2) / 12; 0; 0; (load_w * member_length) / 2 + (1.5 / member_length) * (load_w * member_length^2) / 12; -(load_w * member_length^2) / 12 - 0.5 * (load_w * member_length^2) / 12];
        case 3
            kloc = ((inertia * elasticity) / member_length^3) * ...
                [(area * member_length^2) / inertia, 0, 0, -(area * member_length^2) / inertia, 0, 0;
                 0, 3, 3 * member_length, 0, -3, 0;
                 0, 3 * member_length, 3 * member_length^2, 0, -3 * member_length, 0;
                 -(area * member_length^2) / inertia, 0, 0, (area * member_length^2) / inertia, 0, 0;
                 0, -3, -3 * member_length, 0, 3, 0;
                 0, 0, 0, 0, 0, 0];
            Q = [0; (load_w * member_length) / 2 - (1.5 / member_length) * (-load_w * member_length^2) / 12; (load_w * member_length^2) / 12 - 0.5 * (-load_w * member_length^2) / 12; 0; (load_w * member_length) / 2 + (1.5 / member_length) * (-load_w * member_length^2) / 12; 0];
        case 4
            kloc = (elasticity * area) / member_length * ...
                [1, 0, 0, -1, 0, 0;
                 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0;
                 -1, 0, 0, 1, 0, 0;
                 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0];
            Q = [0; (load_w * member_length) / 2; 0; 0; (load_w * member_length) / 2; 0];
        otherwise
            error('Invalid joint case specified.');
    end
    
    % Transformation matrix
    T = [cos_theta, sin_theta, 0, 0, 0, 0;
        -sin_theta, cos_theta, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0;
         0, 0, 0, cos_theta, sin_theta, 0;
         0, 0, 0, -sin_theta, cos_theta, 0;
         0, 0, 0, 0, 0, 1];
    
    % Global DOF indices
    global_dof = [3 * (node_i - 1) + 1, 3 * (node_i - 1) + 2, 3 * (node_i - 1) + 3, 3 * (node_j - 1) + 1, 3 * (node_j - 1) + 2, 3 * (node_j - 1) + 3];
    
    % Extract displacements for member DOFs
    displacements_member = displacements_full(global_dof);
    
    % Compute member end forces
    member_end_forces = kloc * (T * displacements_member) + Q;
    
    % Prepare data for writing
    member_forces_output = round(0.001 * member_end_forces, 2);
    data_row = {member, member_forces_output(1), member_forces_output(2), member_forces_output(3), member_forces_output(4), member_forces_output(5), member_forces_output(6), member_length, load_w};
    
    % Write data to Excel
    writecell(data_row, 'Output DRSR.xlsx', 'Range', strcat('A', num2str(member + 1)));
end

% Plot undeformed and deformed structures
x_coords_undeformed = Node_Coordinates(:, 2);
y_coords_undeformed = Node_Coordinates(:, 3);
x_coords_deformed = x_coords_undeformed + displacements_full(1:3:end) * applied_scale_factor;
y_coords_deformed = y_coords_undeformed + displacements_full(2:3:end) * applied_scale_factor;

figure;
hold on;
for member = 1:num_members
    node_i = Member_Connectivity(member, 2);
    node_j = Member_Connectivity(member, 3);
    
    % Plot undeformed member
    plot([x_coords_undeformed(node_i), x_coords_undeformed(node_j)], [y_coords_undeformed(node_i), y_coords_undeformed(node_j)], 'r');
    
    % Plot deformed member
    plot([x_coords_deformed(node_i), x_coords_deformed(node_j)], [y_coords_deformed(node_i), y_coords_deformed(node_j)], 'g');
end
hold off;
legend('Undeformed', 'Deformed');
title('Structural Deformation');
xlabel('X Coordinate');
ylabel('Y Coordinate');
disp('The deformed and undeformed structures have been plotted. The deformation is scaled for visibility.');
fprintf('The deformed shape has been exaggerated %d times for visualization purposes.\n', applied_scale_factor);

