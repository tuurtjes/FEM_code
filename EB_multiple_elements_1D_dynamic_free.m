% euler bernoulli beam implementation in matlab
% AS 22/11/17
% 1 dimensional
% Equidistan mesh refinement is implemented in between two points

%%%%%%%%%%%%%%%%%%%%%
%%% Preprocessing %%%
%%%%%%%%%%%%%%%%%%%%%

%clear all           % clear the memory

%% Variables

% Material properties
E   = 210e9;         % [Pa] Youngs modulus
rho = 7.8e3;         % [Kg*m^3]  density of the material    


% Beam dimensions
h = 0.1;            % height
w = 0.1;            % width

A = h*w;            % Area of the beam
I = 1/12*h^3*w;     % moment of inertia

% Number of elements in between the nodes
n_div = 5;

% Points at the endpoints of the structure
n = [ 0 0 0; 1 0 0 ];

% check length in between the points
dn = n(2,:)-n(1,:);

% calculate the length of an element
L_e = norm(dn)/n_div;

% nodal degrees of freedom
nodal_dof = 2;

% define matrix with node coordinates
n_nodes = n_div+1;%number of nodes
nodes = zeros(n_nodes,3);%matrix for the coordinates of the nodes
for i = 0:n_div
    % x-coordinates of the nodes
    nodes(i+1,1) = i*L_e;
    % not implemented yet
    % y-coordinates of the nodes
    % z-coordinates of the nodes
end

% nodal connectivity in this case it can be generated fairly easily
nodal_con = [ 1:n_div ; 2:n_div+1 ]';
% nodal_con(end,1) = 1;

% No matrix transforamtions are needed due to the 1 dimensionality of the
% problem. 
% So far the Local coordinates are equal to the global coordinates

L = L_e;

% definition of the element stiffnes matrix
K_element = E*I/(L^3) * ...
    [   12      6*L     -12     6*L     ; ...
        6*L     4*L^2   -6*L    2*L^2   ; ...
        -12     -6*L    12      -6*L    ; ...
        6*L     2*L^2   -6*L    4*L^2   ];

            
% definition of the element mass matrix 
M_element = rho*A*L/420* ...
    [   156     22*L    54      -13*L   ;...
        22*L    4*L^2   13*L    -3*L^2  ;...
        54      13*L    156     -22*L   ;...
        -13*L   -3*L^2  -22*L   4*L^2   ];

% some sanity checks on the matrices    
% check if element matrices are symmetric and print warning if not
sym_K_el = isequal(K_element,K_element');
sym_M_el = isequal(M_element,M_element');

if sym_K_el ~= 1
    fprintf('!!!\n!!!\tElement stiffness matrix is not symmetric\n!!!\n');
end
if sym_M_el ~= 1
    fprintf('!!!\n!!!\tElement mass matrix is not symmetric\n!!!\n');
end

% check if element matrices are positive definite
[ ~ , PSD_K] = chol(K_element);
[ ~ , PSD_M] = chol(M_element);

if PSD_K ~= 0
    fprintf('!!!\n!!!\tElement Stiffness matrix is not PSD\n!!!\n')
end

if PSD_M ~= 0
    fprintf('!!!\n!!!\tElement Mass matrix is not PSD\n!!!\n')
end



%% Assembling the global matrices

%%% STIFFNESS MATRIX %%%
% size of the global stiffness matrix is equal to the number of degrees of
% freedom of the model. n_dof = n_nodes * nodal_dof
n_dof   = n_nodes * nodal_dof;
K_global= zeros(n_dof); % initialize the stiffness matrix

% Partition the element stiffness matrix 
% !!! WHEN ASSEMBLING A STRUCTURE WHERE THE LOCAL AND GLOBAL COORDINATES
% DIFFER THIS MATRIX NEEDS TO BE ADAPTED FOR EACH ELEMENT ORIENTATION
K_part = mat2cell(K_element,[nodal_dof nodal_dof],[nodal_dof nodal_dof]);

% Actual assembly of K global
for i = 1:size(nodal_con,1)
    con_1 = nodal_con(i,1); %get the nodal connectivity information
    con_2 = nodal_con(i,2); %for each element
    
    % Fill up the global stiffness matrix with partitioned element
    % Stiffness matrices in 4 steps
    i11 = con_1*2;
    i22 = con_2*2;
    % Top left block of the partitioned matrix
    K_global(i11-1:i11,i11-1:i11) = K_global(i11-1:i11,i11-1:i11) +...
                                    K_part{1,1};
    % Bottom right block
    K_global(i22-1:i22,i22-1:i22) = K_global(i22-1:i22,i22-1:i22) +...
                                    K_part{2,2};
    % Bottom left block
	K_global(i22-1:i22,i11-1:i11) = K_global(i22-1:i22,i11-1:i11) +...
                                    K_part{2,1};
    % Top right block
	K_global(i11-1:i11,i22-1:i22) = K_global(i11-1:i11,i22-1:i22) +...
                                    K_part{1,2};
end

symmetry = isequal(K_global,K_global');

if symmetry ~= 1
    fprintf('!!!\n!!!\tGlobal stiffness matrix is not symmetric\n!!!\n');
end

%%% MASS MATRIX %%%
% The size of the global stiffness matrix is equal to the number of degrees of
% freedom of the model. n_dof = n_nodes * nodal_dof
M_global= zeros(n_dof); % initialize the stiffness matrix
% Partition the element mass matrix 
% !!! WHEN ASSEMBLING A STRUCTURE WHERE THE LOCAL AND GLOBAL COORDINATES
% DIFFER THIS MATRIX NEEDS TO BE ADAPTED FOR EACH ELEMENT ORIENTATION
M_part = mat2cell(M_element,[nodal_dof nodal_dof],[nodal_dof nodal_dof]);

% Actual assembly of M global !!! COULD BE INTERATED IN THE K LOOP ASWELL
for i = 1:size(nodal_con,1)
    con_1 = nodal_con(i,1); %get the nodal connectivity information
    con_2 = nodal_con(i,2); %for each element
    
    % Fill up the global stiffness matrix with partitioned element
    % Stiffness matrices in 4 steps
    i11 = con_1*2;
    i22 = con_2*2;
    % Top left block of the partitioned matrix
    M_global(i11-1:i11,i11-1:i11) = M_global(i11-1:i11,i11-1:i11) +...
                                    M_part{1,1};
    % Bottom right block
    M_global(i22-1:i22,i22-1:i22) = M_global(i22-1:i22,i22-1:i22) +...
                                    M_part{2,2};
    % Bottom left block
	M_global(i22-1:i22,i11-1:i11) = M_global(i22-1:i22,i11-1:i11) +...
                                    M_part{2,1};
    % Top right block
	M_global(i11-1:i11,i22-1:i22) = M_global(i11-1:i11,i22-1:i22) +...
                                    M_part{1,2};
end

symmetry = isequal(M_global,M_global');

if symmetry ~= 1
    fprintf('!!!\n!!!\tGlobal mass matrix is not symmetric\n!!!\n');
end


[ ~ , PSD_Kg] = chol(K_global);
[ ~ , PSD_Mg] = chol(M_global);

if PSD_Kg ~= 0
    fprintf('!!!\n!!!\tGlobal Stiffness matrix is not PSD\n!!!\n')
end

if PSD_Mg ~= 0
    fprintf('!!!\n!!!\tGlobal Mass matrix is not PSD\n!!!\n')
end

%% Defining the boundary conditions for the Eigenvalue problem
%beam is free in space, OK for the simple eigenvalue problem
% Constrained DOF's
%constrained_dof = [1 2];

% Force vector
%force = zeros(1,n_dof);
%force(end-1) = 1000;

% Set rows and columns to zero where BC's are applied in stiffnes matrix
% K_element_sol = K_global;
% K_element_sol(constrained_dof,:) = 0;
% K_element_sol(:,constrained_dof) = 0;
% K_element_sol(constrained_dof,constrained_dof) = eye(length(constrained_dof));



%%%%%%%%%%%%%%%%%%
%%% Processing %%%
%%%%%%%%%%%%%%%%%%

%% solving the eigenvalue problem

% Find the unknown diplacements
%displacement = K_element_sol\force';

% Find the unknown forces
%forces = K_global*displacement;

%Eigen_matrix = K_global\M_global;

%[V,D] = eig(Eigen_matrix);
[U,C] = eig(K_global,M_global);
%%%%%%%%%%%%%%%%%%%%%%%
%%% Post processing %%%
%%%%%%%%%%%%%%%%%%%%%%%

freq = sqrt(diag(C));
%plot([1:6],freq(1:6))



%% plotting the results

figure(2)
for i = 3:6
    %figure(i)
    hold on
    plot(nodes(:,1),U(1:2:end,i),'--x')
    title(strcat('bending mode ', num2str(i-2))) 
end


% % Cubic spline interpolation to get shape of beam element 
% pp = csape(n(:,1),[displacement(2) displacement(1) displacement(end-1) displacement(end)]);
% 
% %clf
% figure(n_div)
% hold on
% plot(nodes(:,1),nodes(:,2),'--xb')
% plot(nodes(:,1),[displacement(1:2:end-1)],'-or')
% fnplt(pp)



