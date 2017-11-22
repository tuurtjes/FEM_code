% euler bernoulli beam implementation in matlab
% AS 22/11/17
% 1 dimensional
%%% Preprocessing


clear all           % clear the memory

%% Variables

%constants
E = 210e9;            % youngs modulus

%beam dimensions
h = 0.1;             % height
w = 0.1;            % width

A = h*w;            % Area of the beam
I = 1/12*h^3*w;     % moment of inertia

% Number of elements in between the nodes
n_div = 1;

% Points at the endpoints of the structure
n = [ 0 0 0; 1 0 0 ];

%check length in between the points
dn = n(2,:)-n(1,:);
%calculate the length of an element
l_element = norm(dn)/n_div;

%nodal degrees of freedom
nodal_dof = 2;

%define matrix with node coordinates
n_nodes = n_div+1;%number of nodes
nodes = zeros(n_nodes,3);%matrix for the coordinates of the nodes
for i = 0:n_div
    % x-coordinates of the nodes
    nodes(i+1,1) = i*l_element;
    % not implemented yet
    % y-coordinates of the nodes
    % z-coordinates of the nodes
end

% nodal connectivity in this case it can be generated fairly easily
nodal_con = [ 1:n_div ; 2:n_div+1 ]';
% nodal_con(end,1) = 1;

% no matrix transforamtions are needed due to the 1 dimensionality of the
% problem. 
% Local coordinates are equal to the global coordinates

% an euler bernoulli beam element contains three 'basic' matrix elements
K_1 = E*I/l_element^3;
K_2 = E*I/l_element^2;
K_3 = E*I/l_element;

% definition of the element matrix from theory
K_element = [   12*K_1  6*K_2   -12*K_1 6*K_2 ; ...
                6*K_2   4*K_3   -6*K_2  2*K_3 ; ...
                -12*K_1 -6*K_2  12*K_1  -6*K_2; ...
                6*K_2   2*K_3   -6*K_2  4*K_3 ];
                
% check if element matrix is symmetric if not print

symmetry = isequal(K_element,K_element');

if symmetry ~= 1
    fprintf('!!!\n!!!\tElement stiffness matrix is not symmetric\n!!!\n');
end


%% Assembling the global stiffness matrix

%size of the global stiffness matrix is equal to the number of degrees of
%freedom of the model. n_dof = n_nodes * nodal_dof
n_dof   = n_nodes * nodal_dof;
K_global= zeros(n_dof); % initialize the stiffness matrix

K_part = mat2cell(K_element,[nodal_dof nodal_dof],[nodal_dof nodal_dof]);

for i = 1:size(nodal_con,1)
    con_1 = nodal_con(i,1); %get the nodal connectivity information
    con_2 = nodal_con(i,2); %for each element
    
    %fill up the global stiffness matrix with partitioned element
    %stiffness matrices in 4 steps
    i11 = con_1*2
    i22 = con_2*2
    %top left block of the partitioned matrix
    K_global(i11-1:i11,i11-1:i11) = K_global(i11-1:i11,i11-1:i11) +...
                                    K_part{1,1};
    %bottom right block
    K_global(i22-1:i22,i22-1:i22) = K_global(i22-1:i22,i22-1:i22) +...
                                    K_part{2,2};
    %bottom left
	K_global(i22-1:i22,i11-1:i11) = K_global(i22-1:i22,i11-1:i11) +...
                                    K_part{2,1};
    %top right
	K_global(i11-1:i11,i22-1:i22) = K_global(i11-1:i11,i22-1:i22) +...
                                    K_part{1,2};
end


symmetry = isequal(K_global,K_global');

if symmetry ~= 1
    fprintf('!!!\n!!!\tElement stiffness matrix is not symmetric\n!!!\n');
end


%% Defining the boundary conditions

% constrained DOF's
constrained_dof = [1 2];

%force vector
force = zeros(1,n_dof);
force(end-1) = 1000;

% set rows and columns to zero where BC's are applied in stiffnes matrix
K_element_sol = K_global;
K_element_sol(constrained_dof,:) = 0;
K_element_sol(:,constrained_dof) = 0;
K_element_sol(constrained_dof,constrained_dof) = eye(length(constrained_dof));

%% processing

%find the unknown diplacements
displacement = K_element_sol\force';

%find the unknown forces
forces = K_global*displacement;

%% Post processing

%cubic spline interpolation to get shape of beam element 
pp = csape(n(:,1),[displacement(2) displacement(1) displacement(end-1) displacement(end)]);

%clf
figure(n_div)
hold on
plot(nodes(:,1),nodes(:,2),'--xb')
plot(nodes(:,1),[displacement(1:2:end-1)],'-or')
fnplt(pp)



