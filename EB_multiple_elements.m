% euler bernoulli beam implementation in matlab
% AS 22/11/17
% 1 dimensional
%%% Preprocessing


Clear all           % clear the memory

%% Variables

%constants
E = 2e9;            % youngs modulus

%beam dimensions
h = 0.1;             % height
w = 0.1;            % width

A = h*w;            % Area of the beam
I = 1/12*h^3*w;     % moment of inertia

% node locations already in 3d coordinates only two nodes
n = [ 0 0 0; 1 0 0 ];

dn = n(2,:)-n(1,:);

l_element = norm(dn);

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
    fprintf('!!!\n!!!\tElement matrix is not symmetric\n!!!\n');
end

% constrained DOF's
constrained_dof = [1 2];

%force vector
force = [ 0 0 100 0 ];

% set rows and columns to zero where BC's are applied in stiffnes matrix
K_element_sol = K_element;
K_element_sol(constrained_dof,:) = 0;
K_element_sol(:,constrained_dof) = 0;
K_element_sol(constrained_dof,constrained_dof) = eye(length(constrained_dof));

displacement = K_element_sol\force';

K_element*displacement 

%% plot solution


%cubic spline interpolation to get shape of beam element 
pp = csape(n(:,1),[displacement(2) displacement(1) displacement(3) displacement(4)]);

clf
figure(1)
hold on
plot(n(:,1),n(:,2),'--xb')
plot(n(:,1),[displacement(1) displacement(3)],'-or')
fnplt(pp)








%solve the problem



