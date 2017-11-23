%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Euler Bernoulli eigenfrequencies analytical solution
%%% Arthur Schout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% free free eigen frequencies for the first modes

E = 210e9;          % [Pa] Youngs modulus
rho = 7.8e3;         % [Kg*m^3]  density of the material    


% Beam dimensions
h = 0.1;            % height
w = 0.1;            % width

A = h*w;            % Area of the beam
I = 1/12*h^3*w;     % moment of inertia

L_beam = 1;
n_int = [ 1.50562  2.49975  3.50001 4.50000];

for i = 1:length(n_int)
    Beta(i) = n_int(i) * pi / L_beam;
    Freq(i) = (Beta(i)^2*sqrt(E*I/(A*rho)));
end

Freq

omega_1 = 22.3733/L_beam^2*sqrt(E*I/(A*rho))
figure(20)
hold on
plot([1:6],[ 0 0 Freq ],'--x')
plot([1:6],freq(1:6))
title('Free-free eigenfrequencies')
legend('analytical','fem')

%% clamped free eigen frequencies for the first modes

n_int_2 = [ 0.59686   1.49418 2.50025 3.49999];

for i = 1:length(n_int_2)
    Beta(i) = n_int_2(i) * pi / L;
    Freq(i) = (Beta(i)^2*sqrt(E*I/(A*rho)));
end

Freq

hold on
plot([1:6],[ 0 0 Freq ],'--x')
plot([1:6],freq(1:6))
title('Clamped free eigenfrequencies')
legend('analytical','fem')