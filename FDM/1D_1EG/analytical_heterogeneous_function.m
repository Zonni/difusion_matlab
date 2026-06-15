%% Analytical Solution: Heterogeneous 1D Slab Reactor (One Energy Group)
% Translation from Python for Nuclear Engineering Internship
clear; clc; close all;

% ==========================================
% 1. Define Material Properties & Dimensions
% ==========================================
N = 100;           % Number of points for plotting
% Core properties
Dc = 1.0;           % Diffusion coefficient in core (cm)
Sigma_ac = 0.05;    % Absorption cross-section in core (cm^-1)
nu_Sigma_f = 0.06;  % Production cross-section in core (cm^-1)

% Reflector properties
Dr = 1.0;           % Diffusion coefficient in reflector (cm)
Sigma_ar = 0.02;    % Absorption cross-section in reflector (cm^-1)

% Dimensions
a = 150.0;          % Core half-width (cm)
b = 25.0;           % Reflector thickness (cm)

% ==========================================
% 2. Reflector Parameters
% ==========================================
% Inverse diffusion length in the reflector
kappa_r = sqrt(Sigma_ar / Dr);

% ==========================================
% 3. Solve Criticality Condition for k_eff
% ==========================================
% Initial guess using Bare Reactor
d = 2.13 * Dc;      % Extrapolation distance
a_ext = a + d;
Bg_bare = pi / (2 * a_ext);
k_bare = nu_Sigma_f / (Dc * Bg_bare^2 + Sigma_ac);

k_inf = nu_Sigma_f / Sigma_ac;

fprintf('Initial guess (Bare k_eff): %.5f\n', k_bare);

% Define the transcendental equation as an anonymous function
% We want to find k such that f(k) = 0
criticality_eq = @(k) (Dc * sqrt((nu_Sigma_f/k - Sigma_ac)/Dc) * ...
                tan(sqrt((nu_Sigma_f/k - Sigma_ac)/Dc) * a)) - ...
                (Dr * kappa_r * (1 / tanh(kappa_r * b)));

% Use fzero (MATLAB's equivalent to fsolve for scalar equations)
options = optimset('Display','off');
k_eff = fzero(criticality_eq, [k_bare, k_inf - 1e-5], options);

fprintf('--- Reactor Parameters ---\n');
fprintf('Fixed core half-width (a) : %.2f cm\n', a);
fprintf('Fixed reflector width (b) : %.2f cm\n', b);
fprintf('Calculated k_eff          : %.7f\n\n', k_eff);

% ==========================================
% 4. Calculate Final Buckling and Flux Constants
% ==========================================
Bc = sqrt((nu_Sigma_f / k_eff - Sigma_ac) / Dc);

% Flux normalization
A = 1.0; % Arbitrary flux at center (x=0)
% Continuity at x = a: A*cos(Bc*a) = C*sinh(kappa_r*(a+b-a)) -> C*sinh(kappa_r*b)
C = A * cos(Bc * a) / sinh(kappa_r * b);

% ==========================================
% 5. Generate Flux Vector
% ==========================================
% Define x from -(a+b) to (a+b)
x_plot = linspace(-(a + b), a + b, N);
phi = zeros(size(x_plot));

for i = 1:N
    abs_x = abs(x_plot(i));
    if abs_x <= a
        phi(i) = A * cos(Bc * abs_x);
    else
        % Flux in reflector: sinh(kappa_r * (distance to vacuum boundary))
        phi(i) = C * sinh(kappa_r * (a + b - abs_x));
    end
end

disp(phi)

% Shift x so it starts at 0 (Total width L = 2a + 2b)
x_shifted = x_plot + (a + b);

% ==========================================
% 6. Plotting
% ==========================================
figure('Color', 'w');
plot(x_shifted, phi, 'Color', [0.5 0 0.5], 'LineWidth', 2);
hold on; grid on;

% Interface lines
xline(b, '--', 'HandleVisibility', 'on', 'Label', 'Interface');
xline(2*a + b, '--', 'HandleVisibility', 'off');

title('Analytical Neutron Flux in a Reflected Slab Reactor');
xlabel('Position x (cm)');
ylabel('Normalized Flux \phi(x)');
legend('Neutron Flux', 'Core-Reflector Interface');