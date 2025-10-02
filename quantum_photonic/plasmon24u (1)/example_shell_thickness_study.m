%% Shell Thickness Study
% This script runs the Mie simulation for different shell thicknesses
% to investigate the plasmon hybridization effect
%
% Based on example_mie_spectrum_3d_shell_eV.m

%% cleanup
clear
close all
clc

%% Fixed parameters (as recommended)
a = 0.010;        % Outer radius = 10 nm
nc = 1;           % Core index (air)
nb = 1;           % Background index (air)
metal = 'Au';     % Metal type

%% Shell thickness values to test
f_values = [0.1, 0.3, 0.5, 0.7, 0.9];
colors = lines(length(f_values));

%% Pre-calculate material data
[et,Lt] = feval(['eps',metal]);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0,1,1e3).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

%% Setup figure
figure(1)
clf
hold on
xlabel('Energy (eV)')
ylabel('C / V / k')
title(sprintf('%s shell: Varying shell thickness (a = %g nm)', metal, a*1e3))

%% Storage for peak positions
peak_data = zeros(length(f_values), 6); % f, abs_peak_eV, abs_peak_val, sca_peak_eV, sca_peak_val, separation

%% Calculate for each shell thickness
N = 10; % Number of orders

for idx = 1:length(f_values)
    f = f_values(idx);
    
    fprintf('Calculating for f = %.2f (shell thickness = %.2f nm)\n', ...
        f, a*1e3*(1-(1-f)^(1/3)));
    
    %% Setup geometry
    m = [nc*ones(size(L)), sqrt(feval(['eps',metal],L)), nb*ones(size(L))];
    mu = ones(size(m));
    av = a*[(1-f)^(1/3), 1];
    vol = 4*pi/3*a^3;
    
    %% Calculate spectrum
    C = nan([length(L),3]);
    for k = 1:length(L)
        L1 = L(k);
        coeffs = fmiecoeff(N, L1, av, m(k,:), mu(k,:));
        C1 = crosssection(coeffs, L1, av, m(k,:));
        C1 = C1/vol/(2*pi/L1); % Normalize
        C(k,:) = C1;
    end
    
    %% Plot results
    figure(1)
    plot(eV, C(:,2), 'Color', colors(idx,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('f=%.2f (abs)', f));
    plot(eV, C(:,3), '--', 'Color', colors(idx,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('f=%.2f (sca)', f));
    
    %% Find peaks
    [abs_peak_val, abs_idx] = max(C(:,2));
    [sca_peak_val, sca_idx] = max(C(:,3));
    abs_peak_eV = eV(abs_idx);
    sca_peak_eV = eV(sca_idx);
    
    % Store data
    peak_data(idx, :) = [f, abs_peak_eV, abs_peak_val, sca_peak_eV, sca_peak_val, ...
                         abs(sca_peak_eV - abs_peak_eV)];
    
    % Mark peaks
    plot(abs_peak_eV, abs_peak_val, 'o', 'Color', colors(idx,:), ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(idx,:), 'HandleVisibility', 'off');
    plot(sca_peak_eV, sca_peak_val, 's', 'Color', colors(idx,:), ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(idx,:), 'HandleVisibility', 'off');
end

legend('Location', 'best')
grid on

%% Summary figure
figure(2)
clf

subplot(2,2,1)
plot(f_values, peak_data(:,2), 'o-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Volume fraction f')
ylabel('Peak energy (eV)')
title('Absorption Peak Position')
grid on

subplot(2,2,2)
plot(f_values, peak_data(:,4), 's-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Volume fraction f')
ylabel('Peak energy (eV)')
title('Scattering Peak Position')
grid on

subplot(2,2,3)
plot(f_values, 1.24./peak_data(:,2), 'o-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(f_values, 1.24./peak_data(:,4), 's-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Volume fraction f')
ylabel('Peak wavelength (μm)')
legend('Absorption', 'Scattering')
title('Peak Wavelengths')
grid on

subplot(2,2,4)
plot(f_values, peak_data(:,6), 'd-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Volume fraction f')
ylabel('Peak separation (eV)')
title('Energy Splitting Between Modes')
grid on

%% Print summary table
fprintf('\n=== SUMMARY TABLE ===\n');
fprintf('f    | Abs Peak (eV) | Abs Peak (nm) | Sca Peak (eV) | Sca Peak (nm) | Splitting (eV)\n');
fprintf('-----|---------------|---------------|---------------|---------------|---------------\n');
for idx = 1:length(f_values)
    fprintf('%.2f | %13.3f | %13.1f | %13.3f | %13.1f | %14.3f\n', ...
        peak_data(idx,1), peak_data(idx,2), 1240/peak_data(idx,2), ...
        peak_data(idx,4), 1240/peak_data(idx,4), peak_data(idx,6));
end
fprintf('\n');

%% Physical interpretation guide
fprintf('=== INTERPRETATION GUIDE ===\n');
fprintf('• Thin shell (low f): Strong coupling → Large splitting\n');
fprintf('• Thick shell (high f): Weak coupling → Small splitting\n');
fprintf('• As f→1 (solid sphere): Modes converge to dipole resonance\n');
fprintf('• Peak shift: Related to plasmon hybridization model\n');
fprintf('\n');
