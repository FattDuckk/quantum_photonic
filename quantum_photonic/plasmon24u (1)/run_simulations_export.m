%% MATLAB Data Export Script for Python Analysis
% This script runs simulations and saves results in JSON format
% for easy import into Python notebook

%% Setup
clear
close all
clc

results = struct();

%% 1. Basic Example - Single Shell
fprintf('========================================\n');
fprintf('Running Basic Example (K shell)...\n');
fprintf('========================================\n\n');

% Parameters
a = 0.020;  % 20 nm radius
f = 0.3;    % volume fraction
nc = 1.5;   % core index
nb = 1.5;   % background index
metal = 'K';

% Get material data
[et,Lt] = feval(['eps',metal]);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0,1,1e3).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

% Setup geometry
N = 10;
m = [nc*ones(size(L)), sqrt(feval(['eps',metal],L)), nb*ones(size(L))];
mu = ones(size(m));
av = a*[(1-f)^(1/3), 1];
vol = 4*pi/3*a^3;

% Calculate spectrum
C = nan([length(L),3]);
for k = 1:length(L)
    L1 = L(k);
    coeffs = fmiecoeff(N, L1, av, m(k,:), mu(k,:));
    C1 = crosssection(coeffs, L1, av, m(k,:));
    C1 = C1/vol/(2*pi/L1);  % Normalize
    C(k,:) = C1;
end

% Find peaks
[pks_abs, locs_abs] = findpeaks(C(:,2), 'NPeaks', 2, 'SortStr', 'descend');
[pks_sca, locs_sca] = findpeaks(C(:,3), 'NPeaks', 2, 'SortStr', 'descend');

% Sort by energy (location)
[~, idx_abs] = sort(locs_abs);
[~, idx_sca] = sort(locs_sca);
locs_abs = locs_abs(idx_abs);
pks_abs = pks_abs(idx_abs);
locs_sca = locs_sca(idx_sca);
pks_sca = pks_sca(idx_sca);

% Store results
results.basic.metal = metal;
results.basic.outer_radius_nm = a * 1e3;
results.basic.volume_fraction = f;
results.basic.core_index = nc;
results.basic.background_index = nb;

if length(locs_abs) >= 2
    results.basic.peak1.energy_eV = eV(locs_abs(1));
    results.basic.peak1.wavelength_nm = 1240 / eV(locs_abs(1));
    results.basic.peak1.C_abs = pks_abs(1);
    results.basic.peak2.energy_eV = eV(locs_abs(2));
    results.basic.peak2.wavelength_nm = 1240 / eV(locs_abs(2));
    results.basic.peak2.C_abs = pks_abs(2);
end

if length(locs_sca) >= 2
    results.basic.peak1.C_sca = C(locs_abs(1), 3);
    results.basic.peak2.C_sca = C(locs_abs(2), 3);
end

fprintf('Peak 1 (bonding): %.3f eV (%.1f nm)\n', ...
    results.basic.peak1.energy_eV, results.basic.peak1.wavelength_nm);
fprintf('Peak 2 (antibonding): %.3f eV (%.1f nm)\n', ...
    results.basic.peak2.energy_eV, results.basic.peak2.wavelength_nm);
fprintf('Energy splitting: %.3f eV\n\n', ...
    results.basic.peak2.energy_eV - results.basic.peak1.energy_eV);

% Plot
figure(1);
clf;
plot(eV, C(:,2), 'b-', 'LineWidth', 2);
hold on;
plot(eV, C(:,3), 'r--', 'LineWidth', 2);
plot(eV(locs_abs), pks_abs, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('Energy (eV)');
ylabel('Cross-section (normalized)');
title(sprintf('%s shell: f=%.2f, a=%.0f nm', metal, f, a*1e3));
legend('Absorption', 'Scattering', 'Peaks');
grid on;

%% 2. Shell Thickness Study
fprintf('========================================\n');
fprintf('Running Shell Thickness Study (Au)...\n');
fprintf('========================================\n\n');

metal = 'Au';
a = 0.010;  % 10 nm
nc = 1;
nb = 1;
f_values = [0.1, 0.3, 0.5, 0.7, 0.9];

[et,Lt] = feval(['eps',metal]);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0,1,1e3).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

thickness_data = struct();

for idx = 1:length(f_values)
    f = f_values(idx);
    
    fprintf('f = %.2f... ', f);
    
    % Setup
    m = [nc*ones(size(L)), sqrt(feval(['eps',metal],L)), nb*ones(size(L))];
    mu = ones(size(m));
    av = a*[(1-f)^(1/3), 1];
    vol = 4*pi/3*a^3;
    
    % Calculate
    C = nan([length(L),3]);
    for k = 1:length(L)
        L1 = L(k);
        coeffs = fmiecoeff(N, L1, av, m(k,:), mu(k,:));
        C1 = crosssection(coeffs, L1, av, m(k,:));
        C1 = C1/vol/(2*pi/L1);
        C(k,:) = C1;
    end
    
    % Find peaks
    [~, loc_abs1] = max(C(:,2));
    [~, loc_sca2] = max(C(:,3));
    
    % Store
    thickness_data(idx).f = f;
    thickness_data(idx).shell_thickness_nm = a * 1e3 * (1 - (1-f)^(1/3));
    thickness_data(idx).peak1_energy_eV = eV(loc_abs1);
    thickness_data(idx).peak1_wavelength_nm = 1240 / eV(loc_abs1);
    thickness_data(idx).peak1_C_abs = C(loc_abs1, 2);
    thickness_data(idx).peak1_C_sca = C(loc_abs1, 3);
    thickness_data(idx).peak2_energy_eV = eV(loc_sca2);
    thickness_data(idx).peak2_wavelength_nm = 1240 / eV(loc_sca2);
    thickness_data(idx).peak2_C_abs = C(loc_sca2, 2);
    thickness_data(idx).peak2_C_sca = C(loc_sca2, 3);
    thickness_data(idx).energy_splitting_eV = abs(eV(loc_sca2) - eV(loc_abs1));
    
    fprintf('E1=%.3f eV, E2=%.3f eV, ΔE=%.3f eV\n', ...
        thickness_data(idx).peak1_energy_eV, ...
        thickness_data(idx).peak2_energy_eV, ...
        thickness_data(idx).energy_splitting_eV);
end

results.thickness_study = thickness_data;

%% 3. Metal Comparison
fprintf('\n========================================\n');
fprintf('Running Metal Comparison Study...\n');
fprintf('========================================\n\n');

metals = {'Au', 'Ag', 'Al', 'Cu'};
f = 0.5;
a = 0.010;
nc = 1;
nb = 1;

metal_data = struct();

for m_idx = 1:length(metals)
    metal = metals{m_idx};
    fprintf('%s... ', metal);
    
    [et,Lt] = feval(['eps',metal]);
    Lmin = min(Lt(real(et)<0));
    eV = 1.24*linspace(0,1,1e3).'/Lmin;
    eV = eV(eV>0);
    L = 1.24./eV;
    
    m = [nc*ones(size(L)), sqrt(feval(['eps',metal],L)), nb*ones(size(L))];
    mu = ones(size(m));
    av = a*[(1-f)^(1/3), 1];
    vol = 4*pi/3*a^3;
    
    C = nan([length(L),3]);
    for k = 1:length(L)
        L1 = L(k);
        coeffs = fmiecoeff(N, L1, av, m(k,:), mu(k,:));
        C1 = crosssection(coeffs, L1, av, m(k,:));
        C1 = C1/vol/(2*pi/L1);
        C(k,:) = C1;
    end
    
    [~, loc_abs1] = max(C(:,2));
    [~, loc_sca2] = max(C(:,3));
    
    metal_data(m_idx).metal = metal;
    metal_data(m_idx).peak1_energy_eV = eV(loc_abs1);
    metal_data(m_idx).peak1_wavelength_nm = 1240 / eV(loc_abs1);
    metal_data(m_idx).peak2_energy_eV = eV(loc_sca2);
    metal_data(m_idx).peak2_wavelength_nm = 1240 / eV(loc_sca2);
    metal_data(m_idx).energy_splitting_eV = abs(eV(loc_sca2) - eV(loc_abs1));
    
    fprintf('E1=%.3f eV, E2=%.3f eV\n', ...
        metal_data(m_idx).peak1_energy_eV, ...
        metal_data(m_idx).peak2_energy_eV);
end

results.metal_comparison = metal_data;

%% Save results as JSON
fprintf('\n========================================\n');
fprintf('Saving results to JSON...\n');
fprintf('========================================\n');

json_str = jsonencode(results);
fid = fopen('simulation_results.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Results saved to: simulation_results.json\n');
fprintf('Import this file in Python notebook for analysis!\n');

%% Print summary table
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n\n');

fprintf('THICKNESS STUDY (Au, a=10nm):\n');
fprintf('f    | d(nm) | E1(eV) | E2(eV) | ΔE(eV) |\n');
fprintf('-----|-------|--------|--------|--------|\n');
for idx = 1:length(results.thickness_study)
    d = results.thickness_study(idx);
    fprintf('%.2f | %5.2f | %6.3f | %6.3f | %6.3f |\n', ...
        d.f, d.shell_thickness_nm, d.peak1_energy_eV, ...
        d.peak2_energy_eV, d.energy_splitting_eV);
end

fprintf('\nMETAL COMPARISON (f=0.5, a=10nm):\n');
fprintf('Metal | E1(eV) | λ1(nm) | E2(eV) | λ2(nm) | ΔE(eV) |\n');
fprintf('------|--------|--------|--------|--------|--------|\n');
for idx = 1:length(results.metal_comparison)
    d = results.metal_comparison(idx);
    fprintf('%5s | %6.3f | %6.1f | %6.3f | %6.1f | %6.3f |\n', ...
        d.metal, d.peak1_energy_eV, d.peak1_wavelength_nm, ...
        d.peak2_energy_eV, d.peak2_wavelength_nm, ...
        d.energy_splitting_eV);
end

fprintf('\nDone! Open plasmon_analysis.ipynb to continue analysis.\n');
