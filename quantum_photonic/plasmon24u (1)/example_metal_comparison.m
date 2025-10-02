%% Metal Comparison Study
% This script compares different plasmonic metals
% to see how material properties affect plasmon resonances
%
% Based on example_mie_spectrum_3d_shell_eV.m

%% cleanup
clear
close all
clc

%% Fixed parameters
a = 0.010;        % Outer radius = 10 nm
f = 0.5;          % Shell volume fraction (medium shell)
nc = 1;           % Core index (air)
nb = 1;           % Background index (air)

%% Metals to compare
metals = {'Ag', 'Au', 'Al', 'Cu'};
metal_names = {'Silver (Ag)', 'Gold (Au)', 'Aluminum (Al)', 'Copper (Cu)'};
colors = [0.75 0.75 0.75;  % Silver (gray)
          1.0 0.84 0.0;     % Gold (gold)
          0.5 0.5 1.0;      % Aluminum (light blue)
          0.72 0.45 0.20];  % Copper (copper)

%% Setup common energy range
% Need to find overlapping range for all metals
eV_min = 1.0;  % Start from 1 eV
eV_max = 5.0;  % Go up to 5 eV
eV = linspace(eV_min, eV_max, 500)';
L = 1.24./eV;

%% Setup figure
figure(1)
clf
set(gcf, 'Position', [100 100 1400 900])

%% Storage for comparison
peak_data = zeros(length(metals), 7); % metal_idx, abs_peak_eV, abs_val, abs_width, sca_peak_eV, sca_val, sca_width

%% Calculate for each metal
N = 10; % Number of orders

for idx = 1:length(metals)
    metal = metals{idx};
    
    fprintf('Calculating for %s...\n', metal_names{idx});
    
    %% Check material range
    [et, Lt] = feval(['eps', metal]);
    if max(Lt) < max(L) || min(Lt) > min(L)
        warning('%s data may not cover full wavelength range', metal);
    end
    
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
    
    %% Plot extinction
    subplot(2,2,1)
    hold on
    plot(eV, C(:,1), 'Color', colors(idx,:), 'LineWidth', 2, ...
        'DisplayName', metal_names{idx});
    xlabel('Energy (eV)')
    ylabel('Extinction (C/V/k)')
    title('Extinction Spectra')
    grid on
    
    %% Plot absorption
    subplot(2,2,2)
    hold on
    plot(eV, C(:,2), 'Color', colors(idx,:), 'LineWidth', 2, ...
        'DisplayName', metal_names{idx});
    xlabel('Energy (eV)')
    ylabel('Absorption (C/V/k)')
    title('Absorption Spectra')
    grid on
    
    %% Plot scattering
    subplot(2,2,3)
    hold on
    plot(eV, C(:,3), 'Color', colors(idx,:), 'LineWidth', 2, ...
        'DisplayName', metal_names{idx});
    xlabel('Energy (eV)')
    ylabel('Scattering (C/V/k)')
    title('Scattering Spectra')
    grid on
    
    %% Find peaks and widths
    [abs_peak_val, abs_idx] = max(C(:,2));
    [sca_peak_val, sca_idx] = max(C(:,3));
    abs_peak_eV = eV(abs_idx);
    sca_peak_eV = eV(sca_idx);
    
    % Estimate FWHM (Full Width at Half Maximum)
    abs_half = abs_peak_val / 2;
    abs_half_idx = find(C(:,2) > abs_half);
    if length(abs_half_idx) > 1
        abs_width = eV(abs_half_idx(end)) - eV(abs_half_idx(1));
    else
        abs_width = NaN;
    end
    
    sca_half = sca_peak_val / 2;
    sca_half_idx = find(C(:,3) > sca_half);
    if length(sca_half_idx) > 1
        sca_width = eV(sca_half_idx(end)) - eV(sca_half_idx(1));
    else
        sca_width = NaN;
    end
    
    % Store data
    peak_data(idx, :) = [idx, abs_peak_eV, abs_peak_val, abs_width, ...
                         sca_peak_eV, sca_peak_val, sca_width];
    
    % Mark peaks
    subplot(2,2,2)
    plot(abs_peak_eV, abs_peak_val, 'o', 'Color', colors(idx,:), ...
        'MarkerSize', 10, 'MarkerFaceColor', colors(idx,:), 'HandleVisibility', 'off');
    
    subplot(2,2,3)
    plot(sca_peak_eV, sca_peak_val, 's', 'Color', colors(idx,:), ...
        'MarkerSize', 10, 'MarkerFaceColor', colors(idx,:), 'HandleVisibility', 'off');
end

% Add legends
subplot(2,2,1), legend('Location', 'best')
subplot(2,2,2), legend('Location', 'best')
subplot(2,2,3), legend('Location', 'best')

%% Comparison bar chart
subplot(2,2,4)
x_pos = 1:length(metals);
bar_width = 0.35;
b1 = bar(x_pos - bar_width/2, peak_data(:,2), bar_width, 'FaceColor', [0.2 0.6 0.8]);
hold on
b2 = bar(x_pos + bar_width/2, peak_data(:,5), bar_width, 'FaceColor', [0.8 0.4 0.2]);
set(gca, 'XTick', x_pos, 'XTickLabel', metals)
ylabel('Peak Energy (eV)')
title('Peak Positions Comparison')
legend('Absorption', 'Scattering', 'Location', 'best')
grid on

%% Print summary table
fprintf('\n=== METAL COMPARISON SUMMARY ===\n');
fprintf('Metal | Abs Peak | Abs λ | Abs Width | Sca Peak | Sca λ | Sca Width | Quality\n');
fprintf('------|----------|-------|-----------|----------|-------|-----------|--------\n');
for idx = 1:length(metals)
    Q_factor = peak_data(idx,2) / peak_data(idx,4); % Peak position / width
    fprintf('%-5s | %7.3f  | %5.0f | %8.3f  | %7.3f  | %5.0f | %8.3f  | %7.2f\n', ...
        metals{idx}, peak_data(idx,2), 1240/peak_data(idx,2), peak_data(idx,4), ...
        peak_data(idx,5), 1240/peak_data(idx,5), peak_data(idx,7), Q_factor);
end
fprintf('\n');

%% Create wavelength version
figure(2)
clf
subplot(1,2,1)
hold on
for idx = 1:length(metals)
    metal = metals{idx};
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
    
    plot(L*1e3, C(:,1), 'Color', colors(idx,:), 'LineWidth', 2, ...
        'DisplayName', metal_names{idx});
end
xlabel('Wavelength (nm)')
ylabel('Extinction (C/V/k)')
title('Extinction Spectra (Wavelength Scale)')
legend('Location', 'best')
grid on
set(gca, 'XDir', 'reverse') % Reverse to match energy scale

subplot(1,2,2)
% Show absorption/scattering ratio
for idx = 1:length(metals)
    metal = metals{idx};
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
    
    ratio = C(:,2) ./ C(:,3); % absorption / scattering
    plot(eV, ratio, 'Color', colors(idx,:), 'LineWidth', 2, ...
        'DisplayName', metal_names{idx});
    hold on
end
xlabel('Energy (eV)')
ylabel('Absorption / Scattering Ratio')
title('Loss Character')
legend('Location', 'best')
grid on
yline(1, 'k--', 'Equal', 'HandleVisibility', 'off')

%% Physical interpretation guide
fprintf('=== INTERPRETATION GUIDE ===\n');
fprintf('• Silver (Ag): Lowest loss, sharpest peaks (best Q-factor)\n');
fprintf('• Gold (Au): Good for red/NIR, more lossy than Ag\n');
fprintf('• Aluminum (Al): UV plasmons, higher damping\n');
fprintf('• Copper (Cu): Similar to Au, interband transitions add loss\n');
fprintf('\n');
fprintf('Peak Width: Narrower = less damping = longer plasmon lifetime\n');
fprintf('Peak Position: Depends on plasma frequency and band structure\n');
fprintf('Quality Factor Q = E_peak / FWHM (higher is better)\n');
fprintf('\n');
