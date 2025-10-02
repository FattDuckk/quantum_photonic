%% Metal comparison: Different metals vs Plasmon resonance energy

clear; close all; clc;

%% Fixed parameters
a = 0.010;      % sphere radius (microns)
f = 0.5;        % fixed volume fraction (50% shell)
nc = 1.5;       % core index 
nb = 1;       % background index
N = 10;         % Mie expansion orders

%% Metals to compare
metals = {'Au', 'Ag', 'Cu', 'Al', 'K'};
metal_names = {'Gold', 'Silver', 'Copper', 'Aluminum', 'Potassium'};

%% Calculate resonance for each metal
resonance_energies = zeros(size(metals));

for i = 1:length(metals)
    metal = metals{i};
    fprintf('Calculating %s...\n', metal_names{i});
    
    % Energy setup (specific to each metal)
    [et,Lt] = feval(['eps',metal]);
    Lmin = min(Lt(real(et)<0));
    eV = 1.24*linspace(0,1,1e3).'/Lmin;
    eV = eV(eV>0);
    L = 1.24./eV;
    
    % Setup
    fm = ['eps',metal];
    m = [nc*ones(size(L)), sqrt(feval(fm,L)), nb*ones(size(L))];
    av = a*[(1-f)^(1/3), 1];
    vol = 4*pi/3*a^3;
    
    % Calculate spectrum
    C = nan([length(L),3]);
    Ls = nan(size(L));
    
    for k = 1:length(L)
        L1 = L(k);
        coeffs = fmiecoeff(N, L1, av, m(k,:));
        C1 = crosssection(coeffs, L1, av, m(k,:));
        C1 = C1/vol/(2*pi/L1);
        C(k,:) = C1;
        Ls(k) = 1.24/L(k);
    end
    
    % Find extinction peak
    [~, peak_idx] = max(C(:,1));
    resonance_energies(i) = Ls(peak_idx);
    
    fprintf('  Peak at %.2f eV\n', resonance_energies(i));
end

%% Plot
figure;
bar(resonance_energies);
set(gca, 'XTickLabel', metal_names);
xlabel('Metal Type');
ylabel('Plasmon Resonance Energy (eV)');
title('LSP Resonance Energy by Metal Type');
grid on;

% Add values on bars
for i = 1:length(resonance_energies)
    text(i, resonance_energies(i) + 0.05, sprintf('%.2f eV', resonance_energies(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

%% Summary
fprintf('\n=== METAL COMPARISON RESULTS ===\n');
fprintf('Metal      | Resonance Energy (eV)\n');
fprintf('-----------|---------------------\n');
for i = 1:length(metals)
    fprintf('%-10s | %.2f\n', metal_names{i}, resonance_energies(i));
end

% Find range
energy_range = max(resonance_energies) - min(resonance_energies);
fprintf('\nEnergy range across metals: %.2f eV\n', energy_range);

[~, highest_idx] = max(resonance_energies);
[~, lowest_idx] = min(resonance_energies);
fprintf('Highest energy: %s (%.2f eV)\n', metal_names{highest_idx}, resonance_energies(highest_idx));
fprintf('Lowest energy: %s (%.2f eV)\n', metal_names{lowest_idx}, resonance_energies(lowest_idx));