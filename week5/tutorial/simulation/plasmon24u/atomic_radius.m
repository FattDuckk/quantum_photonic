%% Particle size study: Total radius vs Plasmon resonance energy (fixed shell ratio)

clear; close all; clc;

%% Fixed parameters
f = 0.5;        % fixed volume fraction (50% shell)
nc = 1.5;       % core index 
nb = 1.5;       % background index
metal = 'Ag';   % silver
N = 10;         % Mie expansion orders

%% Variable parameter: total particle radius
radius_nm = 10:5:50;  % from 10 nm to 50 nm radius in 5 nm steps
radius_um = radius_nm * 1e-3;  % convert to microns

%% Energy setup (same for all sizes)
[et,Lt] = feval(['eps',metal]);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0,1,1e3).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

%% Calculate resonance for each particle size
resonance_energies = zeros(size(radius_um));
shell_thickness_nm = zeros(size(radius_um));

for i = 1:length(radius_um)
    a = radius_um(i);
    fprintf('Calculating radius = %.0f nm...\n', a*1e3);
    
    % Setup
    fm = ['eps',metal];
    m = [nc*ones(size(L)), sqrt(feval(fm,L)), nb*ones(size(L))];
    av = a*[(1-f)^(1/3), 1];
    vol = 4*pi/3*a^3;
    
    % Calculate actual shell thickness for this size
    core_radius = a * (1-f)^(1/3);
    shell_thickness_nm(i) = (a - core_radius) * 1e3;
    
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
    
    fprintf('  Peak at %.2f eV (shell = %.1f nm)\n', resonance_energies(i), shell_thickness_nm(i));
end

%% Plot
figure;
plot(radius_nm, resonance_energies, 'o-');
xlabel('Total Particle Radius (nm)');
ylabel('Plasmon Resonance Energy (eV)');
title('Particle Size vs Resonance Energy (f = 0.5)');
grid on;

%% Summary
fprintf('\n=== PARTICLE SIZE RESULTS ===\n');
fprintf('Radius (nm) | Shell (nm) | Resonance (eV)\n');
fprintf('------------|------------|---------------\n');
for i = 1:length(radius_nm)
    fprintf('%8.0f    | %8.1f   | %8.2f\n', radius_nm(i), shell_thickness_nm(i), resonance_energies(i));
end

energy_range = max(resonance_energies) - min(resonance_energies);
fprintf('\nResonance energy range: %.2f eV\n', energy_range);
fprintf('Size range: %d to %d nm radius\n', min(radius_nm), max(radius_nm));