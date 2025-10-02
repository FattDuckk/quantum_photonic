%% Simple plot: Shell thickness vs Scattering peak energy

clear; close all; clc;

%% Parameters
a = 0.020; nc = 1.5; nb = 1.5; metal = 'K'; N = 10;
f_values = 0.1:0.1:0.9;

%% Energy setup
[et,Lt] = feval(['eps',metal]);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0,1,1e3).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

%% Calculate peaks
scatter_peaks = zeros(size(f_values));
shell_thickness_nm = zeros(size(f_values));

for i = 1:length(f_values)
    f = f_values(i);
    
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
    
    % Find scattering peak
    [~, peak_idx] = max(C(:,3));
    scatter_peaks(i) = Ls(peak_idx);
    
    % Shell thickness in nm
    core_radius = a * (1-f)^(1/3);
    shell_thickness_nm(i) = (a - core_radius) * 1e3;
end

%% Plot
figure;
plot(shell_thickness_nm, scatter_peaks, 'o-');
xlabel('Shell Thickness (nm)');
ylabel('Plasmon Resonance Energy (eV)');
title('Shell Thickness vs Plasmon Resonance Energy ');
grid on;