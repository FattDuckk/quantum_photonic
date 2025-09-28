metals = {'K', 'Au', 'Ag', 'Cu'};
compare_metals(metals, 0.020, 0.9, 1.5, 1.5);
%%
sizes = [0.005, 0.010, 0.020, 0.030, 0.050]; % 5-50nm radii
study_size_effects('Au', sizes, 0.9, 1.5, 1.5);
%%
backgrounds = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0];
study_environment_effects('Au', 0.020, 0.9, 1.5, backgrounds);
%%



function [C, Ls, peak_wavelengths] = calculate_lsp_spectrum(metal, a, f, nc, nb)
% Calculate LSP spectrum for given parameters
% Inputs:
%   metal - string: 'Au', 'Ag', 'K', 'Cu', 'Al'
%   a - sphere radius in microns
%   f - metal fraction (0-1)
%   nc - core refractive index
%   nb - background refractive index
% Outputs:
%   C - cross-sections [ext, abs, sca]
%   Ls - energy scale in eV
%   peak_wavelengths - resonance wavelengths

% Get material data
[et,Lt]=feval(['eps',metal]);
Lmin = min(Lt);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0,1,1e3).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

% Setup particle
N=10;
fm=['eps',metal];
m=[nc*ones(size(L)), sqrt(feval(fm,L)), nb*ones(size(L))];
av = a*[(1-f)^(1/3), 1];
vol=4*pi/3*a^3;
mu=ones(size(m));

% Calculate spectrum
C=nan([length(L),3]);
Ls=nan(size(L));
fx=@(x) 1.24./x; % convert to eV

for k=1:length(L)
    L1=L(k);
    coeffs=fmiecoeff(N, L1, av, m(k,:), mu(k,:));
    C1=crosssection(coeffs, L1, av, m(k,:));
    C1=C1/vol/(2*pi/L1); % normalize by volume and k
    C(k,:)=C1;
    Ls(k)=fx(L(k));
end

% Find peak wavelengths
[~,ind]=max(C);
peak_wavelengths = 1.24./Ls(ind); % convert back to wavelengths in microns

end

function compare_metals(metals, a, f, nc, nb)
% Compare different metals on same plot
% metals - cell array: {'Au', 'Ag', 'K'}

figure;
colors = ['b', 'r', 'g', 'm', 'c'];
peak_data = [];

for i = 1:length(metals)
    [C, Ls, peaks] = calculate_lsp_spectrum(metals{i}, a, f, nc, nb);
    
    subplot(2,1,1)
    plot(Ls, C(:,1), 'Color', colors(i), 'LineWidth', 2); % extinction
    hold on;
    
    subplot(2,1,2)
    plot(Ls, C(:,2), '--', 'Color', colors(i), 'LineWidth', 2); % absorption
    plot(Ls, C(:,3), ':', 'Color', colors(i), 'LineWidth', 2); % scattering
    hold on;
    
    peak_data = [peak_data; {metals{i}, peaks(1)*1000}]; % wavelength in nm
end

% Format plots
subplot(2,1,1)
xlabel('Energy (eV)');
ylabel('Extinction C/V/k');
title('LSP Extinction Spectra');
legend(metals, 'Location', 'best');
grid on;

subplot(2,1,2)
xlabel('Energy (eV)');
ylabel('Cross-section C/V/k');
title('Absorption (solid) vs Scattering (dotted)');
legend([strcat(metals, '-abs'); strcat(metals, '-sca')], 'Location', 'best');
grid on;

% Print peak summary
fprintf('\nPeak Wavelengths:\n');
fprintf('Metal\tPeak (nm)\n');
for i = 1:size(peak_data,1)
    fprintf('%s\t%.1f\n', peak_data{i,1}, peak_data{i,2});
end

end

function study_size_effects(metal, sizes, f, nc, nb)
% Study how particle size affects LSP resonance
% sizes - array of radii in microns

figure;
colors = jet(length(sizes));

for i = 1:length(sizes)
    [C, Ls, peaks] = calculate_lsp_spectrum(metal, sizes(i), f, nc, nb);
    
    plot(Ls, C(:,1), 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
end

xlabel('Energy (eV)');
ylabel('Extinction C/V/k');
title(sprintf('%s LSP: Size Effects', metal));
legend(arrayfun(@(x) sprintf('r=%.0fnm', x*1000), sizes, 'UniformOutput', false), 'Location', 'best');
grid on;

end

function study_environment_effects(metal, a, f, nc, backgrounds)
% Study how background medium affects LSP
% backgrounds - array of background indices

figure;
colors = jet(length(backgrounds));
shifts = [];

for i = 1:length(backgrounds)
    [C, Ls, peaks] = calculate_lsp_spectrum(metal, a, f, nc, backgrounds(i));
    
    plot(Ls, C(:,1), 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
    
    shifts = [shifts, peaks(1)*1000]; % peak wavelength in nm
end

xlabel('Energy (eV)');
ylabel('Extinction C/V/k');
title(sprintf('%s LSP: Environment Effects', metal));
legend(arrayfun(@(x) sprintf('n=%.2f', x), backgrounds, 'UniformOutput', false), 'Location', 'best');
grid on;

% Plot wavelength shift
figure;
plot(backgrounds, shifts, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Background Refractive Index');
ylabel('Peak Wavelength (nm)');
title('LSP Wavelength vs Environment');
grid on;

end