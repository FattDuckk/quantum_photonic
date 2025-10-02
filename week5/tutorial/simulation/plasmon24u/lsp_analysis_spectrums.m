%%
compare_materials();
%%
function [C, Ls, peak_data] = calculate_lsp_spectrum(metal, a, f, nc, nb)
% Calculate LSP spectrum for given parameters
% Returns: C (cross-sections), Ls (energies in eV), peak_data (resonance info)

% Get material data
[et,Lt] = feval(['eps',metal]);
Lmin = min(Lt(real(et)<0));
eV = 1.24*linspace(0.1,1,800).'/Lmin;
eV = eV(eV>0);
L = 1.24./eV;

% Setup particle geometry
N = 10;
m = [nc*ones(size(L)), sqrt(feval(['eps',metal],L)), nb*ones(size(L))];
av = a*[(1-f)^(1/3), 1];
vol = 4*pi/3*a^3;
mu = ones(size(m));

% Calculate spectrum
C = nan([length(L),3]);
Ls = nan(size(L));

for k = 1:length(L)
    L1 = L(k);
    coeffs = fmiecoeff(N, L1, av, m(k,:), mu(k,:));
    C1 = crosssection(coeffs, L1, av, m(k,:));
    C1 = C1/vol/(2*pi/L1); % normalize
    C(k,:) = C1;
    Ls(k) = 1.24/L(k); % energy in eV
end

% Find peak information
[peak_ext, idx_ext] = max(C(:,1));
[peak_abs, idx_abs] = max(C(:,2));
[peak_sca, idx_sca] = max(C(:,3));

peak_data = struct();
peak_data.ext_energy = Ls(idx_ext);
peak_data.ext_wavelength = 1240/Ls(idx_ext); % nm
peak_data.ext_intensity = peak_ext;
peak_data.abs_energy = Ls(idx_abs);
peak_data.abs_wavelength = 1240/Ls(idx_abs);
peak_data.sca_energy = Ls(idx_sca);
peak_data.sca_wavelength = 1240/Ls(idx_sca);

end

function compare_materials()
% Week 5/6: Systematic material comparison with expanded range

metals = {'K', 'Au', 'Ag', 'Cu', 'Al'};
colors = {'b', 'r', 'g', 'm', 'c'};

% Fixed parameters
a = 0.020; % 20nm radius
f = 0.9;   % 90% metal fraction
nc = 1.5;  % core index
nb = 1.5;  % background index

figure('Position', [100, 100, 1200, 800]);

% Store results for analysis
results = [];

for i = 1:length(metals)
    [C, Ls, peak_data] = calculate_lsp_spectrum(metals{i}, a, f, nc, nb);
    
    % Plot extinction spectra
    subplot(2,3,1)
    plot(Ls, C(:,1), colors{i}, 'LineWidth', 2);
    hold on;
    
    % Plot absorption vs scattering
    subplot(2,3,2)
    plot(Ls, C(:,2), colors{i}, 'LineWidth', 2); % absorption
    hold on;
    plot(Ls, C(:,3), [colors{i} '--'], 'LineWidth', 1.5); % scattering
    
    % Store data
    results = [results; {metals{i}, peak_data.ext_energy, peak_data.ext_wavelength, ...
                        peak_data.abs_energy, peak_data.sca_energy}];
    
    fprintf('%s: Peak at %.2f eV (%.0f nm)\n', metals{i}, peak_data.ext_energy, peak_data.ext_wavelength);
end

% Format plots
subplot(2,3,1)
xlabel('Energy (eV)');
ylabel('Extinction Cross-section');
title('LSP Extinction Spectra - Material Comparison');
legend(metals, 'Location', 'best');
grid on;

subplot(2,3,2)
xlabel('Energy (eV)');
ylabel('Cross-section');
title('Absorption (solid) vs Scattering (dashed)');
grid on;

% Peak wavelength comparison
subplot(2,3,3)
wavelengths = cellfun(@(x) x{3}, results);
bar(wavelengths);
set(gca, 'XTickLabel', metals);
xlabel('Metal');
ylabel('Peak Wavelength (nm)');
title('LSP Peak Wavelengths');
grid on;

% Energy comparison
subplot(2,3,4)
energies = cellfun(@(x) x{2}, results);
bar(energies);
set(gca, 'XTickLabel', metals);
xlabel('Metal');
ylabel('Peak Energy (eV)');
title('LSP Peak Energies');
grid on;

% Theoretical comparison
subplot(2,3,5)
% Plasma frequencies from material files
plasma_freqs = [3.95, 8.4, 9.07, 8.95, 11.49]; % K, Au, Ag, Cu, Al in eV
theory_energies = plasma_freqs / sqrt(3); % simple prediction
scatter(theory_energies, energies, 100, 'filled');
hold on;
plot([0 8], [0 8], 'k--', 'LineWidth', 1);
xlabel('Theoretical Energy (eV)');
ylabel('Experimental Energy (eV)');
title('Theory vs Experiment');
for i = 1:length(metals)
    text(theory_energies(i), energies(i), ['  ' metals{i}]);
end
grid on;

% Summary table
subplot(2,3,6)
axis off;
table_str = sprintf('Metal   Energy(eV)  Wavelength(nm)\n');
table_str = [table_str '--------------------------------\n'];
for i = 1:length(results)
    table_str = [table_str sprintf('%-6s  %8.2f  %12.0f\n', ...
                 results{i}{1}, results{i}{2}, results{i}{3})];
end
text(0.1, 0.8, table_str, 'FontSize', 10);
title('Summary Results');

end

function study_size_range()
% Week 6: Size effects with expanded range

metal = 'Au'; % Focus on one well-known material
sizes = [0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050]; % 5-50nm range
f = 0.9;
nc = 1.5;
nb = 1.5;

figure('Position', [100, 100, 1000, 700]);
colors = jet(length(sizes));

peak_wavelengths = [];
peak_energies = [];
abs_sca_ratios = [];

for i = 1:length(sizes)
    [C, Ls, peak_data] = calculate_lsp_spectrum(metal, sizes(i), f, nc, nb);
    
    % Plot spectra
    subplot(2,2,1)
    plot(Ls, C(:,1), 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
    
    % Store peak data
    peak_wavelengths(i) = peak_data.ext_wavelength;
    peak_energies(i) = peak_data.ext_energy;
    abs_sca_ratios(i) = peak_data.ext_intensity / max(C(:,3)); % ratio at peak
end

% Format main plot
subplot(2,2,1)
xlabel('Energy (eV)');
ylabel('Extinction Cross-section');
title([metal ' LSP: Size Effects']);
legend(arrayfun(@(x) sprintf('%.0fnm', x*1000), sizes, 'UniformOutput', false), ...
       'Location', 'eastoutside');
grid on;

% Peak wavelength vs size
subplot(2,2,2)
plot(sizes*1000, peak_wavelengths, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Particle Radius (nm)');
ylabel('Peak Wavelength (nm)');
title('LSP Position Independence');
grid on;

% Size parameter analysis
subplot(2,2,3)
size_params = 2*pi*sizes./(peak_wavelengths*1e-3); % ka at resonance
plot(size_params, peak_energies, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Size Parameter (ka)');
ylabel('Peak Energy (eV)');
title('Quasistatic Validity');
xline(0.2, 'k--', 'Quasistatic Limit');
grid on;

% Cross-section scaling
subplot(2,2,4)
volumes = (4/3)*pi*(sizes*1e-9).^3; % in m³
geometric_areas = pi*(sizes*1e-9).^2; % in m²
loglog(sizes*1000, volumes*1e27, 'b-', 'LineWidth', 2); % Volume in nm³
hold on;
loglog(sizes*1000, geometric_areas*1e18, 'r--', 'LineWidth', 2); % Area in nm²
xlabel('Radius (nm)');
ylabel('Cross-section (nm^2 or nm^3)');
title('Scaling Laws');
legend('Volume ∝ r³', 'Area ∝ r²', 'Location', 'best');
grid on;

% Print size analysis
fprintf('\n=== SIZE EFFECT ANALYSIS ===\n');
fprintf('Metal: %s\n', metal);
fprintf('Radius(nm)  Wavelength(nm)  Energy(eV)  Size_param\n');
for i = 1:length(sizes)
    fprintf('%8.0f    %10.0f    %8.2f    %8.3f\n', ...
            sizes(i)*1000, peak_wavelengths(i), peak_energies(i), size_params(i));
end
fprintf('Wavelength std dev: %.1f nm (%.1f%% variation)\n', ...
        std(peak_wavelengths), std(peak_wavelengths)/mean(peak_wavelengths)*100);

end

function study_environment_effects()
% Week 6: Environment/sensing effects

metal = 'Au';
a = 0.020;
f = 0.9;
nc = 1.5;
backgrounds = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2]; % air to high-index

figure('Position', [100, 100, 1000, 700]);
colors = jet(length(backgrounds));

peak_shifts = [];
sensitivity = [];

for i = 1:length(backgrounds)
    [C, Ls, peak_data] = calculate_lsp_spectrum(metal, a, f, nc, backgrounds(i));
    
    % Plot spectra
    subplot(2,2,1)
    plot(Ls, C(:,1), 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
    
    peak_shifts(i) = peak_data.ext_wavelength;
end

% Format main plot
subplot(2,2,1)
xlabel('Energy (eV)');
ylabel('Extinction Cross-section');
title([metal ' LSP: Environment Effects']);
legend(arrayfun(@(x) sprintf('n=%.1f', x), backgrounds, 'UniformOutput', false), ...
       'Location', 'eastoutside');
grid on;

% Wavelength shift vs refractive index
subplot(2,2,2)
plot(backgrounds, peak_shifts, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Background Refractive Index');
ylabel('Peak Wavelength (nm)');
title('LSP Wavelength Shift');
grid on;

% Sensitivity analysis
subplot(2,2,3)
if length(backgrounds) > 1
    sensitivity = diff(peak_shifts)./diff(backgrounds);
    plot(backgrounds(2:end), sensitivity, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Background Index');
    ylabel('Sensitivity (nm/RIU)');
    title('Sensing Sensitivity');
    grid on;
end

% Theoretical comparison
subplot(2,2,4)
% Theory: omega ∝ 1/sqrt(1+2*n²)
theory_shifts = 1240 ./ (peak_shifts(1) * sqrt(1 + 2*backgrounds.^2) / sqrt(1 + 2*backgrounds(1)^2));
plot(backgrounds, peak_shifts, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(backgrounds, theory_shifts, 'b--', 'LineWidth', 2);
xlabel('Background Index');
ylabel('Peak Wavelength (nm)');
title('Theory vs Experiment');
legend('Experiment', 'Theory', 'Location', 'best');
grid on;

% Print environment analysis
fprintf('\n=== ENVIRONMENT EFFECT ANALYSIS ===\n');
fprintf('Metal: %s, Radius: %.0f nm\n', metal, a*1000);
fprintf('Background_n  Wavelength(nm)  Shift(nm)\n');
reference_wavelength = peak_shifts(1);
for i = 1:length(backgrounds)
    shift = peak_shifts(i) - reference_wavelength;
    fprintf('%10.1f    %12.0f    %8.1f\n', backgrounds(i), peak_shifts(i), shift);
end

end