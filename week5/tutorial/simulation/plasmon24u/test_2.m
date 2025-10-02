%% LSP Core-Shell Resonance Investigation
% Systematic study of plasmonic resonance tuning in metal-dielectric core-shell nanoparticles
% Investigates: shell thickness, particle size, metal type, background medium effects

%% Cleanup and Setup
clear
close all
clc

%% Investigation Parameters
% Core parameters
metals = {'Ag', 'Au', 'Cu', 'Al', 'K'}; % All available metals
core_fractions = 0.1:0.1:0.9; % Shell thickness variation (f = core_volume/total_volume)
particle_radii = [10, 20, 30, 40, 50] * 1e-3; % Total radius in microns
background_indices = [1.0, 1.33, 1.5]; % Air, water, glass
background_names = {'Air', 'Water', 'Glass'};

% Fixed parameters
nc = 1.5; % Core index (silica)
N = 15; % Number of Mie expansion orders
wavelength_range = [0.3, 1.2]; % Wavelength range in microns
n_wavelengths = 200;

% Create wavelength array
L = linspace(wavelength_range(1), wavelength_range(2), n_wavelengths);

%% Initialize Results Storage
results = struct();

% Pre-allocate storage arrays
n_metals = length(metals);
n_fractions = length(core_fractions);
n_radii = length(particle_radii);
n_backgrounds = length(background_indices);

% Storage for peak analysis
peak_wavelengths = nan(n_metals, n_fractions, n_radii, n_backgrounds, 2); % 2 peaks max
peak_intensities = nan(n_metals, n_fractions, n_radii, n_backgrounds, 2);
peak_quality_factors = nan(n_metals, n_fractions, n_radii, n_backgrounds, 2);
resonance_splitting = nan(n_metals, n_fractions, n_radii, n_backgrounds);

% Storage for cross-sections
extinction_spectra = nan(n_metals, n_fractions, n_radii, n_backgrounds, n_wavelengths);
absorption_spectra = nan(n_metals, n_fractions, n_radii, n_backgrounds, n_wavelengths);
scattering_spectra = nan(n_metals, n_fractions, n_radii, n_backgrounds, n_wavelengths);

%% Main Investigation Loop
fprintf('Starting LSP Core-Shell Investigation...\n');
fprintf('Progress: [');

total_calculations = n_metals * n_fractions * n_radii * n_backgrounds;
calc_count = 0;

for m_idx = 1:n_metals
    metal = metals{m_idx};
    fprintf('\nProcessing %s...', metal);
    
    for f_idx = 1:n_fractions
        f = core_fractions(f_idx);
        
        for r_idx = 1:n_radii
            a = particle_radii(r_idx);
            
            for nb_idx = 1:n_backgrounds
                nb = background_indices(nb_idx);
                
                calc_count = calc_count + 1;
                if mod(calc_count, round(total_calculations/50)) == 0
                    fprintf('=');
                end
                
                try
                    % Calculate shell geometry
                    av = a * [(1-f)^(1/3), 1]; % [core_radius, total_radius]
                    
                    % Calculate spectra
                    C_ext = zeros(size(L));
                    C_abs = zeros(size(L));
                    C_sca = zeros(size(L));
                    
                    for w_idx = 1:length(L)
                        L1 = L(w_idx);
                        
                        % Get metal permittivity
                        [et, ~] = feval(['eps', metal], L1);
                        m_array = [nc, sqrt(et), nb]; % [core, shell, background]
                        mu_array = ones(size(m_array));
                        
                        % Calculate Mie coefficients
                        coeffs = fmiecoeff(N, L1, av, m_array, mu_array);
                        
                        % Calculate cross-sections
                        C_result = crosssection(coeffs, L1, av, m_array);
                        C_ext(w_idx) = C_result(1);
                        C_abs(w_idx) = C_result(2);
                        C_sca(w_idx) = C_result(3);
                    end
                    
                    % Store spectra
                    extinction_spectra(m_idx, f_idx, r_idx, nb_idx, :) = C_ext;
                    absorption_spectra(m_idx, f_idx, r_idx, nb_idx, :) = C_abs;
                    scattering_spectra(m_idx, f_idx, r_idx, nb_idx, :) = C_sca;
                    
                    % Find peaks in absorption spectrum
                    [peaks, peak_locs] = findpeaks(C_abs, 'MinPeakHeight', max(C_abs)*0.1, ...
                                                   'MinPeakDistance', 20);
                    
                    % Store peak information
                    n_peaks = min(length(peaks), 2); % Store up to 2 peaks
                    if n_peaks > 0
                        [~, sort_idx] = sort(peaks, 'descend');
                        for p = 1:n_peaks
                            peak_idx = peak_locs(sort_idx(p));
                            peak_wavelengths(m_idx, f_idx, r_idx, nb_idx, p) = L(peak_idx);
                            peak_intensities(m_idx, f_idx, r_idx, nb_idx, p) = peaks(sort_idx(p));
                            
                            % Calculate quality factor (rough estimate)
                            peak_val = C_abs(peak_idx);
                            half_max = peak_val / 2;
                            
                            % Find FWHM
                            left_idx = find(C_abs(1:peak_idx) <= half_max, 1, 'last');
                            right_idx = find(C_abs(peak_idx:end) <= half_max, 1, 'first') + peak_idx - 1;
                            
                            if ~isempty(left_idx) && ~isempty(right_idx)
                                fwhm = L(right_idx) - L(left_idx);
                                peak_quality_factors(m_idx, f_idx, r_idx, nb_idx, p) = L(peak_idx) / fwhm;
                            end
                        end
                        
                        % Calculate resonance splitting if 2 peaks exist
                        if n_peaks == 2
                            resonance_splitting(m_idx, f_idx, r_idx, nb_idx) = ...
                                abs(peak_wavelengths(m_idx, f_idx, r_idx, nb_idx, 1) - ...
                                    peak_wavelengths(m_idx, f_idx, r_idx, nb_idx, 2));
                        end
                    end
                    
                catch ME
                    fprintf('\nWarning: Calculation failed for %s, f=%.1f, a=%.0fnm, nb=%.2f\n', ...
                            metal, f, a*1e3, nb);
                    fprintf('Error: %s\n', ME.message);
                end
            end
        end
    end
end

fprintf(']\nCalculations complete!\n');

%% Store Results
results.metals = metals;
results.core_fractions = core_fractions;
results.particle_radii = particle_radii;
results.background_indices = background_indices;
results.background_names = background_names;
results.wavelengths = L;
results.peak_wavelengths = peak_wavelengths;
results.peak_intensities = peak_intensities;
results.peak_quality_factors = peak_quality_factors;
results.resonance_splitting = resonance_splitting;
results.extinction_spectra = extinction_spectra;
results.absorption_spectra = absorption_spectra;
results.scattering_spectra = scattering_spectra;

%% Generate Summary Plots
fprintf('Generating summary plots...\n');

% 1. Peak wavelength vs core fraction for all metals (20nm, air)
figure('Position', [100, 100, 1200, 800]);
r_plot_idx = 2; % 20nm particles
nb_plot_idx = 1; % Air background

subplot(2,3,1);
colors = lines(n_metals);
for m_idx = 1:n_metals
    peak_data = squeeze(peak_wavelengths(m_idx, :, r_plot_idx, nb_plot_idx, 1));
    valid_idx = ~isnan(peak_data);
    if any(valid_idx)
        plot(core_fractions(valid_idx), peak_data(valid_idx)*1000, 'o-', ...
             'Color', colors(m_idx,:), 'LineWidth', 2, 'DisplayName', metals{m_idx});
        hold on;
    end
end
xlabel('Core Fraction f');
ylabel('Peak Wavelength (nm)');
title('Primary Resonance vs Shell Thickness');
legend('Location', 'best');
grid on;

% 2. Peak intensity vs core fraction
subplot(2,3,2);
for m_idx = 1:n_metals
    intensity_data = squeeze(peak_intensities(m_idx, :, r_plot_idx, nb_plot_idx, 1));
    valid_idx = ~isnan(intensity_data);
    if any(valid_idx)
        semilogy(core_fractions(valid_idx), intensity_data(valid_idx), 'o-', ...
                 'Color', colors(m_idx,:), 'LineWidth', 2, 'DisplayName', metals{m_idx});
        hold on;
    end
end
xlabel('Core Fraction f');
ylabel('Peak Intensity');
title('Resonance Strength vs Shell Thickness');
legend('Location', 'best');
grid on;

% 3. Resonance splitting vs core fraction
subplot(2,3,3);
for m_idx = 1:n_metals
    splitting_data = squeeze(resonance_splitting(m_idx, :, r_plot_idx, nb_plot_idx));
    valid_idx = ~isnan(splitting_data);
    if any(valid_idx)
        plot(core_fractions(valid_idx), splitting_data(valid_idx)*1000, 'o-', ...
             'Color', colors(m_idx,:), 'LineWidth', 2, 'DisplayName', metals{m_idx});
        hold on;
    end
end
xlabel('Core Fraction f');
ylabel('Resonance Splitting (nm)');
title('Mode Splitting vs Shell Thickness');
legend('Location', 'best');
grid on;

% 4. Size effect on peak wavelength (Ag, air, f=0.5)
subplot(2,3,4);
m_plot_idx = 1; % Silver
f_plot_idx = 5; % f = 0.5
for r_idx = 1:n_radii
    peak_data = squeeze(peak_wavelengths(m_plot_idx, f_plot_idx, r_idx, nb_plot_idx, 1));
    if ~isnan(peak_data)
        plot(particle_radii(r_idx)*1e3, peak_data*1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        hold on;
    end
end
xlabel('Particle Radius (nm)');
ylabel('Peak Wavelength (nm)');
title('Size Effect on Resonance (Ag, f=0.5)');
grid on;

% 5. Background effect on peak wavelength (Ag, 20nm, f=0.5)
subplot(2,3,5);
for nb_idx = 1:n_backgrounds
    peak_data = squeeze(peak_wavelengths(m_plot_idx, f_plot_idx, r_plot_idx, nb_idx, 1));
    if ~isnan(peak_data)
        bar(nb_idx, peak_data*1000, 'FaceColor', colors(nb_idx,:));
        hold on;
    end
end
set(gca, 'XTick', 1:n_backgrounds, 'XTickLabel', background_names);
ylabel('Peak Wavelength (nm)');
title('Background Effect on Resonance (Ag, 20nm, f=0.5)');
grid on;

% 6. Example absorption spectra (Ag, 20nm, air)
subplot(2,3,6);
f_examples = [1, 3, 5, 7, 9]; % f = 0.1, 0.3, 0.5, 0.7, 0.9
for i = 1:length(f_examples)
    f_idx = f_examples(i);
    spectrum = squeeze(absorption_spectra(m_plot_idx, f_idx, r_plot_idx, nb_plot_idx, :));
    plot(L*1000, spectrum, 'LineWidth', 2, ...
         'DisplayName', sprintf('f = %.1f', core_fractions(f_idx)));
    hold on;
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Absorption Spectra Evolution (Ag, 20nm)');
legend('Location', 'best');
grid on;

sgtitle('LSP Core-Shell Investigation Summary', 'FontSize', 16, 'FontWeight', 'bold');

%% Save Results
fprintf('Saving results...\n');
save('lsp_investigation_results.mat', 'results');

%% Generate Data Summary Report
fprintf('\n=== LSP INVESTIGATION SUMMARY ===\n');
fprintf('Metals investigated: %s\n', strjoin(metals, ', '));
fprintf('Core fractions: %.1f to %.1f (%.0f values)\n', ...
        min(core_fractions), max(core_fractions), length(core_fractions));
fprintf('Particle radii: %.0f to %.0f nm (%.0f values)\n', ...
        min(particle_radii)*1e3, max(particle_radii)*1e3, length(particle_radii));
fprintf('Background media: %s\n', strjoin(background_names, ', '));
fprintf('Total calculations: %d\n', total_calculations);
fprintf('Results saved to: lsp_investigation_results.mat\n');

%% Quick Analysis Tips
fprintf('\n=== ANALYSIS SUGGESTIONS ===\n');
fprintf('1. Examine peak_wavelengths array for trends vs core fraction\n');
fprintf('2. Look for resonance splitting patterns in different metals\n');
fprintf('3. Compare quality factors between metals\n');
fprintf('4. Analyze size-dependent effects on resonance position\n');
fprintf('5. Study background medium influence on spectral tuning\n');
fprintf('6. Plot field enhancement maps for interesting cases\n');

fprintf('\nInvestigation complete! Use the results structure for detailed analysis.\n');