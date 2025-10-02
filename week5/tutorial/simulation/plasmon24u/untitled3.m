% Compare metals
[L1, ~, C_abs1, ~] = calculate_lsp_spectrum(0.020, 0.3, 'Ag', 1.5, 1.5);
[L2, ~, C_abs2, ~] = calculate_lsp_spectrum(0.020, 0.3, 'Au', 1.5, 1.5);
[L3, ~, C_abs3, ~] = calculate_lsp_spectrum(0.020, 0.3, 'Al', 1.5, 1.5);
[L4, ~, C_abs4, ~] = calculate_lsp_spectrum(0.020, 0.3, 'Cu', 1.5, 1.5);
[L5, ~, C_abs5, ~] = calculate_lsp_spectrum(0.020, 0.3, 'K', 1.5, 1.5);


figure;
plot(L1*1000, C_abs1, 'b-', 'DisplayName', 'Ag');
hold on;
plot(L2*1000, C_abs2, 'r-', 'DisplayName', 'Au');
setup_spectrum_plot('Metal Comparison');


plot(L3*1000, C_abs3, 'g-', 'DisplayName', 'Al');
hold on;
plot(L4*1000, C_abs4, 'y-', 'DisplayName', 'Cu');
setup_spectrum_plot('Metal Comparison');

plot(L5*1000, C_abs5, 'c-', 'DisplayName', 'K');
hold on;

%% LSP Core-Shell Calculation Functions
% Modular functions extracted from example_mie_spectrum_3d_shell_eV.m
% Allows easy plotting of multiple cases on same figure

function [L, C_ext, C_abs, C_sca] = calculate_lsp_spectrum(a, f, metal, nc, nb, N)
    % Calculate LSP spectrum for given parameters
    %
    % Inputs:
    %   a - sphere radius (microns)
    %   f - fraction of volume occupied by metal
    %   metal - metal type string ('Ag', 'Au', 'Cu', 'Al', 'K')
    %   nc - core index
    %   nb - background index  
    %   N - number of Mie expansion orders
    %
    % Outputs:
    %   L - wavelength array (microns)
    %   C_ext, C_abs, C_sca - cross-section arrays
    
    if nargin < 6, N = 10; end
    
    % Get material wavelength bounds (following example)
    fm = ['eps', metal];
    [et, Lt] = feval(fm); % check material table
    Lmin = min(Lt); % ensure valid min L
    Lmin = min(Lt(real(et) < 0)); % only include et<0
    
    % Create energy/wavelength scale (following example)
    eV = 1.24 * linspace(0, 1, 500).' / Lmin; % use frequency scale (bounded)
    eV = eV(eV > 0); % avoid zero
    L = 1.24 ./ eV; % wavelength
    
    % Filter to reasonable range
    valid_range = (L >= 0.3) & (L <= 1.2);
    L = L(valid_range);
    
    % Calculate geometry (following example)
    av = a * [(1-f)^(1/3), 1]; % radii of surfaces
    
    % Initialize output arrays
    C_ext = zeros(size(L));
    C_abs = zeros(size(L));
    C_sca = zeros(size(L));
    
    % Main calculation loop (following example)
    for k = 1:length(L)
        L1 = L(k);
        
        try
            % Calculate material properties (following example)
            et_metal = feval(fm, L1);
            m = [nc, sqrt(et_metal), nb]; % refractive indices, inner to outer
            mu = ones(size(m)); % permeabilities
            
            % Calculate Mie coefficients (following example)
            coeffs = fmiecoeff(N, L1, av, m, mu); % expansion coefficients
            
            % Calculate cross-sections (following example)
            C1 = crosssection(coeffs, L1, av, m); % calculate cross-sections
            C_ext(k) = C1(1);
            C_abs(k) = C1(2);
            C_sca(k) = C1(3);
            
        catch
            % Handle calculation errors
            C_ext(k) = NaN;
            C_abs(k) = NaN;
            C_sca(k) = NaN;
        end
    end
end

function [peak_wavelengths, peak_intensities] = find_spectrum_peaks(L, spectrum)
    % Find peaks in spectrum
    %
    % Inputs:
    %   L - wavelength array (microns)
    %   spectrum - intensity array
    %
    % Outputs:
    %   peak_wavelengths - wavelengths of peaks (microns)
    %   peak_intensities - intensities of peaks
    
    % Remove NaN values
    valid_idx = ~isnan(spectrum);
    if sum(valid_idx) < 10
        peak_wavelengths = [];
        peak_intensities = [];
        return;
    end
    
    L_valid = L(valid_idx);
    spectrum_valid = spectrum(valid_idx);
    
    % Find peaks
    [peaks, peak_locs] = findpeaks(spectrum_valid, 'MinPeakHeight', max(spectrum_valid)*0.1, ...
                                   'MinPeakDistance', 10);
    
    if isempty(peaks)
        peak_wavelengths = [];
        peak_intensities = [];
    else
        peak_wavelengths = L_valid(peak_locs);
        peak_intensities = peaks;
    end
end

function plot_lsp_spectrum(L, C_ext, C_abs, C_sca, label_str, line_style)
    % Plot LSP spectrum with customizable styling
    %
    % Inputs:
    %   L - wavelength array (microns)
    %   C_ext, C_abs, C_sca - cross-section arrays
    %   label_str - label for legend
    %   line_style - line style struct with fields: Color, LineStyle, LineWidth
    
    if nargin < 6
        line_style = struct('Color', 'b', 'LineStyle', '-', 'LineWidth', 2);
    end
    
    % Convert wavelength to nm for plotting
    L_nm = L * 1000;
    
    % Plot extinction
    plot(L_nm, C_ext, 'Color', line_style.Color, 'LineStyle', line_style.LineStyle, ...
         'LineWidth', line_style.LineWidth, 'DisplayName', [label_str, ' (ext)']);
    hold on;
    
    % Plot absorption (dashed)
    plot(L_nm, C_abs, 'Color', line_style.Color, 'LineStyle', '--', ...
         'LineWidth', line_style.LineWidth, 'DisplayName', [label_str, ' (abs)']);
    
    % Plot scattering (dotted)  
    plot(L_nm, C_sca, 'Color', line_style.Color, 'LineStyle', ':', ...
         'LineWidth', line_style.LineWidth, 'DisplayName', [label_str, ' (sca)']);
end

function setup_spectrum_plot(title_str)
    % Setup standard spectrum plot formatting
    %
    % Input:
    %   title_str - plot title
    
    xlabel('Wavelength (nm)');
    ylabel('Cross-section');
    title(title_str);
    grid on;
    legend('Location', 'best');
    xlim([300, 1200]); % Standard visible range
end

%% Example usage and testing
if false % Set to true to run examples
    
    %% Example 1: Compare different metals
    figure;
    
    % Base parameters
    a = 0.020; % 20 nm
    f = 0.3;   % 30% metal shell
    nc = 1.5;  % silica core
    nb = 1.5;  % glass background
    
    metals = {'Ag', 'Au', 'Cu'};
    colors = {'b', 'r', 'g'};
    
    for i = 1:length(metals)
        [L, C_ext, C_abs, C_sca] = calculate_lsp_spectrum(a, f, metals{i}, nc, nb);
        
        line_style = struct('Color', colors{i}, 'LineStyle', '-', 'LineWidth', 2);
        plot_lsp_spectrum(L, C_ext, C_abs, C_sca, metals{i}, line_style);
    end
    
    setup_spectrum_plot('Metal Comparison (20nm, f=0.3)');
    
    %% Example 2: Compare different shell thicknesses for Ag
    figure;
    
    f_values = [0.1, 0.3, 0.5, 0.7, 0.9];
    colors = lines(length(f_values));
    
    for i = 1:length(f_values)
        [L, C_ext, C_abs, C_sca] = calculate_lsp_spectrum(a, f_values(i), 'Ag', nc, nb);
        
        line_style = struct('Color', colors(i,:), 'LineStyle', '-', 'LineWidth', 2);
        label_str = sprintf('f = %.1f', f_values(i));
        
        % Plot only absorption for clarity
        plot(L*1000, C_abs, 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', label_str);
        hold on;
    end
    
    setup_spectrum_plot('Shell Thickness Effect (Ag, 20nm)');
    
    %% Example 3: Peak analysis
    [L, C_ext, C_abs, C_sca] = calculate_lsp_spectrum(0.020, 0.3, 'Ag', 1.5, 1.5);
    [peak_waves, peak_ints] = find_spectrum_peaks(L, C_abs);
    
    fprintf('Found %d peaks:\n', length(peak_waves));
    for i = 1:length(peak_waves)
        fprintf('  Peak %d: %.0f nm, intensity %.2e\n', i, peak_waves(i)*1000, peak_ints(i));
    end
    
end%% LSP Systematic Investigation
% Based on example_mie_spectrum_3d_shell_eV.m structure
% Systematic study of core-shell plasmonic resonances

%% cleanup
clear
close all
clc

%% Investigation Parameters (following example structure)
% Base parameters from example
a_base = 0.020; % base sphere radius (in microns) = 20nm%% LSP Systematic Investigation
% Based on example_mie_spectrum_3d_shell_eV.m structure
% Systematic study of core-shell plasmonic resonances

%% cleanup
clear
close all
clc

%% Investigation Parameters (following example structure)
% Base parameters from example
a_base = 0.020; % base sphere radius (in microns) = 20nm
nc = 1.5; % core index (silica)
nb_base = 1.5; % base background index

% Parameters to vary systematically
metals = {'Ag', 'Au', 'Cu', 'Al', 'K'}; % All available metals
f_values = [0.1, 0.3, 0.5, 0.7, 0.9]; % Metal volume fractions to test
a_values = [0.010, 0.020, 0.030] * 1; % Particle sizes (10, 20, 30 nm)
nb_values = [1.0, 1.33, 1.5]; % Background indices (air, water, glass)
nb_names = {'Air', 'Water', 'Glass'};

% Fixed calculation parameters (from example)
N = 10; % number of orders

%% Energy/wavelength setup (following example approach)
% Check material availability and set wavelength bounds
[et_test, Lt_test] = epsAg(); % Get Ag data to determine bounds
Lmin = min(Lt_test(real(et_test) < 0)); % only include plasmonic region
Lmax = max(Lt_test);

% Create energy scale (following example)
n_points = 500; % Reduced from 1000 for faster calculation
eV = 1.24 * linspace(0.1, 1, n_points).' / Lmin;
eV = eV(eV > 0);
L = 1.24 ./ eV; % wavelength array

% Filter to reasonable range
valid_range = (L >= 0.3) & (L <= 1.2); % 300-1200 nm
L = L(valid_range);
eV = eV(valid_range);

fprintf('Wavelength range: %.0f - %.0f nm (%d points)\n', ...
        min(L)*1000, max(L)*1000, length(L));

%% Initialize storage
results = struct();
n_metals = length(metals);
n_f = length(f_values);
n_a = length(a_values);
n_nb = length(nb_values);
n_wavelengths = length(L);

% Pre-allocate arrays (following example output structure)
extinction_data = zeros(n_metals, n_f, n_a, n_nb, n_wavelengths);
absorption_data = zeros(n_metals, n_f, n_a, n_nb, n_wavelengths);
scattering_data = zeros(n_metals, n_f, n_a, n_nb, n_wavelengths);

% Peak analysis storage
peak_info = struct();
peak_info.wavelengths = {};
peak_info.intensities = {};
peak_info.indices = {};

%% Main calculation loop (following example structure)
fprintf('Starting systematic calculation...\n');
total_cases = n_metals * n_f * n_a * n_nb;
case_count = 0;

for metal_idx = 1:n_metals
    metal = metals{metal_idx};
    fprintf('\nProcessing %s:\n', metal);
    
    % Get material function handle (following example)
    fm = ['eps', metal];
    
    for f_idx = 1:n_f
        f = f_values(f_idx);
        
        for a_idx = 1:n_a
            a = a_values(a_idx);
            
            for nb_idx = 1:n_nb
                nb = nb_values(nb_idx);
                
                case_count = case_count + 1;
                fprintf('  Case %d/%d: f=%.1f, a=%.0fnm, nb=%.2f\n', ...
                        case_count, total_cases, f, a*1e3, nb);
                
                % Calculate geometry (following example)
                av = a * [(1-f)^(1/3), 1]; % radii of surfaces
                vol = 4*pi/3 * a^3;
                area = pi * a^2;
                
                % Initialize cross-section arrays
                C = nan(length(L), 3); % [ext, abs, sca]
                
                % Wavelength loop (following example structure)
                for k = 1:length(L)
                    L1 = L(k);
                    
                    try
                        % Get material properties (following example)
                        et_metal = feval(fm, L1);
                        m = [nc, sqrt(et_metal), nb]; % refractive indices
                        mu = ones(size(m)); % permeabilities
                        
                        % Calculate Mie coefficients (following example)
                        coeffs = fmiecoeff(N, L1, av, m, mu);
                        
                        % Calculate cross-sections (following example)
                        C1 = crosssection(coeffs, L1, av, m);
                        C(k, :) = C1;
                        
                    catch ME
                        % Handle calculation errors gracefully
                        C(k, :) = [NaN, NaN, NaN];
                    end
                end
                
                % Store results
                extinction_data(metal_idx, f_idx, a_idx, nb_idx, :) = C(:, 1);
                absorption_data(metal_idx, f_idx, a_idx, nb_idx, :) = C(:, 2);
                scattering_data(metal_idx, f_idx, a_idx, nb_idx, :) = C(:, 3);
                
                % Find peaks (following example approach)
                abs_spectrum = C(:, 2);
                valid_data = ~isnan(abs_spectrum);
                
                if sum(valid_data) > 10 % Need sufficient data points
                    try
                        [peaks, peak_locs] = findpeaks(abs_spectrum(valid_data), ...
                                                       'MinPeakHeight', max(abs_spectrum)*0.1, ...
                                                       'MinPeakDistance', 10);
                        if ~isempty(peaks)
                            % Store peak information
                            case_key = sprintf('%s_f%.1f_a%.0f_nb%.2f', ...
                                               metal, f, a*1e3, nb);
                            peak_info.wavelengths{end+1} = L(valid_data);
                            peak_info.intensities{end+1} = peaks;
                            peak_info.indices{end+1} = peak_locs;
                            peak_info.case_keys{end+1} = case_key;
                        end
                    catch
                        % Peak finding failed - continue
                    end
                end
            end
        end
    end
end

%% Store all results
results.metals = metals;
results.f_values = f_values;
results.a_values = a_values;
results.nb_values = nb_values;
results.nb_names = nb_names;
results.wavelengths = L;
results.extinction = extinction_data;
results.absorption = absorption_data;
results.scattering = scattering_data;
results.peak_info = peak_info;

%% Generate key analysis plots (following example visualization style)
fprintf('\nGenerating analysis plots...\n');

%% Plot 1: Shell thickness effect for different metals (base case: 20nm, glass)
figure('Position', [100, 100, 1400, 900]);

% Reference case indices
a_ref_idx = 2; % 20 nm
nb_ref_idx = 3; % glass background

subplot(2, 3, 1);
colors = lines(n_metals);
for metal_idx = 1:n_metals
    % Plot peak wavelength vs shell thickness
    peak_wavelengths = [];
    f_plot = [];
    
    for f_idx = 1:n_f
        spectrum = squeeze(absorption_data(metal_idx, f_idx, a_ref_idx, nb_ref_idx, :));
        if ~all(isnan(spectrum))
            [peak_val, peak_loc] = max(spectrum);
            if peak_val > 0
                peak_wavelengths(end+1) = L(peak_loc) * 1000; % Convert to nm
                f_plot(end+1) = f_values(f_idx);
            end
        end
    end
    
    if ~isempty(peak_wavelengths)
        plot(f_plot, peak_wavelengths, 'o-', 'Color', colors(metal_idx,:), ...
             'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', metals{metal_idx});
        hold on;
    end
end
xlabel('Metal Volume Fraction f');
ylabel('Peak Wavelength (nm)');
title('Resonance Tuning vs Shell Thickness');
legend('Location', 'best');
grid on;

%% Plot 2: Spectral comparison for one metal (Ag) at different shell thicknesses
subplot(2, 3, 2);
metal_plot_idx = 1; % Ag
for f_idx = 1:length(f_values)
    spectrum = squeeze(absorption_data(metal_plot_idx, f_idx, a_ref_idx, nb_ref_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, ...
             'DisplayName', sprintf('f = %.1f', f_values(f_idx)));
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Ag Spectra: Shell Thickness Effect');
legend('Location', 'best');
grid on;

%% Plot 3: Size effect comparison (Ag, f=0.5, glass)
subplot(2, 3, 3);
f_ref_idx = 3; % f = 0.5
size_colors = copper(n_a);
for a_idx = 1:n_a
    spectrum = squeeze(absorption_data(metal_plot_idx, f_ref_idx, a_idx, nb_ref_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, 'Color', size_colors(a_idx,:), ...
             'DisplayName', sprintf('%.0f nm', a_values(a_idx)*1e3));
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Size Effect (Ag, f=0.5)');
legend('Location', 'best');
grid on;

%% Plot 4: Background medium effect (Ag, f=0.5, 20nm)
subplot(2, 3, 4);
bg_colors = autumn(n_nb);
for nb_idx = 1:n_nb
    spectrum = squeeze(absorption_data(metal_plot_idx, f_ref_idx, a_ref_idx, nb_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, 'Color', bg_colors(nb_idx,:), ...
             'DisplayName', nb_names{nb_idx});
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Background Medium Effect (Ag, f=0.5, 20nm)');
legend('Location', 'best');
grid on;

%% Plot 5: Metal comparison (f=0.5, 20nm, glass)
subplot(2, 3, 5);
for metal_idx = 1:n_metals
    spectrum = squeeze(absorption_data(metal_idx, f_ref_idx, a_ref_idx, nb_ref_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, 'Color', colors(metal_idx,:), ...
             'DisplayName', metals{metal_idx});
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Metal Comparison (f=0.5, 20nm, Glass)');
legend('Location', 'best');
grid on;

%% Plot 6: Peak intensity heatmap for Ag
subplot(2, 3, 6);
peak_matrix = zeros(n_f, n_a);
for f_idx = 1:n_f
    for a_idx = 1:n_a
        spectrum = squeeze(absorption_data(metal_plot_idx, f_idx, a_idx, nb_ref_idx, :));
        if ~all(isnan(spectrum))
            peak_matrix(f_idx, a_idx) = max(spectrum);
        end
    end
end

imagesc(a_values*1e3, f_values, peak_matrix);
colorbar;
xlabel('Particle Radius (nm)');
ylabel('Metal Volume Fraction f');
title('Peak Absorption Intensity (Ag, Glass)');
set(gca, 'YDir', 'normal');

sgtitle('LSP Core-Shell Systematic Investigation', 'FontSize', 16, 'FontWeight', 'bold');

%% Save results
fprintf('Saving results...\n');
save('lsp_systematic_results.mat', 'results');

%% Summary report
fprintf('\n=== INVESTIGATION COMPLETE ===\n');
fprintf('Cases calculated: %d\n', case_count);
fprintf('Metals: %s\n', strjoin(metals, ', '));
fprintf('Shell fractions: [%.1f-%.1f]\n', min(f_values), max(f_values));
fprintf('Particle sizes: [%.0f-%.0f] nm\n', min(a_values)*1e3, max(a_values)*1e3);
fprintf('Background media: %s\n', strjoin(nb_names, ', '));
fprintf('Results saved to: lsp_systematic_results.mat\n');
fprintf('\nAnalyze results.absorption, results.extinction, results.scattering arrays\n');
fprintf('Dimensions: [metal, f, size, background, wavelength]
nc = 1.5; % core index (silica)
nb_base = 1.5; % base background index

% Parameters to vary systematically
metals = {'Ag', 'Au', 'Cu', 'Al', 'K'}; % All available metals
f_values = [0.1, 0.3, 0.5, 0.7, 0.9]; % Metal volume fractions to test
a_values = [0.010, 0.020, 0.030] * 1; % Particle sizes (10, 20, 30 nm)
nb_values = [1.0, 1.33, 1.5]; % Background indices (air, water, glass)
nb_names = {'Air', 'Water', 'Glass'};

% Fixed calculation parameters (from example)
N = 10; % number of orders

%% Energy/wavelength setup (following example approach)
% Check material availability and set wavelength bounds
[et_test, Lt_test] = epsAg(); % Get Ag data to determine bounds
Lmin = min(Lt_test(real(et_test) < 0)); % only include plasmonic region
Lmax = max(Lt_test);

% Create energy scale (following example)
n_points = 500; % Reduced from 1000 for faster calculation
eV = 1.24 * linspace(0.1, 1, n_points).' / Lmin;
eV = eV(eV > 0);
L = 1.24 ./ eV; % wavelength array

% Filter to reasonable range
valid_range = (L >= 0.3) & (L <= 1.2); % 300-1200 nm
L = L(valid_range);
eV = eV(valid_range);

fprintf('Wavelength range: %.0f - %.0f nm (%d points)\n', ...
        min(L)*1000, max(L)*1000, length(L));

%% Initialize storage
results = struct();
n_metals = length(metals);
n_f = length(f_values);
n_a = length(a_values);
n_nb = length(nb_values);
n_wavelengths = length(L);

% Pre-allocate arrays (following example output structure)
extinction_data = zeros(n_metals, n_f, n_a, n_nb, n_wavelengths);
absorption_data = zeros(n_metals, n_f, n_a, n_nb, n_wavelengths);
scattering_data = zeros(n_metals, n_f, n_a, n_nb, n_wavelengths);

% Peak analysis storage
peak_info = struct();
peak_info.wavelengths = {};
peak_info.intensities = {};
peak_info.indices = {};

%% Main calculation loop (following example structure)
fprintf('Starting systematic calculation...\n');
total_cases = n_metals * n_f * n_a * n_nb;
case_count = 0;

for metal_idx = 1:n_metals
    metal = metals{metal_idx};
    fprintf('\nProcessing %s:\n', metal);
    
    % Get material function handle (following example)
    fm = ['eps', metal];
    
    for f_idx = 1:n_f
        f = f_values(f_idx);
        
        for a_idx = 1:n_a
            a = a_values(a_idx);
            
            for nb_idx = 1:n_nb
                nb = nb_values(nb_idx);
                
                case_count = case_count + 1;
                fprintf('  Case %d/%d: f=%.1f, a=%.0fnm, nb=%.2f\n', ...
                        case_count, total_cases, f, a*1e3, nb);
                
                % Calculate geometry (following example)
                av = a * [(1-f)^(1/3), 1]; % radii of surfaces
                vol = 4*pi/3 * a^3;
                area = pi * a^2;
                
                % Initialize cross-section arrays
                C = nan(length(L), 3); % [ext, abs, sca]
                
                % Wavelength loop (following example structure)
                for k = 1:length(L)
                    L1 = L(k);
                    
                    try
                        % Get material properties (following example)
                        et_metal = feval(fm, L1);
                        m = [nc, sqrt(et_metal), nb]; % refractive indices
                        mu = ones(size(m)); % permeabilities
                        
                        % Calculate Mie coefficients (following example)
                        coeffs = fmiecoeff(N, L1, av, m, mu);
                        
                        % Calculate cross-sections (following example)
                        C1 = crosssection(coeffs, L1, av, m);
                        C(k, :) = C1;
                        
                    catch ME
                        % Handle calculation errors gracefully
                        C(k, :) = [NaN, NaN, NaN];
                    end
                end
                
                % Store results
                extinction_data(metal_idx, f_idx, a_idx, nb_idx, :) = C(:, 1);
                absorption_data(metal_idx, f_idx, a_idx, nb_idx, :) = C(:, 2);
                scattering_data(metal_idx, f_idx, a_idx, nb_idx, :) = C(:, 3);
                
                % Find peaks (following example approach)
                abs_spectrum = C(:, 2);
                valid_data = ~isnan(abs_spectrum);
                
                if sum(valid_data) > 10 % Need sufficient data points
                    try
                        [peaks, peak_locs] = findpeaks(abs_spectrum(valid_data), ...
                                                       'MinPeakHeight', max(abs_spectrum)*0.1, ...
                                                       'MinPeakDistance', 10);
                        if ~isempty(peaks)
                            % Store peak information
                            case_key = sprintf('%s_f%.1f_a%.0f_nb%.2f', ...
                                               metal, f, a*1e3, nb);
                            peak_info.wavelengths{end+1} = L(valid_data);
                            peak_info.intensities{end+1} = peaks;
                            peak_info.indices{end+1} = peak_locs;
                            peak_info.case_keys{end+1} = case_key;
                        end
                    catch
                        % Peak finding failed - continue
                    end
                end
            end
        end
    end
end

%% Store all results
results.metals = metals;
results.f_values = f_values;
results.a_values = a_values;
results.nb_values = nb_values;
results.nb_names = nb_names;
results.wavelengths = L;
results.extinction = extinction_data;
results.absorption = absorption_data;
results.scattering = scattering_data;
results.peak_info = peak_info;

%% Generate key analysis plots (following example visualization style)
fprintf('\nGenerating analysis plots...\n');

%% Plot 1: Shell thickness effect for different metals (base case: 20nm, glass)
figure('Position', [100, 100, 1400, 900]);

% Reference case indices
a_ref_idx = 2; % 20 nm
nb_ref_idx = 3; % glass background

subplot(2, 3, 1);
colors = lines(n_metals);
for metal_idx = 1:n_metals
    % Plot peak wavelength vs shell thickness
    peak_wavelengths = [];
    f_plot = [];
    
    for f_idx = 1:n_f
        spectrum = squeeze(absorption_data(metal_idx, f_idx, a_ref_idx, nb_ref_idx, :));
        if ~all(isnan(spectrum))
            [peak_val, peak_loc] = max(spectrum);
            if peak_val > 0
                peak_wavelengths(end+1) = L(peak_loc) * 1000; % Convert to nm
                f_plot(end+1) = f_values(f_idx);
            end
        end
    end
    
    if ~isempty(peak_wavelengths)
        plot(f_plot, peak_wavelengths, 'o-', 'Color', colors(metal_idx,:), ...
             'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', metals{metal_idx});
        hold on;
    end
end
xlabel('Metal Volume Fraction f');
ylabel('Peak Wavelength (nm)');
title('Resonance Tuning vs Shell Thickness');
legend('Location', 'best');
grid on;

%% Plot 2: Spectral comparison for one metal (Ag) at different shell thicknesses
subplot(2, 3, 2);
metal_plot_idx = 1; % Ag
for f_idx = 1:length(f_values)
    spectrum = squeeze(absorption_data(metal_plot_idx, f_idx, a_ref_idx, nb_ref_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, ...
             'DisplayName', sprintf('f = %.1f', f_values(f_idx)));
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Ag Spectra: Shell Thickness Effect');
legend('Location', 'best');
grid on;

%% Plot 3: Size effect comparison (Ag, f=0.5, glass)
subplot(2, 3, 3);
f_ref_idx = 3; % f = 0.5
size_colors = copper(n_a);
for a_idx = 1:n_a
    spectrum = squeeze(absorption_data(metal_plot_idx, f_ref_idx, a_idx, nb_ref_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, 'Color', size_colors(a_idx,:), ...
             'DisplayName', sprintf('%.0f nm', a_values(a_idx)*1e3));
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Size Effect (Ag, f=0.5)');
legend('Location', 'best');
grid on;

%% Plot 4: Background medium effect (Ag, f=0.5, 20nm)
subplot(2, 3, 4);
bg_colors = autumn(n_nb);
for nb_idx = 1:n_nb
    spectrum = squeeze(absorption_data(metal_plot_idx, f_ref_idx, a_ref_idx, nb_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, 'Color', bg_colors(nb_idx,:), ...
             'DisplayName', nb_names{nb_idx});
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Background Medium Effect (Ag, f=0.5, 20nm)');
legend('Location', 'best');
grid on;

%% Plot 5: Metal comparison (f=0.5, 20nm, glass)
subplot(2, 3, 5);
for metal_idx = 1:n_metals
    spectrum = squeeze(absorption_data(metal_idx, f_ref_idx, a_ref_idx, nb_ref_idx, :));
    if ~all(isnan(spectrum))
        plot(L*1000, spectrum, 'LineWidth', 2, 'Color', colors(metal_idx,:), ...
             'DisplayName', metals{metal_idx});
        hold on;
    end
end
xlabel('Wavelength (nm)');
ylabel('Absorption Cross-section');
title('Metal Comparison (f=0.5, 20nm, Glass)');
legend('Location', 'best');
grid on;

%% Plot 6: Peak intensity heatmap for Ag
subplot(2, 3, 6);
peak_matrix = zeros(n_f, n_a);
for f_idx = 1:n_f
    for a_idx = 1:n_a
        spectrum = squeeze(absorption_data(metal_plot_idx, f_idx, a_idx, nb_ref_idx, :));
        if ~all(isnan(spectrum))
            peak_matrix(f_idx, a_idx) = max(spectrum);
        end
    end
end

imagesc(a_values*1e3, f_values, peak_matrix);
colorbar;
xlabel('Particle Radius (nm)');
ylabel('Metal Volume Fraction f');
title('Peak Absorption Intensity (Ag, Glass)');
set(gca, 'YDir', 'normal');

sgtitle('LSP Core-Shell Systematic Investigation', 'FontSize', 16, 'FontWeight', 'bold');

%% Save results
fprintf('Saving results...\n');
save('lsp_systematic_results.mat', 'results');

%% Summary report
fprintf('\n=== INVESTIGATION COMPLETE ===\n');
fprintf('Cases calculated: %d\n', case_count);
fprintf('Metals: %s\n', strjoin(metals, ', '));
fprintf('Shell fractions: [%.1f-%.1f]\n', min(f_values), max(f_values));
fprintf('Particle sizes: [%.0f-%.0f] nm\n', min(a_values)*1e3, max(a_values)*1e3);
fprintf('Background media: %s\n', strjoin(nb_names, ', '));
fprintf('Results saved to: lsp_systematic_results.mat\n');
fprintf('\nAnalyze results.absorption, results.extinction, results.scattering arrays\n');
fprintf('Dimensions: [metal, f, size, background, wavelength');