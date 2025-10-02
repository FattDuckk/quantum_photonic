compare_metals_easy()  % Uses defaults: 20nm, f=0.3
compare_metals_easy(0.020, 0.5)  % 20nm particles, f=0.5
%% LSP Easy-to-Use Functions
% Simplified functions for easy plotting and comparison

function [L, C_ext, C_abs, C_sca, eV, Q_abs] = calculate_lsp_spectrum(a, f, metal, nc, nb, N)
    % Calculate LSP spectrum for given parameters
    %
    % Inputs:
    %   a - sphere radius (microns) 
    %   f - fraction of volume occupied by metal
    %   metal - metal type string ('Ag', 'Au', 'Cu', 'Al', 'K')
    %   nc - core index (default 1.5)
    %   nb - background index (default 1.5)
    %   N - number of Mie expansion orders (default 10)
    %
    % Outputs:
    %   L - wavelength array (microns)
    %   C_ext, C_abs, C_sca - cross-section arrays
    %   eV - photon energy array (eV)
    %   Q_abs - normalized absorption efficiency C/V/k (dimensionless)
    
    % Set defaults
    if nargin < 4, nc = 1.5; end
    if nargin < 5, nb = 1.5; end  
    if nargin < 6, N = 10; end
    
    % Get material wavelength bounds (following example)
    fm = ['eps', metal];
    [et, Lt] = feval(fm); % check material table
    Lmin = min(Lt); % ensure valid min L
    Lmin = min(Lt(real(et) < 0)); % only include et<0
    
    % Create energy/wavelength scale (following example exactly)
    eV = 1.24 * linspace(0, 1, 500).' / Lmin; % use frequency scale (bounded)
    eV = eV(eV > 0); % avoid zero
    L = 1.24 ./ eV; % wavelength
    
    % Filter to reasonable range
    valid_range = (L >= 0.3) & (L <= 1.2);
    L = L(valid_range);
    eV = eV(valid_range);
    
    % Calculate geometry (following example)
    av = a * [(1-f)^(1/3), 1]; % radii of surfaces
    vol = 4*pi/3 * a^3; % total volume
    
    % Initialize output arrays
    C_ext = zeros(size(L));
    C_abs = zeros(size(L));
    C_sca = zeros(size(L));
    
    % Main calculation loop (following example exactly)
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
    
    % Calculate C/V/k normalization (following example)
    k_wave = 2*pi ./ L; % wavenumber
    Q_abs = C_abs ./ vol ./ k_wave; % normalized efficiency (dimensionless)
end

function compare_metals_easy(a, f, nc, nb, plot_type)
    % Easy function to compare all metals
    %
    % Inputs:
    %   a - particle radius (microns), default 0.020 (20nm)
    %   f - metal fraction, default 0.3  
    %   nc - core index, default 1.5
    %   nb - background index, default 1.5
    %   plot_type - 'wavelength' or 'energy' or 'both', default 'both'
    
    % Set defaults
    if nargin < 1, a = 0.020; end
    if nargin < 2, f = 0.3; end
    if nargin < 3, nc = 1.5; end
    if nargin < 4, nb = 1.5; end
    if nargin < 5, plot_type = 'both'; end
    
    metals = {'Ag', 'Au', 'Cu', 'Al', 'K'};
    colors = lines(length(metals));
    
    % Calculate all spectra
    data = cell(length(metals), 1);
    for i = 1:length(metals)
        [L, C_ext, C_abs, C_sca, eV, Q_abs] = calculate_lsp_spectrum(a, f, metals{i}, nc, nb);
        data{i} = struct('L', L, 'C_ext', C_ext, 'C_abs', C_abs, 'C_sca', C_sca, ...
                        'eV', eV, 'Q_abs', Q_abs, 'metal', metals{i});
    end
    
    % Create plots
    if strcmp(plot_type, 'wavelength') || strcmp(plot_type, 'both')
        figure('Position', [100, 100, 1200, 400]);
        
        % Wavelength plot - Absorption cross-section
        subplot(1, 2, 1);
        for i = 1:length(metals)
            plot(data{i}.L*1000, data{i}.C_abs, 'Color', colors(i,:), ...
                 'LineWidth', 2, 'DisplayName', metals{i});
            hold on;
        end
        xlabel('Wavelength (nm)');
        ylabel('Absorption Cross-section');
        title(sprintf('Metal Comparison: Cross-section\na=%.0fnm, f=%.1f', a*1000, f));
        legend('Location', 'best');
        grid on;
        xlim([300, 1200]);
        
        % Wavelength plot - Normalized efficiency C/V/k
        subplot(1, 2, 2);
        for i = 1:length(metals)
            plot(data{i}.L*1000, data{i}.Q_abs, 'Color', colors(i,:), ...
                 'LineWidth', 2, 'DisplayName', metals{i});
            hold on;
        end
        xlabel('Wavelength (nm)');
        ylabel('C/V/k (dimensionless)');
        title('Metal Comparison: Normalized Efficiency');
        legend('Location', 'best');
        grid on;
        xlim([300, 1200]);
    end
    
    if strcmp(plot_type, 'energy') || strcmp(plot_type, 'both')
        figure('Position', [100, 600, 1200, 400]);
        
        % Energy plot - Absorption cross-section  
        subplot(1, 2, 1);
        for i = 1:length(metals)
            plot(data{i}.eV, data{i}.C_abs, 'Color', colors(i,:), ...
                 'LineWidth', 2, 'DisplayName', metals{i});
            hold on;
        end
        xlabel('Photon Energy (eV)');
        ylabel('Absorption Cross-section');
        title(sprintf('Metal Comparison vs Energy\na=%.0fnm, f=%.1f', a*1000, f));
        legend('Location', 'best');
        grid on;
        
        % Energy plot - Normalized efficiency C/V/k
        subplot(1, 2, 2);
        for i = 1:length(metals)
            plot(data{i}.eV, data{i}.Q_abs, 'Color', colors(i,:), ...
                 'LineWidth', 2, 'DisplayName', metals{i});
            hold on;
        end
        xlabel('Photon Energy (eV)');
        ylabel('C/V/k (dimensionless)');
        title('Normalized Efficiency vs Energy');
        legend('Location', 'best');
        grid on;
    end
end

function compare_shell_thickness(metal, a, nc, nb)
    % Easy function to compare shell thicknesses for one metal
    %
    % Inputs:
    %   metal - metal type ('Ag', 'Au', etc.), default 'Ag'
    %   a - particle radius (microns), default 0.020
    %   nc - core index, default 1.5  
    %   nb - background index, default 1.5
    
    % Set defaults
    if nargin < 1, metal = 'Ag'; end
    if nargin < 2, a = 0.020; end
    if nargin < 3, nc = 1.5; end
    if nargin < 4, nb = 1.5; end
    
    f_values = [0.1, 0.3, 0.5, 0.7, 0.9];
    colors = parula(length(f_values));
    
    figure('Position', [100, 100, 1200, 400]);
    
    % Wavelength plot
    subplot(1, 2, 1);
    for i = 1:length(f_values)
        [L, ~, C_abs, ~] = calculate_lsp_spectrum(a, f_values(i), metal, nc, nb);
        plot(L*1000, C_abs, 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('f = %.1f', f_values(i)));
        hold on;
    end
    xlabel('Wavelength (nm)');
    ylabel('Absorption Cross-section');
    title(sprintf('%s Shell Thickness Effect\na=%.0fnm', metal, a*1000));
    legend('Location', 'best');
    grid on;
    
    % Peak wavelength vs shell thickness
    subplot(1, 2, 2);
    peak_wavelengths = [];
    for i = 1:length(f_values)
        [L, ~, C_abs, ~] = calculate_lsp_spectrum(a, f_values(i), metal, nc, nb);
        [~, peak_idx] = max(C_abs);
        if ~isnan(C_abs(peak_idx))
            peak_wavelengths(i) = L(peak_idx) * 1000;
        else
            peak_wavelengths(i) = NaN;
        end
    end
    
    plot(f_values, peak_wavelengths, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Metal Volume Fraction f');
    ylabel('Peak Wavelength (nm)');
    title('Resonance Tuning');
    grid on;
end

%% Quick examples for testing
fprintf('=== LSP Easy Functions Loaded ===\n');
fprintf('Quick usage examples:\n');
fprintf('1. compare_metals_easy()  %% Compare all metals with defaults\n');
fprintf('2. compare_metals_easy(0.020, 0.5)  %% 20nm, f=0.5\n');
fprintf('3. compare_shell_thickness(''Ag'')  %% Ag shell thickness effect\n');
fprintf('4. [L,~,C_abs] = calculate_lsp_spectrum(0.020, 0.3, ''Au'')  %% Single calc\n');

%% Run a quick demo
if false % Set to true to run demo
    fprintf('\nRunning demo...\n');
    compare_metals_easy();
    compare_shell_thickness('Ag');
end