clear; clc;
% NOx Formation Simulation for Different Engine Cycles
% Constants
R = 8.314; % J/mol.K
% Initial mole fractions
N2_0 = 0.79;
O2_0 = 0.21 * 0.07;
O_0 = 1e-10; 
N_0 = 0;
NO_0 = 0;
% Time array (0 to 4 ms)
t = linspace(0, 0.004, 1000);
% Temperature profiles for each cycle
T_profiles.Dual = @(t_scalar) 2489.61 * (0.1152 / (0.1152 + (1 - 0.1152) * (t_scalar / 0.005))) ^ 0.4;
T_profiles.Otto = @(t_scalar) 2890.93 * (0.0686 / (0.0686 + (1 - 0.0686) * (t_scalar / 0.005))) ^ 0.4;
T_profiles.Diesel = @(t_scalar) 2314.75 * (0.18 / (0.18 + (1 - 0.18) * (t_scalar / 0.005))) ^ 0.4;
T_profiles.Atkinson = @(t_scalar) 2890.93 * (0.0686 / (0.0686 + (1.17 - 0.0686) * (t_scalar / 0.005))) ^ 0.4;

% Reaction rate function
reactionRates = @(y, t, Tfunc) [
    -1.47e13 * (Tfunc(t)^0.3) * exp(-75286.81/(R*Tfunc(t))) * y(1) * y(3);       % dN2/dt
    -6.4e9 * Tfunc(t) * exp(-6285.5/(R*Tfunc(t))) * y(4) * y(2);                 % dO2/dt
    -1.47e13 * (Tfunc(t)^0.3) * exp(-75286.81/(R*Tfunc(t))) * y(1) * y(3) + ...
     6.4e9 * Tfunc(t) * exp(-6285.5/(R*Tfunc(t))) * y(4) * y(2);                 % dO/dt
     1.47e13 * (Tfunc(t)^0.3) * exp(-75286.81/(R*Tfunc(t))) * y(1) * y(3) - ...
     6.4e9 * Tfunc(t) * exp(-6285.5/(R*Tfunc(t))) * y(4) * y(2);                 % dN/dt
     1.47e13 * (Tfunc(t)^0.3) * exp(-75286.81/(R*Tfunc(t))) * y(1) * y(3) + ...
     6.4e9 * Tfunc(t) * exp(-6285.5/(R*Tfunc(t))) * y(4) * y(2)                  % dNO/dt
];

% Simulation for all cycles
results = struct();
cycles = fieldnames(T_profiles);
final_ppms = zeros(length(cycles),1);

fprintf('Final NO Mole Fractions:\n');
for i = 1:length(cycles)
    name = cycles{i};
    y0 = [N2_0, O2_0, O_0, N_0, NO_0];
    Tfunc = T_profiles.(name);
    [T_sim, Y] = ode15s(@(t,y) reactionRates(y,t,Tfunc), t, y0);
    NO_conc = Y(:,5);
    results.(name).t = T_sim;
    results.(name).NO = NO_conc;
    results.(name).final_ppm = NO_conc(end) * 1e6;
    final_ppms(i) = results.(name).final_ppm;
    fprintf('%s: %.2f ppm\n', name, results.(name).final_ppm);
end

%% Plot 1: NO Concentration Over Time
figure;
hold on;
colors = lines(length(cycles));
for i = 1:length(cycles)
    plot(results.(cycles{i}).t * 1000, results.(cycles{i}).NO * 1e6, ...
        'DisplayName', [cycles{i} ' Cycle'], 'Color', colors(i,:));
end
xlabel('Time (ms)'); ylabel('[NO] (ppm)');
title('NO Formation Across Engine Cycles');
grid on; legend; hold off;

%% Plot 2: Temperature Profiles
figure;
for i = 1:length(cycles)
    subplot(2,2,i);
    T_vals = arrayfun(T_profiles.(cycles{i}), t);
    plot(t*1000, T_vals, 'b');
    title([cycles{i} ' Cycle']);
    xlabel('Time (ms)'); ylabel('Temperature (K)');
    grid on;
end
sgtitle('Temperature Profiles of Engine Cycles');

%% Plot 3: Contour Plot of NO vs. CR and Tpeak
CR = linspace(8, 20, 100);        % Compression Ratio range
Tpeak = linspace(2000, 3000, 100);% Peak temperature range
[CR_grid, T_grid] = meshgrid(CR, Tpeak);
a = 1e-6; b = 0.002; c = 1.5;
NO_map = a * exp(b .* T_grid) .* CR_grid .^ c;

figure;
contourf(CR_grid, T_grid, NO_map, 20, 'LineColor','none');
colorbar; colormap hot;
xlabel('Compression Ratio'); ylabel('Peak Temperature (K)');
title('NO Concentration vs. Compression Ratio and Peak Temperature');
grid on;

%% Plot 4: Bar Chart of Final NO Emissions
figure;
bar(final_ppms);
set(gca, 'XTickLabel', cycles, 'FontSize', 10);
ylabel('Final [NO] (ppm)');
title('Comparison of Final NO Emissions by Cycle');
grid on;

%% Plot 5: 3D Surface Plot of NO vs CR and Tpeak
figure;
surf(CR_grid, T_grid, NO_map, 'EdgeColor', 'none');
xlabel('Compression Ratio'); ylabel('Peak Temperature (K)'); zlabel('[NO] ppm');
title('NO Formation Surface');
colormap turbo;
colorbar;
view(45, 30);

