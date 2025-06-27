% -----DATAS -----
clc; clear;close all
P1 = 100e3; % Initial pressure in Pa
T1 = 300; % Initial temperature in K
rc = 0.05; % Cut-off ratio
P3_P2 = 1.7; % Pressure ratio for constant volume
gamma = 1.4;
R = 287;
Tmax_limit = 2500; % Maximum allowable cycle temperature in K
expansion_ratio_atkinson = 17; % Atkinson expansion ratio
compression_ratios = 12:0.5:20;

% Initialize arrays
eff_dual = [];
cr_dual = [];

% ---------- Dual Cycle: Search for Optimal CR Under Tmax Constraint ----------
fprintf('Dual Cycle Under Tmax ≤ %.0f K:\n', Tmax_limit);
fprintf('CR | Efficiency (%%) | Tmax (K)\n');
for r = compression_ratios
    [Qin, Qout, eff, Tmax, ~, ~] = dual_cycle_analysis(r, P1, T1, gamma, R, rc, P3_P2);
    if Tmax <= Tmax_limit
        eff_dual(end+1) = eff * 100;
        cr_dual(end+1) = r;
        fprintf('%4.1f | %12.2f | %8.2f\n', r, eff * 100, Tmax);
    end
end

% Get optimal CR for Dual cycle
[max_eff, idx] = max(eff_dual);
optimal_CR = cr_dual(idx);
fprintf('\nOptimal Compression Ratio for Dual Cycle (Tmax ≤ %.0f K): %.1f\n', Tmax_limit, optimal_CR);
fprintf('Max Efficiency: %.2f%%\n\n', max_eff);

% ---------- Otto, Diesel, Dual Cycles at Optimal CR ----------
[Qin_dual, Qout_dual, eff_dual_at_opt, Tmax_dual, V_dual, P_dual] = dual_cycle_analysis(optimal_CR, P1, T1, gamma, R, rc, P3_P2);
Qin_best = Qin_dual;

[Qin_d, Qout_d, eff_d, Tmax_d, V_d, P_d] = diesel_cycle(optimal_CR, P1, T1, gamma, R, rc, Qin_best);
[Qin_o, Qout_o, eff_o, Tmax_o, V_o, P_o] = otto_cycle(optimal_CR, P1, T1, gamma, R, P3_P2, Qin_best);
[Qin_a, Qout_a, eff_a, Tmax_a, V_a, P_a] = atkinson_cycle(optimal_CR, expansion_ratio_atkinson, P1, T1, gamma, R, Qin_best);

% ---------- Display Summary ----------
fprintf('Cycle | Qin (kJ/kg) | Qout (kJ/kg) | Efficiency (%%) | Tmax (K)\n');
fprintf('---------|-------------|--------------|----------------|----------\n');
fprintf('Diesel | %11.2f | %12.2f | %14.2f | %8.2f\n', Qin_d/1000, Qout_d/1000, eff_d*100, Tmax_d);
fprintf('Otto | %11.2f | %12.2f | %14.2f | %8.2f\n', Qin_o/1000, Qout_o/1000, eff_o*100, Tmax_o);
fprintf('Dual | %11.2f | %12.2f | %14.2f | %8.2f\n', Qin_dual/1000, Qout_dual/1000, eff_dual_at_opt*100, Tmax_dual);
fprintf('Atkinson | %11.2f | %12.2f | %14.2f | %8.2f\n', Qin_a/1000, Qout_a/1000, eff_a*100, Tmax_a);

% ---------- P–V Diagrams ----------
figure; hold on;
plot(V_d, P_d, 'g-', 'LineWidth', 2, 'DisplayName', 'Diesel Cycle');
plot(V_o, P_o, 'r--', 'LineWidth', 2, 'DisplayName', 'Otto Cycle');
plot(V_dual, P_dual, 'b-.', 'LineWidth', 2, 'DisplayName', 'Dual Cycle');
xlabel('Volume [m^3]'); ylabel('Pressure [Pa]');
title(sprintf('P–V Diagram at Optimal CR = %.1f', optimal_CR));
legend; grid on;

figure; hold on;
plot(V_d, P_d, 'g-', 'LineWidth', 2, 'DisplayName', 'Diesel Cycle');
plot(V_o, P_o, 'r--', 'LineWidth', 2, 'DisplayName', 'Otto Cycle');
plot(V_dual, P_dual, 'b-.', 'LineWidth', 2, 'DisplayName', 'Dual Cycle');
plot(V_a, P_a, 'm:', 'LineWidth', 2, 'DisplayName', 'Atkinson Cycle');
xlabel('Volume [m^3]'); ylabel('Pressure [Pa]');
title(sprintf('P–V Diagram at Optimal CR = %.1f', optimal_CR));
legend; grid on;

% Individual P–V Diagrams
figure; plot(V_o, P_o, 'r--', 'LineWidth', 2); xlabel('Volume [m^3]'); ylabel('Pressure [Pa]'); title('Otto Cycle P–V Diagram'); grid on;
figure; plot(V_d, P_d, 'g-', 'LineWidth', 2); xlabel('Volume [m^3]'); ylabel('Pressure [Pa]'); title('Diesel Cycle P–V Diagram'); grid on;
figure; plot(V_dual, P_dual, 'b-.', 'LineWidth', 2); xlabel('Volume [m^3]'); ylabel('Pressure [Pa]'); title('Dual Cycle P–V Diagram'); grid on;
figure; plot(V_a, P_a, 'm:', 'LineWidth', 2); xlabel('Volume [m^3]'); ylabel('Pressure [Pa]'); title('Atkinson Cycle P–V Diagram'); grid on;

% ---------- Bar Chart: Efficiencies ----------
figure;
bar([eff_o, eff_d, eff_dual_at_opt, eff_a]*100);
set(gca, 'xticklabel', {'Otto', 'Diesel', 'Dual', 'Atkinson'});
ylabel('Efficiency [%]');
title(sprintf('Efficiency Comparison at CR = %.1f', optimal_CR));
grid on;

% ---------- Plot: Efficiency vs Compression Ratio (Dual) ----------
figure;
plot(cr_dual, eff_dual, 'b-o', 'LineWidth', 2);
xlabel('Compression Ratio'); ylabel('Efficiency [%]');
title(sprintf('Dual Cycle Efficiency vs CR (Tmax ≤ %.0f K)', Tmax_limit));
grid on;

% ---------------------- FUNCTIONS ----------------------
function [Qin, Qout, eff, Tmax, V, P] = dual_cycle_analysis(r, P1, T1, gamma, R, rc, P3_P2)
cv = R / (gamma - 1);
cp = gamma * cv;
v1 = 1; v2 = v1 / r;
T2 = T1 * r^(gamma - 1); P2 = P1 * r^gamma;
T3 = T2 * P3_P2; P3 = P2 * P3_P2;
V3 = v2; V4 = (v1 - V3) * rc + V3;
T4 = T3 * (V4 / V3);
T5 = T4 * (V4 / v1)^(gamma - 1); P5 = P1;
Qin = cv * (T3 - T2) + cp * (T4 - T3);
Qout = cv * (T5 - T1);
eff = 1 - Qout / Qin;
Tmax = max([T2, T3, T4, T5]);
V = [linspace(v1, v2, 50), repmat(v2,1,50), linspace(v2, V4,50), linspace(V4, v1,50), repmat(v1,1,50)];
P = [P1*(V(1:50)/v1).^(-gamma), linspace(P2,P3,50), repmat(P3,1,50), P3*(V(151:200)/V4).^(-gamma), linspace(P5,P1,50)];
end

function [Qin, Qout, eff, Tmax, V, P] = otto_cycle(r, P1, T1, gamma, R, P3_P2, Qin_best)
cv = R / (gamma - 1);
v1 = 1; v2 = v1 / r;
T2 = T1 * r^(gamma - 1); P2 = P1 * r^gamma;

if nargin < 7
    T3 = T2 * P3_P2;
else
    T3 = T2 + Qin_best / cv;
end

P3 = P2 * (T3 / T2);
T4 = T3 * (v2 / v1)^(gamma - 1); P4 = P3 * (v2 / v1)^gamma;
Qin = cv * (T3 - T2); Qout = cv * (T4 - T1);
eff = 1 - Qout / Qin;
Tmax = max([T2, T3, T4]);
V = [linspace(v1, v2, 50), repmat(v2,1,50), linspace(v2, v1,50), repmat(v1,1,50)];
P = [P1*(V(1:50)/v1).^(-gamma), linspace(P2,P3,50), P3*(V(101:150)/v2).^(-gamma), linspace(P4,P1,50)];
end

function [Qin, Qout, eff, Tmax, V, P] = diesel_cycle(r, P1, T1, gamma, R, rc, Qin_best)
cv = R / (gamma - 1); cp = gamma * cv;
v1 = 1; v2 = v1 / r;
T2 = T1 * r^(gamma - 1); P2 = P1 * r^gamma;

if nargin < 7
    V3 = (v1 - v2) * rc + v2;
    T3 = T2 * (V3 / v2);
else
    T3 = T2 + Qin_best / cp;
    V3 = v2 * (T3 / T2);
end

P3 = P2;
T4 = T3 * (V3 / v1)^(gamma - 1); P4 = P3 * (V3 / v1)^gamma;
Qin = cp * (T3 - T2); Qout = cv * (T4 - T1);
eff = 1 - Qout / Qin;
Tmax = max([T2, T3, T4]);
V = [linspace(v1, v2, 50), linspace(v2, V3,50), linspace(V3, v1,50), repmat(v1,1,50)];
P = [P1*(V(1:50)/v1).^(-gamma), repmat(P3,1,50), P3*(V(101:150)/V3).^(-gamma), linspace(P4,P1,50)];
end

function [Qin, Qout, eff, Tmax, V, P] = atkinson_cycle(r, re, P1, T1, gamma, R, Qin_input)
cv = R / (gamma - 1); cp = gamma * cv;
    v1 = 1; v2 = v1 / r;
    T2 = T1 * r^(gamma - 1); 
    P2 = P1 * r^gamma;
    T3 = T2 + Qin_input / cv;
    P3 = P2 * (T3 / T2);
    v4 = v2 * re;
    T4 = T3 * (v2 / v4)^(gamma - 1); 
    P4 = P3 * (v2 / v4)^gamma;
    P5 = P1; 
    T5 = T4 * (P5 / P4);  
    Qout = cv * (T4 - T5) + cp * (T5 - T1);
    Qin = Qin_input;  
    eff = 1 - Qout / Qin;
    Tmax = max([T2, T3, T4]);
    % P–V Data for plotting
    V = [linspace(v1, v2, 50), repmat(v2,1,50), linspace(v2, v4,50), repmat(v4,1,50)];
    P = [P1*(V(1:50)/v1).^(-gamma), ...
         P2 * (T2 + (0:49)/49 * (T3 - T2)) ./ T2, ...
         P3 * (v2 ./ V(101:150)).^gamma, ...
         linspace(P4, P1, 50)];
end

