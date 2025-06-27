clc; clear; close all
% === Thermodynamic Parameters ===
pres_ref = 105e3;          % Reference pressure [Pa]
temp_ref = 300;            % Reference temperature [K]
rat_spec = 1.4;            % Specific heat ratio (gamma)
gas_const = 287;           % Gas constant [J/kg·K]
area_throat = 5e-4;        % Effective flow area [m^2]
temp_up = 360;             % Upstream total temperature [K]

% === Pressure Range ===
press_start = 100e3;       % Starting upstream pressure [Pa]
press_end = 500e3;         % Ending upstream pressure [Pa]
N_steps = 500;
P_range = linspace(press_start, press_end, N_steps);
P_kPa = P_range / 1000;    % Convert to kPa for plots

% === Critical Pressure Ratio for Sonic Flow ===
rat_crit = ((rat_spec + 1)/2)^(rat_spec / (rat_spec - 1));

% === Preallocate Results ===
vel_exit = zeros(size(P_range));
mdot_result = zeros(size(P_range));

% === Calculate Mass Flow Rate & Velocity ===
for n = 1:length(P_range)
    p_ratio = P_range(n) / pres_ref;
    
    if p_ratio < rat_crit
        % Subsonic Flow
        vel_exit(n) = sqrt(2 * gas_const * temp_up * (rat_spec / (rat_spec - 1)) * ...
                          (1 - (pres_ref / P_range(n))^((rat_spec - 1)/rat_spec)));
                      
        mdot_result(n) = area_throat * P_range(n) * sqrt((2 / (gas_const * temp_up)) * ...
                          (pres_ref / P_range(n))^(2/rat_spec) * ...
                          (rat_spec / (rat_spec - 1)) * ...
                          (1 - (pres_ref / P_range(n))^((rat_spec - 1)/rat_spec)));
    else
        % Sonic Flow (choked)
        vel_exit(n) = sqrt(rat_spec * gas_const * temp_ref);
        
        mdot_result(n) = area_throat * P_range(n) * sqrt((rat_spec / (gas_const * temp_up)) * ...
                          (2 / (rat_spec + 1))^((rat_spec + 1)/(rat_spec - 1)));
    end
end

% === Identify Sonic Transition Point ===
p_trans_index = find(P_range / pres_ref >= rat_crit, 1, 'first');
p_sonic_transition = P_kPa(p_trans_index);

% === Plot: Exit Velocity vs Inlet Pressure ===
figure;
plot(P_kPa, vel_exit, 'b', 'LineWidth', 2); hold on;
xline(p_sonic_transition, '--k', sprintf('Sonic Transition ≈ %.1f kPa', p_sonic_transition), ...
       'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LineWidth',1.2);
xlabel('Inlet Pressure [kPa]');
ylabel('Exit Velocity [m/s]');
title('Exit Velocity vs Inlet Pressure with Sonic Transition');
legend('Exit Velocity','Sonic Limit');
grid on;

% === Plot: Mass Flow Rate vs Inlet Pressure ===
figure;
plot(P_kPa, mdot_result, 'r', 'LineWidth', 2); hold on;
xline(p_sonic_transition, '--k', sprintf('Sonic Transition ≈ %.1f kPa', p_sonic_transition), ...
       'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LineWidth',1.2);
xlabel('Inlet Pressure [kPa]');
ylabel('Mass Flow Rate [kg/s]');
title('Mass Flow Rate vs Inlet Pressure with Flow Regime Indication');
legend('Mass Flow Rate','Sonic Limit');
grid on;

% === Valve Lift Profile Setup ===
angle_vals = linspace(0, 720, 1000);  

% Intake Valve: 340° → 580°
intake_lift = zeros(size(angle_vals));
intake_open = 340; intake_close = 580;
intake_duration = intake_close - intake_open;
for i = 1:length(angle_vals)
    if angle_vals(i) >= intake_open && angle_vals(i) <= intake_close
        theta = angle_vals(i) - intake_open;
        intake_lift(i) = sin(pi * theta / intake_duration)^2;
    end
end

% Exhaust Valve: 110° → 360°
exhaust_lift = zeros(size(angle_vals));
exhaust_open = 110; exhaust_close = 360;
exhaust_duration = exhaust_close - exhaust_open;
for i = 1:length(angle_vals)
    if angle_vals(i) >= exhaust_open && angle_vals(i) <= exhaust_close
        theta = angle_vals(i) - exhaust_open;
        exhaust_lift(i) = sin(pi * theta / exhaust_duration)^2;
    end
end

% Normalize to mm (Max Lift = 10 mm)
max_lift = 0.01;
intake_lift = intake_lift * max_lift * 1000;
exhaust_lift = exhaust_lift * max_lift * 1000;

% === Plot: Valve Lift vs Crank Angle ===
figure;
plot(angle_vals, intake_lift, 'b', 'LineWidth', 2); hold on;
plot(angle_vals, exhaust_lift, 'r--', 'LineWidth', 2);
xlabel('Crank Angle [°]');
ylabel('Valve Lift [mm]');
title('Valve Lift Profiles of Intake and Exhaust Valves');
legend('Intake Valve','Exhaust Valve');
grid on;

