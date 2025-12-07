% =========================================================================
% Flyback Converter Simulation - State Space Average Model
% =========================================================================
% Author: LAOUAR Ouassim
% Institution: CentraleSupélec
% 
% Based on the paper:
% A. S. Raj, A. M. Siddeshwar, Guruswamy K P, Maheshan C M and 
% Vijay Sanekere C, "Modelling of flyback converter using state space 
% averaging technique," 2015 IEEE International Conference on Electronics, 
% Computing and Communication Technologies (CONECCT), Bangalore, India, 
% 2015, pp. 1-5, doi: 10.1109/CONECCT.2015.7383871.
%
% Description:
% This script simulates a flyback converter using state-space averaging
% technique. The converter model is split into two sub-intervals based on
% switch state (ON/OFF) and solved using RK4 numerical integration.
%
% State Variables: [i_Lm; v_C]
% - i_Lm: Magnetizing inductance current
% - v_C: Output capacitor voltage
% =========================================================================

clc;
clear;
close all;

%% ========================================================================
%  1. SIMULATION PARAMETERS
% =========================================================================

% Control parameters
Duty_Cycle = 0.8;      % PWM duty cycle (0.0 to 1.0)
f_sw       = 100e3;    % Switching frequency [Hz]
T_sim      = 0.01;     % Total simulation time [s]

% Circuit parameters
Vin  = 24;             % Input voltage [V]
Vd   = 0.7;            % Diode forward voltage drop [V]
R    = 10;             % Load resistance [Ω]
C    = 47e-6;          % Output capacitance [F]
Lm   = 100e-6;         % Magnetizing inductance [H]
Rsw  = 0.1;            % Switch on-resistance [Ω]
Rc   = 0.05;           % Capacitor ESR [Ω]
n    = 2;              % Transformer turns ratio (Np/Ns)

%% ========================================================================
%  2. STATE-SPACE MATRICES
% =========================================================================
% Two operating modes based on equations from the referenced paper:
% - Mode 1 (Switch ON): Primary side energy storage
% - Mode 2 (Switch OFF): Energy transfer to secondary side
% =========================================================================

% Common denominator terms
denom_time_const = C * (R + Rc);
denom_time_const2 = C * (R + Rc);
denom_L_damping  = Lm * (R + Rc);

% --- Mode 1: Switch ON ---
% Magnetizing inductance charges from input
A1 = [ -Rsw/Lm,              0;
        0,                  -1/denom_time_const ];

B1 = [  1/Lm,                0;
        0,                   0 ];
   
C1 = [  0,                   R/(R+Rc) ];
D1 = [  0,                   0 ];

% --- Mode 2: Switch OFF ---
% Energy transfers to secondary through transformer
A2 = [ -(n^2 * Rc * R)/denom_L_damping,    -(n * R)/denom_L_damping;
        (n * R)/denom_time_const2,         -1/denom_time_const ];

B2 = [  0,                  -n/Lm;
        0,                   0 ];

C2 = [ -(n*R*Rc)/(R+Rc),     R/(R+Rc) ];  
D2 = [  0,                   0 ];

%% ========================================================================
%  3. TIME DISCRETIZATION
% =========================================================================

% Use high resolution for accurate PWM representation
points_per_cycle = 200;
dt = 1/(f_sw * points_per_cycle);
time = 0:dt:T_sim;
N = length(time);

fprintf('=== Simulation Configuration ===\n');
fprintf('Time step: %.2e s\n', dt);
fprintf('Total points: %d\n', N);
fprintf('Points per switching cycle: %d\n', points_per_cycle);
fprintf('Switching period: %.2e s\n\n', 1/f_sw);

%% ========================================================================
%  4. NUMERICAL INTEGRATION (RK4 METHOD)
% =========================================================================

% Initialize state vector [i_Lm; v_C]
x = zeros(2, N); 
x(:,1) = [0; 0];  % Initial conditions: zero current and voltage

% Output and mode tracking arrays
y_out = zeros(1, N);  % Output voltage
mode = zeros(1, N);   % Switch state (1=ON, 0=OFF)

% Input vector
u = [Vin; Vd]; 
T_period = 1/f_sw;

% Main simulation loop
fprintf('Running simulation');

for k = 1:N-1
    t_curr = time(k);
    
    % Progress indicator
    if mod(k, floor(N/10)) == 0
        fprintf('.');
    end
    
    % PWM control logic
    t_in_cycle = mod(t_curr, T_period);
    
    if t_in_cycle < (Duty_Cycle * T_period)
        % Switch ON: Use Mode 1 matrices
        A = A1; 
        B = B1; 
        C_out = C1; 
        D_out = D1;
        mode(k) = 1;
    else
        % Switch OFF: Use Mode 2 matrices
        A = A2; 
        B = B2; 
        C_out = C2; 
        D_out = D2;
        mode(k) = 0;
    end
    
    % RK4 integration for state update
    k1 = A * x(:,k) + B * u;
    k2 = A * (x(:,k) + 0.5*dt*k1) + B * u;
    k3 = A * (x(:,k) + 0.5*dt*k2) + B * u;
    k4 = A * (x(:,k) + dt*k3) + B * u;
    
    x(:,k+1) = x(:,k) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % Calculate output voltage
    y_out(k) = C_out * x(:,k) + D_out * u;
end

% Final point
y_out(N) = C2 * x(:,N) + D2 * u;
mode(N) = mode(N-1);

fprintf(' Complete!\n\n');

%% ========================================================================
%  5. PERFORMANCE METRICS (STEADY-STATE ANALYSIS)
% =========================================================================

% Analyze last 20% of simulation (steady-state region)
steady_idx = floor(0.8*N):N;

% Voltage metrics
V_out_avg = mean(y_out(steady_idx));
V_out_ripple = max(y_out(steady_idx)) - min(y_out(steady_idx));
V_cap_avg = mean(x(2,steady_idx));

% Current metrics
I_Lm_avg = mean(x(1,steady_idx));
I_Lm_peak = max(x(1,steady_idx));
I_Lm_ripple = max(x(1,steady_idx)) - min(x(1,steady_idx));

% Power and efficiency
P_out = V_out_avg^2/R;
P_in = Vin*I_Lm_avg*Duty_Cycle;
efficiency = 100*P_out/P_in;

% Ideal flyback output voltage (for comparison)
V_out_ideal = ((1/n) * Duty_Cycle * Vin) / (1 - Duty_Cycle);

% Display results
fprintf('=== Performance Metrics (Steady-State) ===\n');
fprintf('Output Voltage:\n');
fprintf('  Average: %.3f V\n', V_out_avg);
fprintf('  Ideal (theory): %.3f V\n', V_out_ideal);
fprintf('  Error: %.2f%%\n', 100*abs(V_out_avg-V_out_ideal)/V_out_ideal);
fprintf('  Ripple: %.3f V (%.2f%%)\n\n', V_out_ripple, 100*V_out_ripple/V_out_avg);

fprintf('Magnetizing Current:\n');
fprintf('  Average: %.3f A\n', I_Lm_avg);
fprintf('  Peak: %.3f A\n', I_Lm_peak);
fprintf('  Ripple: %.3f A\n\n', I_Lm_ripple);

fprintf('Power:\n');
fprintf('  Output: %.3f W\n', P_out);
fprintf('  Efficiency (approx): %.1f%%\n\n', efficiency);

%% ========================================================================
%  5b. STATE-SPACE AVERAGING & SMALL-SIGNAL ANALYSIS
% =========================================================================

% --- Large-Signal Averaged Model ---
fprintf('=== Computing State-Space Averaged Model ===\n');

d  = Duty_Cycle;
d1 = d;        % First sub-interval (switch ON)
d2 = 1 - d1;   % Second sub-interval (switch OFF)

% Averaged continuous-time state-space matrices
A_bar = d1*A1 + d2*A2;
B_bar = d1*B1 + d2*B2;
C_bar = d1*C1 + d2*C2;
D_bar = d1*D1 + d2*D2;

% Build LTI averaged model with inputs [Vin; Vd]
sys_avg = ss(A_bar, B_bar, C_bar, D_bar);

% Transfer functions from inputs to output (large-signal)
G_vin = tf(sys_avg(:,1));   % v_out / Vin
G_vd  = tf(sys_avg(:,2));   % v_out / Vd

fprintf('Large-signal averaged model created.\n\n');

% --- Small-Signal Model ---
fprintf('=== Deriving Small-Signal Transfer Functions ===\n');

% Operating point (use steady-state averages from simulation)
x_op = [I_Lm_avg; V_cap_avg];   % [i_Lm_op; v_C_op]
u_op = [Vin; Vd];               % [Vin_op; Vd_op]

% Sensitivity of state and output to duty perturbation d_hat
E_d = (A1 - A2)*x_op + (B1 - B2)*u_op;  % State sensitivity (2x1)
F_d = (C1 - C2)*x_op + (D1 - D2)*u_op;  % Output sensitivity (scalar)

% Build small-signal LTI model: Gvd(s) = v_out_hat(s) / d_hat(s)
A_lin = A_bar;
B_lin = E_d;       % Duty input
C_lin = C_bar;
D_lin = F_d;
sys_gvd = ss(A_lin, B_lin, C_lin, D_lin);

% Control-to-output transfer function
Gvd = tf(sys_gvd);

fprintf('Control-to-Output Transfer Function Gvd(s) = v_out_hat(s) / d_hat(s):\n');
disp(Gvd);

% Small-signal input-to-output: Gvg(s) = v_out_hat(s) / Vin_hat(s)
B_vin = B_bar(:,1);     % Column corresponding to Vin
D_vin = D_bar(1,1);     % Direct term Vin -> v_out
sys_gvg = ss(A_bar, B_vin, C_bar, D_vin);
Gvg = tf(sys_gvg);

fprintf('\nInput-to-Output Transfer Function Gvg(s) = v_out_hat(s) / Vin_hat(s):\n');
disp(Gvg);
fprintf('\n');

%% ========================================================================
%  6. VISUALIZATION
% =========================================================================

% Main analysis figure
figure('Color','w', 'Position', [50 50 1600 900]);

% Define zoom window (last 5 switching cycles)
t_zoom_start = T_sim - 5/f_sw;
zoom_idx = time >= t_zoom_start;

% Color scheme
color_blue = [0 0.4470 0.7410];
color_orange = [0.8500 0.3250 0.0980];
color_purple = [0.4940 0.1840 0.5560];
color_green = [0.4660 0.6740 0.1880];
color_red = [0.6350 0.0780 0.1840];

% --- Row 1: Full time-domain waveforms ---

subplot(3,3,1);
plot(time*1000, x(1,:), 'LineWidth', 1.5, 'Color', color_blue);
grid on;
title('Magnetizing Current i_{Lm}');
ylabel('Current [A]');
xlabel('Time [ms]');
xlim([0 T_sim*1000]);
text(0.02, 0.98, sprintf('Avg: %.3f A\nPeak: %.3f A\nRipple: %.3f A', ...
    I_Lm_avg, I_Lm_peak, I_Lm_ripple), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 8);

subplot(3,3,2);
plot(time*1000, y_out, 'LineWidth', 1.5, 'Color', color_orange);
hold on;
yline(V_out_avg, '--k', 'LineWidth', 1.5);
yline(V_out_ideal, '-.', 'Color', color_green, 'LineWidth', 1.5);
grid on;
title('Output Voltage v_{out}');
ylabel('Voltage [V]');
xlabel('Time [ms]');
xlim([0 T_sim*1000]);
legend('v_{out}', sprintf('Avg: %.2fV', V_out_avg), sprintf('Ideal: %.2fV', V_out_ideal), ...
    'Location', 'southeast', 'FontSize', 7);

subplot(3,3,3);
plot(time*1000, x(2,:), 'LineWidth', 1.5, 'Color', color_purple);
hold on;
yline(V_cap_avg, '--k', 'LineWidth', 1.5);
grid on;
title('Capacitor Voltage v_C');
ylabel('Voltage [V]');
xlabel('Time [ms]');
xlim([0 T_sim*1000]);
text(0.98, 0.98, sprintf('Avg: %.3f V', V_cap_avg), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 8);

% --- Row 2: Zoomed view (last 5 cycles) ---

subplot(3,3,4);
plot(time(zoom_idx)*1000, x(1,zoom_idx), 'LineWidth', 2, 'Color', color_blue);
hold on;
yline(I_Lm_avg, '--r', 'LineWidth', 1.5);
grid on;
title('i_{Lm} - Zoomed (Last 5 Cycles)');
ylabel('Current [A]');
xlabel('Time [ms]');

subplot(3,3,5);
plot(time(zoom_idx)*1000, y_out(zoom_idx), 'LineWidth', 2, 'Color', color_orange);
hold on;
yline(V_out_avg, '--k', 'LineWidth', 1.5);
yline(V_out_avg + V_out_ripple/2, ':', 'Color', color_red, 'LineWidth', 1);
yline(V_out_avg - V_out_ripple/2, ':', 'Color', color_red, 'LineWidth', 1);
grid on;
title('v_{out} - Zoomed (Last 5 Cycles)');
ylabel('Voltage [V]');
xlabel('Time [ms]');
legend('v_{out}', 'Average', 'Ripple bounds', 'Location', 'best', 'FontSize', 7);

subplot(3,3,6);
plot(time(zoom_idx)*1000, mode(zoom_idx), 'LineWidth', 2, 'Color', color_purple);
grid on;
title('Switch State (PWM) - Zoomed');
ylabel('State');
xlabel('Time [ms]');
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'OFF', 'ON'});

% --- Row 3: Analysis plots ---

subplot(3,3,7);
plot(x(1,steady_idx), x(2,steady_idx), 'LineWidth', 1, 'Color', color_blue);
hold on;
plot(x(1,steady_idx(1)), x(2,steady_idx(1)), 'go', 'MarkerSize', 8, 'LineWidth', 2);
plot(x(1,steady_idx(end)), x(2,steady_idx(end)), 'rs', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
title('Phase Plane (Steady-State)');
xlabel('i_{Lm} [A]');
ylabel('v_C [V]');
legend('Trajectory', 'Start', 'End', 'Location', 'best', 'FontSize', 7);
axis equal;

subplot(3,3,8);
P_inst = y_out.^2 / R;
plot(time*1000, P_inst, 'LineWidth', 1.5, 'Color', color_green);
hold on;
yline(P_out, '--k', 'LineWidth', 1.5);
grid on;
title('Instantaneous Output Power');
ylabel('Power [W]');
xlabel('Time [ms]');
xlim([0 T_sim*1000]);
legend('P_{out}(t)', sprintf('Avg: %.2fW', P_out), 'Location', 'southeast', 'FontSize', 7);

% Efficiency over time
subplot(3,3,9);
eta_inst = (y_out.^2 / R) ./ (Vin * x(1,:) * Duty_Cycle + eps) * 100;
eta_inst(eta_inst > 100) = 100; % Cap at 100%
eta_inst(eta_inst < 0) = 0;     % No negative efficiency
plot(time*1000, eta_inst, 'LineWidth', 1.5, 'Color', color_red);
hold on;
yline(efficiency, '--k', 'LineWidth', 1.5);
grid on;
title('Instantaneous Efficiency');
ylabel('Efficiency [%]');
xlabel('Time [ms]');
xlim([0 T_sim*1000]);
ylim([0 min(110, max(eta_inst)*1.1)]);
legend('η(t)', sprintf('Avg: %.1f%%', efficiency), 'Location', 'southeast', 'FontSize', 7);

% Main title
sgtitle(sprintf('Flyback Converter Analysis | V_{in}=%.1fV, D=%.2f, f_{sw}=%.0fkHz', ...
    Vin, Duty_Cycle, f_sw/1e3), 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%  7. PERFORMANCE SUMMARY FIGURE
% =========================================================================

figure('Color','w', 'Position', [100 100 800 600]);
axis off;

% Title
text(0.5, 0.95, '\bf\fontsize{16}Flyback Converter - Performance Summary', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top');

% Circuit schematic info
text(0.5, 0.87, sprintf('State-Space Averaging Model | Based on IEEE Paper doi:10.1109/CONECCT.2015.7383871'), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);

% Left column - Input parameters
text(0.05, 0.78, '\bf\fontsize{12}INPUT PARAMETERS', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'Color', [0 0.4470 0.7410]);
line([0.05 0.45], [0.76 0.76], 'Color', [0 0.4470 0.7410], 'LineWidth', 2);

input_text = {
    sprintf('\\bf Input Voltage (V_{in}):  \\rm%.1f V', Vin);
    sprintf('\\bf Duty Cycle (D):  \\rm%.3f (%.1f%%)', Duty_Cycle, Duty_Cycle*100);
    sprintf('\\bf Switching Frequency (f_{sw}):  \\rm%.0f kHz', f_sw/1e3);
    sprintf('\\bf Load Resistance (R):  \\rm%.1f Ω', R);
    sprintf('\\bf Output Capacitance (C):  \\rm%.1f µF', C*1e6);
    sprintf('\\bf Magnetizing Inductance (L_m):  \\rm%.0f µH', Lm*1e6);
    sprintf('\\bf Turns Ratio (n = N_p/N_s):  \\rm%.2f', n);
    sprintf('\\bf Switch Resistance (R_{sw}):  \\rm%.2f Ω', Rsw);
    sprintf('\\bf Capacitor ESR (R_c):  \\rm%.3f Ω', Rc);
    sprintf('\\bf Diode Voltage Drop (V_d):  \\rm%.2f V', Vd);
};
text(0.08, 0.72, input_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'tex');

% Right column - Output performance
text(0.55, 0.78, '\bf\fontsize{12}OUTPUT PERFORMANCE', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'Color', [0.8500 0.3250 0.0980]);
line([0.55 0.95], [0.76 0.76], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

output_text = {
    sprintf('\\bf Average Output Voltage:  \\rm%.3f V', V_out_avg);
    sprintf('\\bf Ideal Output Voltage:  \\rm%.3f V', V_out_ideal);
    sprintf('\\bf Voltage Error:  \\rm%.2f %%', 100*abs(V_out_avg-V_out_ideal)/V_out_ideal);
    sprintf('\\bf Output Voltage Ripple:  \\rm%.0f mV (%.2f%%)', V_out_ripple*1000, 100*V_out_ripple/V_out_avg);
    sprintf('\\bf Output Power:  \\rm%.3f W', P_out);
    sprintf('\\bf Input Power (approx):  \\rm%.3f W', P_in);
    sprintf('\\bf Efficiency:  \\rm%.1f %%', efficiency);
};
text(0.58, 0.72, output_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'tex');

% Bottom section - Inductor current analysis
text(0.05, 0.38, '\bf\fontsize{12}MAGNETIZING CURRENT ANALYSIS', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'Color', [0.4940 0.1840 0.5560]);
line([0.05 0.95], [0.36 0.36], 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);

current_text = {
    sprintf('\\bf Average Current (I_{Lm,avg}):  \\rm%.3f A', I_Lm_avg);
    sprintf('\\bf Peak Current (I_{Lm,peak}):  \\rm%.3f A', I_Lm_peak);
    sprintf('\\bf Current Ripple (ΔI_{Lm}):  \\rm%.0f mA (%.2f%%)', I_Lm_ripple*1000, 100*I_Lm_ripple/I_Lm_avg);
    sprintf('\\bf RMS Current (approx):  \\rm%.3f A', sqrt(I_Lm_avg^2 + I_Lm_ripple^2/12));
};
text(0.08, 0.32, current_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'tex');

% Capacitor voltage analysis
text(0.55, 0.38, '\bf\fontsize{12}CAPACITOR VOLTAGE', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'Color', [0.4660 0.6740 0.1880]);
line([0.55 0.95], [0.36 0.36], 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);

cap_text = {
    sprintf('\\bf Average Capacitor Voltage:  \\rm%.3f V', V_cap_avg);
    sprintf('\\bf Capacitor Voltage Ripple:  \\rm%.0f mV', (max(x(2,steady_idx))-min(x(2,steady_idx)))*1000);
};
text(0.58, 0.32, cap_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'tex');

% Footer with simulation info
text(0.5, 0.08, sprintf('Simulation: %.0f ms | Time step: %.2e s | Method: RK4', ...
    T_sim*1000, dt), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', [0.5 0.5 0.5]);

text(0.5, 0.03, sprintf('Author: LAOUAR Ouassim | CentraleSupélec'), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 8, 'Color', [0.3 0.3 0.3], 'FontWeight', 'bold');

%% ========================================================================
%  8. FREQUENCY-DOMAIN ANALYSIS
% =========================================================================

% Extract poles and zeros for analysis
poles_gvd = pole(Gvd);
zeros_gvd = zero(Gvd);
poles_gvg = pole(Gvg);
zeros_gvg = zero(Gvg);

% Calculate stability margins
[Gm_gvd, Pm_gvd, Wgm_gvd, Wpm_gvd] = margin(Gvd);
[Gm_gvg, Pm_gvg, Wgm_gvg, Wpm_gvg] = margin(Gvg);

fprintf('=== Stability Analysis ===\n');
fprintf('Control-to-Output Gvd(s):\n');
fprintf('  Gain Margin: %.2f dB at %.2f rad/s\n', 20*log10(Gm_gvd), Wgm_gvd);
fprintf('  Phase Margin: %.2f° at %.2f rad/s\n', Pm_gvd, Wpm_gvd);
fprintf('  Poles: '); disp(poles_gvd');
fprintf('  Zeros: '); disp(zeros_gvd');
fprintf('\nInput-to-Output Gvg(s):\n');
fprintf('  Gain Margin: %.2f dB at %.2f rad/s\n', 20*log10(Gm_gvg), Wgm_gvg);
fprintf('  Phase Margin: %.2f° at %.2f rad/s\n\n', Pm_gvg, Wpm_gvg);

% --- Figure 1: Bode Plots (Enhanced) ---
figure('Color','w', 'Position', [50 50 1400 800]);

% Gvd Bode Plot
subplot(2,2,1);
[mag_gvd, phase_gvd, wout_gvd] = bode(Gvd);
mag_gvd = squeeze(mag_gvd);
phase_gvd = squeeze(phase_gvd);
wout_gvd = squeeze(wout_gvd);

yyaxis left
semilogx(wout_gvd, 20*log10(mag_gvd), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
ylabel('Magnitude [dB]', 'FontSize', 10);
grid on;
hold on;
% Mark crossover frequency
[~, idx_0dB] = min(abs(20*log10(mag_gvd)));
plot(wout_gvd(idx_0dB), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

yyaxis right
semilogx(wout_gvd, phase_gvd, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
ylabel('Phase [deg]', 'FontSize', 10);
yline(-180, '--k', 'LineWidth', 1);
% Mark phase margin
if ~isinf(Wpm_gvd)
    plot(Wpm_gvd, Pm_gvd-180, 'rs', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('Frequency [rad/s]', 'FontSize', 10);
title('G_{vd}(s) - Control-to-Output', 'FontSize', 11, 'Interpreter', 'tex');
legend('Magnitude', '0 dB', 'Phase', '-180°', 'PM', 'Location', 'southwest', 'FontSize', 8);

% Gvg Bode Plot
subplot(2,2,2);
[mag_gvg, phase_gvg, wout_gvg] = bode(Gvg);
mag_gvg = squeeze(mag_gvg);
phase_gvg = squeeze(phase_gvg);
wout_gvg = squeeze(wout_gvg);

yyaxis left
semilogx(wout_gvg, 20*log10(mag_gvg), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
ylabel('Magnitude [dB]', 'FontSize', 10);
grid on;
hold on;
[~, idx_0dB_gvg] = min(abs(20*log10(mag_gvg)));
plot(wout_gvg(idx_0dB_gvg), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

yyaxis right
semilogx(wout_gvg, phase_gvg, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
ylabel('Phase [deg]', 'FontSize', 10);
yline(-180, '--k', 'LineWidth', 1);
if ~isinf(Wpm_gvg)
    plot(Wpm_gvg, Pm_gvg-180, 'rs', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('Frequency [rad/s]', 'FontSize', 10);
title('G_{vg}(s) - Input-to-Output', 'FontSize', 11, 'Interpreter', 'tex');
legend('Magnitude', '0 dB', 'Phase', '-180°', 'PM', 'Location', 'southwest', 'FontSize', 8);

% Nyquist Plot - Gvd
subplot(2,2,3);
nyquist(Gvd);
grid on;
hold on;
% Draw unit circle
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), '--k', 'LineWidth', 1);
% Mark critical point
plot(-1, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 3);
title('Nyquist Diagram - G_{vd}(s)', 'FontSize', 11, 'Interpreter', 'tex');
legend('Nyquist', 'Unit Circle', 'Critical Point (-1,0)', 'Location', 'best', 'FontSize', 8);
xlim([-2 2]);
ylim([-2 2]);

% Nyquist Plot - Gvg
subplot(2,2,4);
nyquist(Gvg);
grid on;
hold on;
plot(cos(theta), sin(theta), '--k', 'LineWidth', 1);
plot(-1, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 3);
title('Nyquist Diagram - G_{vg}(s)', 'FontSize', 11, 'Interpreter', 'tex');
legend('Nyquist', 'Unit Circle', 'Critical Point (-1,0)', 'Location', 'best', 'FontSize', 8);
xlim([-2 2]);
ylim([-2 2]);

sgtitle('Frequency-Domain Analysis - Small-Signal Transfer Functions', ...
    'FontSize', 14, 'FontWeight', 'bold');


%% ========================================================================
%  9. OPTIONAL: EXPORT SIMULATION DATA
% =========================================================================
% Uncomment the following lines to save simulation results to .mat file
%
% results.time = time;
% results.current = x(1,:);
% results.voltage_cap = x(2,:);
% results.voltage_out = y_out;
% results.mode = mode;
% results.params = struct('Vin',Vin,'Duty',Duty_Cycle,'fsw',f_sw,...
%                         'R',R,'C',C,'Lm',Lm,'n',n);
% save('flyback_simulation_results.mat', 'results');
% fprintf('Results saved to flyback_simulation_results.mat\n');