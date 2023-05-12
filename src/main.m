% Matboard Bridge - small-scale box girder bridge
% Copyright © 2023 Nabeth Ghazi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

clear; close all;
% values in N, mm, and MPa
% using tension side as positive for SFD and BMD

%% 0. initialize parameters
global L P
L = 1200; % Length of bridge [mm]
P = 400; % total weight of train [N]
x = 0:L;

%% 1. SFD, BMD under train loading
x_train_max = 240;
V_env = zeros(x_train_max + 1, L + 1);
M_env = zeros(x_train_max + 1, L + 1);
for x_train = 0:x_train_max
    [x_load, SFD, BMD] = gen_diagrams(x_train);
    V_env(x_train + 1, :) = get(x_load, SFD, x);
    M_env(x_train + 1, :) = get(x_load, BMD, x);
end

[~, V_i] = max(abs(V_env));
[~, M_i] = max(abs(M_env));
V = zeros(1, L + 1);
M = zeros(1, L + 1);
for col = x + 1
    V(col) = V_env(V_i(col), col); % shear force envelope
    M(col) = M_env(M_i(col), col); % bending moment envelope
end

%% 2. define cross section parameters
x_params = [0];
params(1) = cross_section(  [75 105 1.27 1.27 10 10], ... % bi
                            [1.27 3.81 94.92 94.92 1.27 1.27], ... % hi
                            [0.635 98.095 48.73 48.73 95.555 95.555], ... % yi
                            100, ... % h_total
                            [96.19 97.46 98.73 1.27], ... % y_glue
                            ...
                            29.7987, ... % y_buck1
                            29.7987, ... % y_buck2
                            25.9887, ... % y_buck3
                            ...
                            3.81, ... % t_buck1
                            3.81, ... % t_buck2
                            1.27, ... % t_buck3
                            1.27, ... % t_buckV
                            ...
                            73.73, ... % b_buck1
                            15.635, ... % b_buck2
                            25.3537, ... % b_buck3
                            94.92, ... % b_buckV
                            ...
                            300); % a_buckV

%% 3. calculate capacity and applied stress
E = 4000;
mu = 0.2;

% applied stress
S_tens_app = zeros(size(x));
S_comp_app = zeros(size(x));
T_max_app = zeros(size(x));
T_glue_app = zeros(size(x));
S_buck1_app = zeros(size(x));
S_buck2_app = zeros(size(x));
S_buck3_app = zeros(size(x));
T_buck_app = zeros(size(x));

% capacity stress
S_tens_ult = 30;
S_comp_ult = -6;
T_max_ult = 4;
T_glue_ult = 2;
S_buck1_ult = zeros(size(x));
S_buck2_ult = zeros(size(x));
S_buck3_ult = zeros(size(x));
T_buck_ult = zeros(size(x));

for i = x + 1
    % cross section properties
    cs = get_piecewise(x_params, params, i);
    I = cs.I;
    ybot = cs.ybot;
    ytop = cs.ytop;
    ybar = cs.ybar;
    Q_max = cs.Q(ybar);
    b_centroid = cs.bV(ybar);

    % applied stress
    S_tens_app(i) = M(i) * ybot / I;
    S_comp_app(i) = M(i) * ytop / I;
    T_max_app(i) = abs(V(i) * Q_max / I / b_centroid);
    for y_glue = cs.y_glue
        Q_glue = cs.Q(y_glue);
        b_glue = cs.bV(y_glue);
        T_glue_app(i) = max(abs(V(i) * Q_glue / I / b_glue), T_glue_app(i));
    end
    S_buck1_app(i) = abs(M(i) * cs.y_buck1 / I);
    S_buck2_app(i) = abs(M(i) * cs.y_buck2 / I);
    S_buck3_app(i) = abs(M(i) * cs.y_buck3 / I);
    T_buck_app(i) = T_max_app(i);

    % local buckling stress capacities
    S_buck1_ult(i) = cross_section.S_buck_ult(4, E, mu, cs.t_buck1, cs.b_buck1);
    S_buck2_ult(i) = cross_section.S_buck_ult(0.425, E, mu, cs.t_buck2, cs.b_buck2);
    S_buck3_ult(i) = cross_section.S_buck_ult(6, E, mu, cs.t_buck3, cs.b_buck3);
    T_buck_ult(i) = cross_section.T_buck_ult(5, E, mu, cs.t_buckV, cs.b_buckV, ...
                                            cs.a_buckV);
end

%% 4. Factor Of Safety
FOS_tens = S_tens_ult ./ S_tens_app;
FOS_comp = S_comp_ult ./ S_comp_app;
FOS_shear = T_max_ult ./ T_max_app;
FOS_glue = T_glue_ult ./ T_glue_app;
FOS_buck1 = S_buck1_ult ./ S_buck1_app;
FOS_buck2 = S_buck2_ult ./ S_buck2_app;
FOS_buck3 = S_buck3_ult ./ S_buck3_app;
FOS_buckV = T_buck_ult ./ T_buck_app;

fprintf("FOS_tens = %g\n", slide_rule(min(FOS_tens)))
fprintf("FOS_comp = %g\n", slide_rule(min(FOS_comp)))
fprintf("FOS_shear = %g\n", slide_rule(min(FOS_shear)))
fprintf("FOS_glue = %g\n", slide_rule(min(FOS_glue)))
fprintf("FOS_buck1 = %g\n", slide_rule(min(FOS_buck1)))
fprintf("FOS_buck2 = %g\n", slide_rule(min(FOS_buck2)))
fprintf("FOS_buck3 = %g\n", slide_rule(min(FOS_buck3)))
fprintf("FOS_buckV = %g\n", slide_rule(min(FOS_buckV)))

%% 5. min FOS and the failure load Pfail
minFOS = min([FOS_tens FOS_comp FOS_shear FOS_glue ...
        FOS_buck1 FOS_buck2 FOS_buck3 FOS_buckV]);
Pf = P * minFOS;

fprintf("\nminFOS = %g\n", slide_rule(minFOS))
fprintf("P_fail = %g N\n", slide_rule(Pf))

%% 6. Vfail and Mfail
Mf_tens = FOS_tens .* M;
Mf_comp = FOS_comp .* M;
Vf_shear = FOS_shear .* V;
Vf_glue = FOS_glue .* V;
Mf_buck1 = FOS_buck1 .* M;
Mf_buck2 = FOS_buck2 .* M;
Mf_buck3 = FOS_buck3 .* M;
Vf_buckV = FOS_buckV .* V;

%% 7. output plots of Vfail and Mfail
figure
subplot (2, 3, 1)
hold on, grid on, grid minor
plot_shear(x, V_env, V, Vf_shear)
legend("Matboard Shear Failure")

subplot (2, 3, 2)
plot_shear(x, V_env, V, Vf_glue)
legend("Glue Shear Failure")
title ("Shear Force Diagram vs. Shear Force Capacities")

subplot (2, 3, 3)
plot_shear(x, V_env, V, Vf_buckV)
legend("Matboard Shear Buckling Failure")

subplot (2, 3, 4)
plot_moment(x, M_env, M, [Mf_comp; Mf_tens]')
legend("Matboard Compression Failure", "Matboard Tension Failure")

subplot (2, 3, 5)
plot_moment(x, M_env, M, [Mf_buck1; Mf_buck2]')
legend( "Matboard Flexural Buckling Case 1 Failure", ...
        "Matboard Flexural Buckling Case 2 Failure" )
title (["Bending Moment Diagram vs.", "Bending Moment Capacities"])

subplot (2, 3, 6)
plot_moment(x, M_env, M, Mf_buck3)
legend("Matboard Flexural Buckling Case 3 Failure")

%% functions

function plot_shear(x, V_env, V, Vf)
    hold on, grid on, grid minor

    plot(x, abs(Vf), "r")
    plot(x, -abs(Vf), "r")
    plot(x, V_env, "Color", "#c8c8c8")
    plot(x, V, "k")

    plot([0, x(end)], [0, 0], "k")
    xlabel("Distance Along bridge (mm)")
    ylabel("Shear Force (N)")
end


function plot_moment(x, M_env, M, Mf)
    hold on, grid on, grid minor

    plot(x, Mf)
    plot(x, M_env, "Color", "#c8c8c8")
    plot(x, M, "k")

    plot([0, x(end)], [0, 0], "k")
    xlabel("Distance Along bridge (mm)")
    ylabel("Bending Moment (Nmm)")
    set(gca, "YDir", "reverse")
end


function [x_load, SFD, BMD] = gen_diagrams(x_train)
    % declare parameters
    global L P

    % loading on truss
    P_train = ones(1, 6) * P / 6;
    x_load = [52 228 392 568 732 908] + x_train; % train Load positions

    rxn_B = sum(P_train .* x_load) / L;
    rxn_A = sum(P_train) - rxn_B;
    w = [rxn_A, -P_train, rxn_B];

    x_load = duplicate([0, x_load, 1200]);

    % tension side positive
    % calculate diagrams
    SFD = cumsum(w);
    SFD = SFD(1:end - 1);
    SFD = duplicate(SFD);
    SFD = [0, SFD, 0];

    BMD = cumtrapz(x_load, SFD);
end


function v = duplicate(v)
    v = [v; v];
    v = [v(:)]';
end


function res = get(x, v, xq)
    x(2:2:end) = x(2:2:end) + eps(x(2:2:end));
    res = interp1(x, v, xq);
end


function res = get_piecewise(x, v, xq)
    v = v(x <= xq);
    res = v(end);
end


function res = slide_rule(num)
    if abs(num) < 1e-9
        res = 0;
        return
    end

    sf = 3;
    if mod(log10(abs(num)), 1) < log10(2)
        sf = 4;
    end
    res = round(num, sf, "significant");
end