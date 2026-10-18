clear all; clc; close all;

%% ============================================================
%  PROIECT COMPLET - NACA 23015
%  Metoda Panourilor Hess-Smith + Teoria Profilelor subtiri(laborator3)
%  Comparatie cu Datele Experimentale
%% ============================================================

%% 1. Parametri Generali si Interval Unghiuri
c = 1;          % Coarda profilului
Vinf = 87.6;    % Viteza curgerii la infinit

% Unghiurile de atac pentru simulari
alpha_deg_hs = -8:2:18; 
alpha_rad_hs = deg2rad(alpha_deg_hs);

alpha_deg_tat = -8:0.5:18; % Interval extins pentru continuitate pe grafic

%% 2. Geometrie Profil NACA 23015 (UIUC Database)
coordonate_uiuc = [
    1.0000     0.0016; 0.9500     0.0112; 0.9000     0.0204; 0.8000     0.0373;
    0.7000     0.0525; 0.6000     0.0661; 0.5000     0.0774; 0.4000     0.0859;
    0.3000     0.0905; 0.2500     0.0908; 0.2000     0.0892; 0.1500     0.0852;
    0.1000     0.0764; 0.0750     0.0690; 0.0500     0.0589; 0.0250     0.0444;
    0.0125     0.0334; 0.0000     0.0000; 0.0125    -0.0154; 0.0250    -0.0225;
    0.0500    -0.0304; 0.0750    -0.0361; 0.1000    -0.0409; 0.1500    -0.0484;
    0.2000    -0.0541; 0.2500    -0.0578; 0.3000    -0.0596; 0.4000    -0.0592;
    0.5000    -0.0550; 0.6000    -0.0481; 0.7000    -0.0391; 0.8000    -0.0283;
    0.9000    -0.0159; 0.9500    -0.0090; 1.0000    -0.0016
];
x_geom = coordonate_uiuc(:,1);
y_geom = coordonate_uiuc(:,2);

%% 3. Date Experimentale Integrate (lab4_cl.txt si lab4_cm.txt)
alfa_cl_exp = [-8.1633, -5.9864, -4.0272, -2.0680, 0.1088, 2.0680, 4.0272, 5.9864, 7.9456, 9.9048, 12.0816, 14.0408, 16.0000, 18.1769];
CL_exp      = [-0.6865, -0.5083, -0.3197, -0.0996, 0.0891, 0.3197, 0.5293, 0.7493, 0.9485, 1.1790,  1.3886,  1.5144,  1.6402,  1.5773];

alfa_cm_exp = [-8.1633, -6.2041, -4.0272, -2.0680, 0.1088, 2.0680, 4.0272, 5.9864, 7.9456, 9.9048, 12.0816, 14.0408, 16.0000, 18.1769];
CM_exp      = [-0.0891, -0.0681, -0.0576, -0.0472, -0.0262, -0.0157, -0.0052, -0.0052, -0.0052, -0.0052, -0.0157, -0.0157, -0.0157, -0.1415];

%% 4. Calcul Teoretic - Teoria Profilelor subtiri(laborator3)
syms x
f1 = 2.6595*(x^3 - 0.6075*x^2 + 0.1147*x); 
f2 = 0.02208*(1 - x);                     
funcs = {f1, f2};
intervals = [0 0.2025; 0.2025 1.0];

fprintf('Se ruleaza calculul pentru Teoria Profilelor Subtiri...\n');
CL_teo = zeros(size(alpha_deg_tat));
CM_teo = zeros(size(alpha_deg_tat));
for i = 1:length(alpha_deg_tat)
    rez = thinAirfoilTheory(alpha_deg_tat(i), funcs, intervals);
    CL_teo(i) = rez.CL;
    CM_teo(i) = rez.CM_c4;
end

%% 5. Calcul Numeric - Metoda Panourilor Hess-Smith
fprintf('Se ruleaza calculul pentru Metoda Panourilor Hess-Smith...\n');
CL_panel = zeros(size(alpha_rad_hs));
Cm_c4_panel = zeros(size(alpha_rad_hs));
Cp_all = cell(length(alpha_rad_hs), 1);
xcp_all = cell(length(alpha_rad_hs), 1);

for k = 1:length(alpha_rad_hs)
    alpha = -alpha_rad_hs(k); 
    sol = hessSmith_academic(x_geom, y_geom, Vinf, alpha);
    
    CL_panel(k) = calc_CL(sol, c);
    Cm_c4_panel(k) = calc_Cm_c4(sol, c, alpha);
    Cp_all{k} = sol.Cp;
    xcp_all{k} = sol.xc/c;
end

%% ============================================================
%  REPREZENTARI GRAFICE
%% ============================================================

%% Grafic 1: CL vs Alpha
figure('Color', 'w', 'Name', 'Comparatie CL - NACA 23015');
plot(alpha_deg_hs, CL_panel, 'd-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
hold on;
plot(alpha_deg_tat, CL_teo, 'b--', 'LineWidth', 2);
plot(alfa_cl_exp, CL_exp, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Unghi de atac \alpha [deg]'); ylabel('C_L');
title('C_L vs \alpha pentru profilul NACA 23015');
legend('Metoda Panourilor Hess-Smith', 'Teoria Profilelor subtiri(laborator3)', 'Date Experimentale (Theory of Wing Sections)', 'Location', 'NorthWest');

%% Grafic 2: Cm vs Alpha
figure('Color', 'w', 'Name', 'Comparatie Cm - NACA 23015');
plot(alpha_deg_hs, Cm_c4_panel, 'd-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
hold on;
plot(alpha_deg_tat, CM_teo, 'b--', 'LineWidth', 2);
plot(alfa_cm_exp, CM_exp, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('Unghi de atac \alpha [deg]'); ylabel('C_{m,c/4}');
title('C_{m,c/4} vs \alpha pentru profilul NACA 23015');
ylim([-0.2 0.1]);
legend('Metoda Panourilor Hess-Smith', 'Teoria Profilelor subtiri(laborator3)', 'Date Experimentale (Theory of Wing Sections)', 'Location', 'SouthWest');

%% Grafic 3: Distributia Coeficientului de Presiune Cp (Hess-Smith)
alpha_Cp_deg = [-4 0 4 8 12];
figure('Color', 'w', 'Name', 'Distributie Cp - Hess-Smith');
hold on;
stiluri = {'-','--','-.',':','-'};
culori = [0 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];

for k = 1:length(alpha_Cp_deg)
    idx = find(alpha_deg_hs == alpha_Cp_deg(k));
    if ~isempty(idx)
        plot(xcp_all{idx}, Cp_all{idx}, 'LineStyle', stiluri{k}, 'Color', culori(k,:), 'LineWidth', 2.0, ...
            'DisplayName', sprintf('\\alpha = %d^\\circ', alpha_Cp_deg(k)));
    end
end
set(gca, 'YDir', 'reverse');
grid on; box on;
xlabel('x/c'); ylabel('C_p');
title('Distributia coeficientului de presiune C_p (Hess-Smith) - NACA 23015');
legend('Location', 'Northeast');
xlim([-0.1 1.1]);

%% Grafic 4: Geometria Profilului
figure('Color', 'w', 'Name', 'Geometrie Profil');
plot(x_geom, y_geom, 'k-', 'LineWidth', 1.5);
axis equal; grid on; box on;
xlabel('x/c'); ylabel('y/c');
title('Profil NACA 23015 (Discretizare panouri)');


%% ============================================================
%% FUNCTII AUXILIARE
%% ============================================================

function results = thinAirfoilTheory(alphaDeg, camberInput, intervals, N)
    if nargin < 3, intervals = []; end
    if nargin < 4, N = 5; end
    
    syms x theta real
    alpha = deg2rad(alphaDeg);
    
    if isa(camberInput, 'sym')
        zc = camberInput;
    elseif iscell(camberInput)
        zc = buildPiecewiseFromCells(camberInput, intervals);
    end
    
    dzdx = diff(zc, x);
    xTheta = (1 - cos(theta))/2;
    dzdx_theta = simplify(subs(dzdx, x, xTheta));
    
    A0 = alpha - (1/pi) * double(vpaintegral(dzdx_theta, theta, 0, pi));
    A = zeros(1, N);
    for n = 1:N
        integrand = dzdx_theta * cos(n*theta);
        A(n) = (2/pi) * double(vpaintegral(integrand, theta, 0, pi));
    end
    
    A1 = A(1); 
    A2 = (N >= 2) * A(2);
    
    results = struct();
    results.CL = pi * (2*A0 + A1);
    results.CM_c4 = (pi/4) * (A2 - A1);
end

function zc = buildPiecewiseFromCells(funcs, intervals)
    syms x real
    n = numel(funcs);
    args = cell(1, 2*n);
    for i = 1:n
        xa = intervals(i, 1); xb = intervals(i, 2);
        cond = (x >= xa) & (x <= xb);
        args{2*i-1} = cond;
        args{2*i} = funcs{i};
    end
    zc = piecewise(args{:});
end

function CL = calc_CL(sol, c)
    CL = 0;
    for j = 1:length(sol.Cp)
        CL = CL + (-sol.Cp(j)) * sol.geom.ny(j) * sol.geom.S(j) / c;
    end
end

function Cm = calc_Cm_c4(sol, c, alpha)
    x_ref = 0.25*c; Mx = 0;
    for j = 1:length(sol.Cp)
        dx = sol.xc(j) - x_ref; dy = sol.yc(j);
        dFx_geom = -sol.Cp(j) * sol.geom.nx(j) * sol.geom.S(j);
        dFy_geom = -sol.Cp(j) * sol.geom.ny(j) * sol.geom.S(j);
        dFx = dFx_geom*cos(alpha) - dFy_geom*sin(alpha);
        dFy = dFy_geom*cos(alpha) + dFx_geom*sin(alpha);
        Mx = Mx + dFx*dy + dFy*dx;
    end
    Cm = -Mx/(c^2);
end

function sol = hessSmith_academic(x, y, Vinf, alpha)
    nPoints = numel(x); N = nPoints - 1;
    S = zeros(N,1); xc = zeros(N,1); yc = zeros(N,1);
    for i = 1:N
        S(i) = sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2);
        xc(i) = 0.5*(x(i) + x(i+1)); yc(i) = 0.5*(y(i) + y(i+1));
    end
    theta = zeros(N,1);
    for i = 1:N
        dx = x(i+1)-x(i); dy = y(i+1)-y(i);
        theta(i) = atan2(dy, dx);
    end
    tx = cos(theta); ty = sin(theta); nx = sin(theta); ny = -cos(theta);
    Vn_Source = zeros(N,N); Vt_Source = zeros(N,N); Vn_Vortex = zeros(N,N); Vt_Vortex = zeros(N,N);
    for i = 1:N
        for j = 1:N
            if i == j
                Vn_Source(i,j) = 0.5; Vt_Source(i,j) = 0.0; Vn_Vortex(i,j) = 0.0; Vt_Vortex(i,j) = 0.5;
            else
                [Ssrc, Svort] = panel_influence(xc(i), yc(i), x(j), y(j), theta(j), S(j));
                Vn_Source(i,j) = Ssrc.u*nx(i) + Ssrc.v*ny(i);
                Vt_Source(i,j) = Ssrc.u*tx(i) + Ssrc.v*ty(i);
                Vn_Vortex(i,j) = Svort.u*nx(i) + Svort.v*ny(i);
                Vt_Vortex(i,j) = Svort.u*tx(i) + Svort.v*ty(i);
            end
        end
    end
    Uinf = Vinf*cos(alpha); Vinf_y = -Vinf*sin(alpha);
    LHS = zeros(N+1, N+1); RHS = zeros(N+1, 1);
    for i = 1:N
        LHS(i,1:N) = Vn_Source(i,:); LHS(i,N+1) = sum(Vn_Vortex(i,:)); RHS(i) = -(Uinf*nx(i) + Vinf_y*ny(i));
    end
    iu = 1; il = N;
    LHS(N+1,1:N) = Vt_Source(iu,:) + Vt_Source(il,:); LHS(N+1,N+1) = sum(Vt_Vortex(iu,:)) + sum(Vt_Vortex(il,:));
    RHS(N+1) = -((Uinf*tx(iu) + Vinf_y*ty(iu)) + (Uinf*tx(il) + Vinf_y*ty(il)));
    unk = LHS \ RHS; sigma = unk(1:N); gamma = unk(N+1);
    Vt = zeros(N,1);
    for i = 1:N
        Vt(i) = (Uinf*tx(i) + Vinf_y*ty(i)) + sum(Vt_Source(i,:).*sigma(:)') + gamma * sum(Vt_Vortex(i,:));
    end
    Cp = 1 - (Vt./Vinf).^2;
    sol.xc = xc; sol.yc = yc; sol.sigma = sigma; sol.gamma = gamma; sol.Vt = Vt; sol.Cp = Cp; sol.iTE_upper = iu; sol.iTE_lower = il;
    sol.geom.x = x; sol.geom.y = y; sol.geom.S = S; sol.geom.phi = theta; sol.geom.tx = tx; sol.geom.ty = ty; sol.geom.nx = nx; sol.geom.ny = ny;
end

function [Ssrc, Svort] = panel_influence(xp, yp, xj, yj, thetaj, Sj)
    dx = xp - xj; dy = yp - yj;
    c = cos(thetaj); s = sin(thetaj);
    xloc = dx*c + dy*s; yloc = -dx*s + dy*c;
    r1sq = xloc^2 + yloc^2; r2sq = (xloc - Sj)^2 + yloc^2;
    th1 = atan2(yloc, xloc); th2 = atan2(yloc, xloc - Sj);
    lnTerm = 0.5*log(r1sq/r2sq); dth = th2 - th1; dth = atan2(sin(dth), cos(dth));
    uS = (1/(2*pi))*lnTerm; vS = (1/(2*pi))*dth; uV = -(1/(2*pi))*dth; vV = (1/(2*pi))*lnTerm;
    uSg = uS*c - vS*s; vSg = uS*s + vS*c; uVg = uV*c - vV*s; vVg = uV*s + vV*c;
    Ssrc = struct('u', uSg, 'v', vSg); Svort = struct('u', uVg, 'v', vVg);
end