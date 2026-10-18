clear all; clc; close all;

%% 1. Definire Geometrie NACA 23015 (Camber Line pe intervale)
syms x
% Functiile specifice seriei NACA 230xx (conform Abbott & von Doenhoff)
f1 = 2.6595*(x^3 - 0.6075*x^2 + 0.1147*x); 
f2 = 0.02208*(1 - x);                     

funcs = {f1, f2};
intervals = [0 0.2025; 0.2025 1.0]; % Punctul de jonctiune standard pentru p=0.15

%% 2. Calcul Teoretic - Teoria Profilelor Subtiri
alfa = -8:0.5:14; % Intervalul cerut in proiect
m = length(alfa);
CL_teo = zeros(1, m);
CM_teo = zeros(1, m);

for i = 1:m
    % Apelam functia principala pentru fiecare unghi
    rez = thinAirfoilTheory(alfa(i), funcs, intervals);
    CL_teo(i) = rez.CL;
    CM_teo(i) = rez.CM_c4;
end

%% 3. Incarcare Date Experimentale (Digitalizate de tine)
% Programul cauta fisierele .txt in folderul curent
if exist('cl_exp.txt','file') == 2
    d = load('cl_exp.txt');
    alfa_cl = d(:,1); CL_exp = d(:,2);
else
    alfa_cl = []; CL_exp = [];
    fprintf('ATENTIE: Nu am gasit fisierul cl_exp.txt!\n');
end

if exist('cm_exp.txt','file') == 2
    d2 = load('cm_exp.txt');
    alfa_cm = d2(:,1); CM_exp = d2(:,2);
else
    alfa_cm = []; CM_exp = [];
    fprintf('ATENTIE: Nu am gasit fisierul cm_exp.txt!\n');
end

%% 4. Reprezentare Grafica (Rezultate Numerice suprapuse peste Experiment)
figure('Color', 'w', 'Name', 'Proiect Teoria Profilelor Subtiri - NACA 23015');

% --- Graficul Coeficientului de Portanta (CL) ---
subplot(2,1,1);
plot(alfa, CL_teo, 'b-', 'LineWidth', 2); hold on;
if ~isempty(CL_exp)
    plot(alfa_cl, CL_exp, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
end
grid on; 
xlabel('Unghi de atac \alpha [deg]'); ylabel('C_L');
title('NACA 23015: Coeficientul de Portanta C_L vs \alpha');
legend('Teorie Profile Subtiri (Numeric)', 'Date Experimentale (Abbott)', 'Location', 'NorthWest');

% --- Graficul Coeficientului de Moment (Cm) ---
subplot(2,1,2);
plot(alfa, CM_teo, 'b-', 'LineWidth', 2); hold on;
if ~isempty(CM_exp)
    plot(alfa_cm, CM_exp, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
end
grid on; 
xlabel('Unghi de atac \alpha [deg]'); ylabel('C_{m,c/4}');
title('NACA 23015: Coeficientul de Moment la sfertul coardei C_{m,c/4} vs \alpha');
ylim([-0.1 0.1]); % Zoom pe zona relevanta pentru moment
legend('Teorie Profile Subtiri (Numeric)', 'Date Experimentale', 'Location', 'SouthEast');

%% ========== FUNCTII AJUTATOARE (Nu modifica aici) ==========

function results = thinAirfoilTheory(alphaDeg, camberInput, intervals, N)
    if nargin < 3, intervals = []; end
    if nargin < 4, N = 5; end
    
    syms x theta real
    alpha = deg2rad(alphaDeg);
    
    % Construire functie zc
    if isa(camberInput, 'sym')
        zc = camberInput;
    elseif iscell(camberInput)
        zc = buildPiecewiseFromCells(camberInput, intervals);
    end
    
    % Derivata dz/dx si mapare pe cercul unitate (theta)
    dzdx = diff(zc, x);
    xTheta = (1 - cos(theta))/2;
    dzdx_theta = simplify(subs(dzdx, x, xTheta));
    
    % Calcul Coeficienti Fourier
    A0 = alpha - (1/pi) * double(vpaintegral(dzdx_theta, theta, 0, pi));
    A = zeros(1, N);
    for n = 1:N
        integrand = dzdx_theta * cos(n*theta);
        A(n) = (2/pi) * double(vpaintegral(integrand, theta, 0, pi));
    end
    
    A1 = A(1); 
    A2 = (N >= 2) * A(2);
    
    % Rezultate finale
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