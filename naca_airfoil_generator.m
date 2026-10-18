clear; clc; close all;

%% ------------------------------------------------------------
% PARAMETRII GENERALI
%% ------------------------------------------------------------

Np = 300;
useCosineSpacing = true;

%% ------------------------------------------------------------
% SELECTIA METODEI
%% ------------------------------------------------------------

method = 'NACA5';

%% ------------------------------------------------------------
% PARAMETRII PROFILULUI (NACA 55010)
%% ------------------------------------------------------------

params.L = 5;        % din 55010
params.PQ = 50;      % din 55010
params.T = 10;       % grosime 10%
params.reflex = false;

%% ------------------------------------------------------------
% GENERARE PROFIL
%% ------------------------------------------------------------

[xu, yu, xl, yl, x, yc, yt] = generate_naca_profile(method, params, Np, useCosineSpacing);

%% ------------------------------------------------------------
% PLOT PROFIL
%% ------------------------------------------------------------

figure('Color','w');
hold on; grid on; axis equal;

plot(xu, yu, 'b-', 'LineWidth', 1.5);
plot(xl, yl, 'r-', 'LineWidth', 1.5);
plot(x, yc, 'k--', 'LineWidth', 1.0);
plot([0 1], [0 0], 'g-.', 'LineWidth', 1.0);

xlabel('x/c');
ylabel('y/c');
title('Profil generat: NACA 55010');
legend('Extrados','Intrados','Linie de curbura','Coarda','Location','Best');

%% ------------------------------------------------------------
% COMPARATIE CU NACA 0012
%% ------------------------------------------------------------

params_ref.M = 0;
params_ref.P = 0;
params_ref.T = 12;

[xu_ref, yu_ref, xl_ref, yl_ref, ~, ~, ~] = ...
    generate_naca_profile('NACA4', params_ref, Np, useCosineSpacing);

figure('Color','w');
hold on; grid on; axis equal;

plot(xu, yu, 'b-', 'LineWidth', 1.5);
plot(xl, yl, 'b-', 'LineWidth', 1.5);

plot(xu_ref, yu_ref, 'r--', 'LineWidth', 1.5);
plot(xl_ref, yl_ref, 'r--', 'LineWidth', 1.5);

plot([0 1], [0 0], 'k:', 'LineWidth', 1.0);

xlabel('x/c');
ylabel('y/c');
title('Comparatie intre NACA 55010 si NACA 0012');
legend('NACA 55010 extrados','NACA 55010 intrados', ...
       'NACA 0012 extrados','NACA 0012 intrados', ...
       'Coarda','Location','Best');

%% ------------------------------------------------------------
% EXPORT DATE
%% ------------------------------------------------------------

writematrix([xu(:), yu(:)], 'NACA55010_upper.dat', 'Delimiter', 'tab');
writematrix([xl(:), yl(:)], 'NACA55010_lower.dat', 'Delimiter', 'tab');

writematrix([xu_ref(:), yu_ref(:)], 'NACA0012_upper.dat', 'Delimiter', 'tab');
writematrix([xl_ref(:), yl_ref(:)], 'NACA0012_lower.dat', 'Delimiter', 'tab');

disp('Calcul terminat.');

%% ============================================================
% FUNCTII
%% ============================================================

function [xu, yu, xl, yl, x, yc, yt] = generate_naca_profile(method, params, Np, useCosineSpacing)

    x = chordwise_distribution(Np, useCosineSpacing);

    switch upper(method)
        case 'NACA4'
            [yc, dyc_dx] = camber_naca4(x, params.M, params.P);
            yt = thickness_naca4(x, params.T, true);

        case 'NACA5'
            [yc, dyc_dx] = camber_naca5(x, params.L, params.PQ, params.reflex);
            yt = thickness_naca4(x, params.T, true);
    end

    theta = atan(dyc_dx);

    xu = x - yt .* sin(theta);
    yu = yc + yt .* cos(theta);

    xl = x + yt .* sin(theta);
    yl = yc - yt .* cos(theta);

    xu(1)=0; xl(1)=0;
    xu(end)=1; xl(end)=1;
end

%% ------------------------------------------------------------

function x = chordwise_distribution(Np, useCosineSpacing)
    if useCosineSpacing
        beta = linspace(0, pi, Np);
        x = 0.5*(1-cos(beta));
    else
        x = linspace(0,1,Np);
    end
end

%% ------------------------------------------------------------

function [yc, dyc_dx] = camber_naca4(x, M, P)

    if M>1, M=M/100; end
    if P>1, P=P/10; end

    yc=zeros(size(x)); dyc_dx=zeros(size(x));

    if M==0 || P==0, return; end

    i1 = x<P;
    i2 = x>=P;

    yc(i1)=(M/P^2)*(2*P*x(i1)-x(i1).^2);
    dyc_dx(i1)=(2*M/P^2)*(P-x(i1));

    yc(i2)=(M/(1-P)^2)*((1-2*P)+2*P*x(i2)-x(i2).^2);
    dyc_dx(i2)=(2*M/(1-P)^2)*(P-x(i2));
end

%% ------------------------------------------------------------

function yt = thickness_naca4(x, T, sharpTE)

    if T>1, T=T/100; end

    if sharpTE
        a4=-0.1036;
    else
        a4=-0.1015;
    end

    yt = 5*T*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3+a4*x.^4);
    yt(1)=0; yt(end)=0;
end

%% ------------------------------------------------------------

function [yc, dyc_dx] = camber_naca5(x, L, PQ, reflex)

    if reflex
        error('Nu este implementat.');
    end

    cl_design = 0.15 * L;

    switch PQ
        case 50
            m = 0.3910;
            k1_ref = 3.230;
            cl_ref = 0.3;
        otherwise
            error('Doar PQ=50 folosit aici.');
    end

    k1 = k1_ref * (cl_design / cl_ref);

    yc=zeros(size(x)); dyc_dx=zeros(size(x));

    i1 = x<m;
    i2 = x>=m;
    yc(i1) = (k1 / 6) * (x(i1).^3 - 3*m*x(i1).^2 + m^2*(3 - m).*x(i1));
    dyc_dx(i1) = (k1/6)*(3.*x(i1).^2 - 6.*m.*x(i1) + (m^2)*(3-m));

    yc(i2) = (k1*m^3/6)*(1-x(i2));
    dyc_dx(i2) = -(k1*m^3/6);
end
