clear ;
clc;
close all;
% =========================
%     Octave Frequency
% =========================

nf=28;                             
f=ones(nf,1);
f(1)=20;                            
for i=2:nf
f(i)=f(i-1).*(2.^(1/3));             
end
f_signal = f;
f=f';
f_lr=f.*(2.^(-1/6));
f_up=f.*(2.^(1/6));
delta_f=f_up-f_lr;

% alpha = 0;                             
% % M =0.95;            % 0.95; 
% % T_infty = 253.725;
% % u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
% % rho_infty = 0.50335;              
% % nu = 0.00003214958;
% 
%  M = 0.95;            % 0.95; 
%  T_infty = 283.2197458;
%  u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
%  rho_infty = 0.84868;              
%  nu = 0.000020802;
%  q = 0.5*rho_infty*u_infty^2;
%  F = 1 + 0.13*M^2;


alpha = 0;                             
M = 1.6;          
T_infty =216.08;       
u_infty = M*sqrt(1.4*287*T_infty);                     
rho_infty = 0.3168841;                     
nu = 4.48E-05;      
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;


  % alpha = 0;                             
  % M = 3.2; 
  % T_infty = 212.3900404;
  % u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
  % rho_infty =0.2944102326;              
  % nu = 0.00004748820111;
  % q = 0.5*rho_infty*u_infty^2;
  % F = 1 + 0.13*M^2;

%% 


disp("Flow regions that can be solved")
disp("  ")

disp(" 1 --> for Expansion fan effected region ")

disp(" ")

n = input("Enter the region number n = ");

if n ==1

disp(" ");

% -----------------------------------------------------------------
%                        Expansion Corner
% -----------------------------------------------------------------

disp("  ")
disp("Reference number for Expansion region")
disp(" ")
disp("1.1 --> Expansion Corner Plateau Region, Transonic and Supersonic Flow")
disp("1.2 -->  Expansion corner Reattachement Region, Transonic Flow")
disp("  ")

i = input("Enter the reference no for the region i = ");

xe = 6.4;
Rexe = (u_infty*xe)/nu;

delta_star_e = (0.0371/(Rexe)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xe;

omega_e = (2*pi*f*delta_star_e)/u_infty;

% ====================================================================
% Expansion Corner Plateau Region, Transonic and Supersonic Flow
% ====================================================================

if i == 1.1

C_exp = [3;6];
factor = C_exp;

Prms_ratio_exp = 0.040/ (1 + M^2);
P = Prms_ratio_exp;

% =========================================================
%   Expansion corner Reattachement Region, Transonic Flow
% ==========================================================

elseif i == 1.2

    
CRT_exp = [9;15];
factor = CRT_exp;

Prms_RT_exp = 0.16/(1 + M^2);
P = Prms_RT_exp;

end

markers = {'o', '*'};
colors = lines(length(factor));

figure(1)
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================
for idx = 1:2
    G(idx,:) = (4*(P)^2*factor(idx)*(F^(1.433)./(1 + factor(idx)^2*F^(2.867).*omega_e.^2)))*(q^2*delta_star_e/u_infty);
    figure(1)
    %plot(f,G(idx,:),"Color","b","LineWidth",0.7,"Marker","o");
    plot(f, G(idx,:), ...
        'Color', colors(idx,:), ...
        'LineWidth', 0.7, ...
        'Marker', markers{idx}, ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('Factor = %d', factor(idx)));
    hold on
    grid on
    ax = gca;
    ax.XMinorGrid = "on";
    ax.YMinorGrid = "on";
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.TickLength = [0.02 0.01];  % [major minor]
    xlabel("Frequency (in Hz)")
    ylabel("(PSD)")
    legend("C = 3(actual)","C = 6")
    if i == 1.1
        title("Expansion Corner Plateau Region, Transonic and Supersonic Flow")
    elseif i ==1.2
        title("Expansion corner Reattachement Region, Transonic Flow")
    end
    pbaspect([1,1,1]);
    %daspect([50,2,1]);
    xlim([-1000 12000]);
end
end

figure(2)
for jdx = 1:2
    PSD = G'; 
    p_ref = 20e-6;                       % Reference sound pressure level (20 µPa)
    %delta_f = delta_f';
    del = delta_f';
    P = PSD.*del;
    SPL(:,jdx) = 10 * log10(P(:,jdx) / p_ref^2);
    OASPL(:,jdx) = 10 * log10(sum(10.^(SPL(:,jdx)./10)))
    semilogx(f,SPL(:,jdx),'Color', colors(jdx,:), ...
        'LineWidth', 0.7, ...
        'Marker', markers{jdx}, ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('Factor = %d', factor(jdx)))
    hold on
    grid on
    ylim([100,170]);
    ax = gca;
    ax.XMinorGrid = "on";
    ax.YMinorGrid = "on";
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.TickLength = [0.02 0.01];  % [major minor]
    xlabel("Frequency (in Hz)")
    ylabel("SPL (in dB)")
    if i == 1.1
        title("Expansion Corner Plateau Region, Transonic and Supersonic Flow")
    elseif i ==1.2
        title("Expansion corner Reattachement Region, Transonic Flow")
    end
end

hold on
if i ==1.1
    load("expansion_plateau.txt")
    semilogx(expansion_plateau(:,1),expansion_plateau(:,2),"o","Color","k")
    legend("C = 3(actual)","C = 6","Literature","location","best")
else

    load("expansion_peak.txt")
    semilogx(expansion_peak(:,1), expansion_peak(:,2), "o", "Color", "k")
    legend("C = 9(actual)", "C = 15", "Literature", "location", "best")
end
