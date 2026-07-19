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
% M = 1.050354675;         %0.8;             % 0.95; 
% T_infty = 244.7785758;       %264.5762469;     %253.72544;
% u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
% rho_infty = 0.44783051;       %0.593246144;        %0.503358164;              
% nu = 0.0000351040472;      %0.00002820623496;          %0.00003214958791;
% q = 0.5*rho_infty*u_infty^2;

 alpha = 0;                             
M = 1.34; 
T_infty = 219.1693957;
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
rho_infty = 0.3282247;              
nu = 0.000043732993;
q = 0.5*rho_infty*u_infty^2;
%q = 30000;
F = 1 + 0.13*M^2;


disp("Flow regions that can be solved")
disp("  ")

disp(" 1 --> for Attached flow region ")

disp(" ")


disp(" ");

% -----------------------------------------------------------------
%                      Attached Region
% -----------------------------------------------------------------

xs = 28.3;
Rexs = (u_infty*xs)/nu;
delta_star_a = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
omega_a = (2*pi*f*delta_star_a)/u_infty;


C_a = [1;2;3;4;9];
factor = C_a;

%Prms_q_a = 0.010/F;                                                    
Prms_q_a = min(0.026,(0.041/1 + 1.606*M^2));                                                    
P1 = Prms_q_a;

markers = {'o', 's', '^', 'd', 'x', '*'};
colors = lines(length(factor));

figure(1)
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================
for idx = 1:5
    G(idx,:) = (4*(P1)^2*factor(idx)*(F^(1.433)./(1 + factor(idx)^2*F^(2.867).*omega_a.^2)))*(q^2*delta_star_a/u_infty);
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
    title(" Attached flow region")
    legend("C = 1(actual)","C = 2","C = 3","C = 4","C = 9")
    pbaspect([1,1,1]);
    xlim([-1000 12000]);
end


figure(2)
for jdx = 1:5
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
    ax = gca;
    ax.XMinorGrid = "on";
    ax.YMinorGrid = "on";
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.TickLength = [0.02 0.01];  % [major minor]
    xlabel("Frequency (in Hz)")
    ylabel("SPL (in dB)")
    legend("C = 1(actual)","C = 2","C = 3","C = 4","C = 9","location","northwest")
    title(" Attached flow region")  
end

figure(3)

plot(OASPL,factor,"Marker","o","MarkerFaceColor","r")
grid on
xlabel("OASPL (in dB)")
ylabel("Frequency (in Hz)")
title("Attached flow region")

