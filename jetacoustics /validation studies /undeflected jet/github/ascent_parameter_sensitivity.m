clear ;
clc;

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

alpha = 0;                             
M = 0.7; 
T_infty = 253.72544;
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
rho_infty = 0.503358164;              
nu = 0.00003214958791;
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;

%% 

disp("Flow regions that can be solved")
disp("  ")

disp(" 1 --> for Attached flow region ")

disp(" ")


disp(" ");

% -----------------------------------------------------------------
%                      Attached Region
% -----------------------------------------------------------------

xs =1.5;
Rexs = (u_infty*xs)/nu;
delta_star_a = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
omega_a = (2*pi*f*delta_star_a)/u_infty;


C_a = [1;2;3;4;9];
factor = C_a;

Prms_q_a = 0.010/F;                                                    
P1 = Prms_q_a;

markers = {'o', 's', '^', 'd', 'x', '*'};
colors = lines(length(factor));

figure("Name","PSD1")
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


figure("Name","SPL1")
for jdx = 1:5
    PSD = G'; 
    p_ref = 20e-6;                       % Reference sound pressure level (20 µPa)
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

figure("Name","OASPL1")

plot(OASPL,factor,"Marker","o","MarkerFaceColor","r")
grid on
xlabel("OASPL (in dB)")
ylabel("Frequency (in Hz)")
title("Attached flow region")
%% 


delta_star_a = (0.0033:0.001:0.0083)';

C_a = 9;
factor = C_a;

Prms_q_a = 0.010/F;                                                    
P1 = Prms_q_a;

markers = {'o', 's', '^', 'd', 'x', '*'};
colors = lines(length(delta_star_a));

figure("Name","PSD2")
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================
for idx = 1:6
    omega_a(idx,:) = (2*pi*f*delta_star_a(idx))/u_infty;
    G(idx,:) = (4*(P1)^2*factor*(F^(1.433)./(1 + factor^2*F^(2.867).*omega_a(idx,:).^2)))*(q^2*delta_star_a(idx)/u_infty);
    figure(1)
    plot(f, G(idx,:), ...
        'Color', colors(idx,:), ...
        'LineWidth', 0.7, ...
        'Marker', markers{idx}, ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('delta_star = %d', delta_star_a(idx)));
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
    legend("delta = 0.0033(actual)","delta = 0.0043","delta = 0.0053","delta = 0.0063","delta = 0.0073","delta = 0.0083")
    pbaspect([1,1,1]);
    xlim([-1000 12000]);
end


figure("Name","SPL2")
for jdx = 1:6
    PSD = G'; 
    p_ref = 20e-6;                       % Reference sound pressure level (20 µPa)
    del = delta_f';
    P = PSD.*del;
    SPL(:,jdx) = 10 * log10(P(:,jdx) / p_ref^2);
    OASPL(:,jdx) = 10 * log10(sum(10.^(SPL(:,jdx)./10)))
    semilogx(f,SPL(:,jdx),'Color', colors(jdx,:), ...
        'LineWidth', 0.7, ...
        'Marker', markers{jdx}, ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('delta = %d', delta_star_a(jdx)))
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
    legend("delta = 0.0033(actual)","delta = 0.0043","delta = 0.0053","delta = 0.0063","delta = 0.0073","delta = 0.0083","location","northwest")
    title(" Attached flow region")  
end

figure("Name","OASPL2")

plot(OASPL,delta_star_a,"Marker","o","MarkerFaceColor","r")
grid on
xlabel("OASPL (in dB)")
ylabel("delta_star_a (in Hz)")
title("Attached flow region")

%% 


