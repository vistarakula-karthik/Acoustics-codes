clear ;
clc;

% =========================
%     Octave Frequency
% =========================

nf=15;                             
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
 M =0.8;            % 0.95; 
 T_infty = 265.9574468;
 u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
 rho_infty = 1.481178627;              
 nu = 0.0000168025518/rho_infty;
 q = 0.5*rho_infty*u_infty^2;
 F = 1 + 0.13*M^2;

 %  alpha = 0;                             
 % M =1.05;            % 0.95; 
 % T_infty = 245.800;
 % u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
 % rho_infty = 1.21628;              
 % nu = 0.000015773/rho_infty;
 % q = 0.5*rho_infty*u_infty^2;
 % F = 1 + 0.13*M^2;
% 
% alpha = 0;                             
%  M = 3.2;            % 0.95; 
%  T_infty = 212.390;
%  u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
%  rho_infty = 0.2944;              
%  nu = 0.000047488;
%  q = 0.5*rho_infty*u_infty^2;
%  F = 1 + 0.13*M^2;

%  alpha = 0;                             
% M = 0.95; 
% T_infty = 253.72544;
% u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
% rho_infty = 0.503358164;              
% nu = 0.00003214958791;
% q = 0.5*rho_infty*u_infty^2;
% F = 1 + 0.13*M^2;


%  alpha = 0;                             
% M = 0.95; 
% T_infty = 286.8607;
% u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
% rho_infty = 0.91498796;              
% nu = 0.0000194882;
% q = 0.5*rho_infty*u_infty^2;
% F = 1 + 0.13*M^2;

  % alpha = 0;                             
  % M = 1.3; 
  % T_infty = 221.46;
  % u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
  % rho_infty = 0.3388;              
  % nu = 0.000042737;
  % q = 0.5*rho_infty*u_infty^2;
  % F = 1 + 0.13*M^2;
  
    % alpha = 0;                             
    % M = 0.95; 
    % T_infty = 254.129;
    % u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
    % rho_infty = 1.3219;              
    % nu = 0.00001225;
    % q = 0.5*rho_infty*u_infty^2;
    % F = 1 + 0.13*M^2;

%% 

disp("Flow regions that can be solved")
disp("  ")

disp(" 1 --> for Attached flow region ")

disp(" ")


disp(" ");

% -----------------------------------------------------------------
%                      Attached Region
% -----------------------------------------------------------------

xs =0.79;
Rexs = (u_infty*xs)/nu;
delta_star_a = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
omega_a = (2*pi*f*delta_star_a)/u_infty;


C_a = 12;
factor = C_a;
Prms_q_a = 0.010/F;
Prms_q_s = max(0.026,(0.041/ (1 + 1.606*M^2)));
%Prms_q_s = 0.04;
P1 = Prms_q_s;

markers = {'o'};
colors = lines(length(factor));

figure(1)
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================
for idx = 1:1
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
for jdx = 1:1
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
    title(" Attached flow region")  
    ylim([100 170])
end


hold on
load("attached.txt")
semilogx(attached(:,1),attached(:,2),"-o","Color","r")
legend("C = 1(actual)","Literature","location","northwest")