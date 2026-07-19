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

% alpha = 0;                             
% M = 0.95; 
% T_infty = 253.72544;
% u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
% rho_infty = 0.503358164;              
% nu = 0.00003214958791;
% q = 0.5*rho_infty*u_infty^2;
% F = 1 + 0.13*M^2;

alpha = 0;                             
M = 0.95; 
T_infty = 286.8607;
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
rho_infty = 0.91498796;              
nu = 0.0000194882;
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;



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


xe = 52.73;
Rexe = (u_infty*xe)/nu;

delta_star_e = (0.0371/(Rexe)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xe;

omega_e = (2*pi*f*delta_star_e)/u_infty;

% ====================================================================
% Expansion Corner Plateau Region, Transonic and Supersonic Flow
% ====================================================================


if i == 1.1

C_exp = [6;20;30;40;50;60;100];
factor = C_exp;

Prms_ratio_exp = 0.040/ (1 + M^2);
P = Prms_ratio_exp;

% =========================================================
%   Expansion corner Reattachement Region, Transonic Flow
% ==========================================================

elseif i == 1.2

    
CRT_exp = [15;20;30;40;50;60;170];
factor = CRT_exp;

Prms_RT_exp = 0.16/(1 + M^2);
P = Prms_RT_exp;

end

markers = {'o', 's', '^', 'd', 'x', '*',"o"};
colors = lines(length(factor));

figure(1)
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================
for idx = 1:7
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
    legend("C = 3(actual)","C = 20","C = 30","C = 40","C = 50","C = 60","C = 170")
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
for jdx = 1:7
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
   legend("C = 3(actual)","C = 20","C = 30","C = 40","C = 50","C = 60","C = 170","location","best")
    if i == 1.1
        title("Expansion Corner Plateau Region, Transonic and Supersonic Flow")
    elseif i ==1.2
        title("Expansion corner Reattachement Region, Transonic Flow")
    end
end

figure(3)

plot(OASPL,factor,"Marker","o","MarkerFaceColor","r")
grid on
xlabel("OASPL (in dB)")
ylabel("factor ")
if i == 1.1
        title("Expansion Corner Plateau Region, Transonic and Supersonic Flow")
    elseif i ==1.2
        title("Expansion corner Reattachement Region, Transonic Flow")
end


