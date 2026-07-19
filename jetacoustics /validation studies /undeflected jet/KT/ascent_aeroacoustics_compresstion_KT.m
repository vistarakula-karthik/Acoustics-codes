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
M = 1.6;          
T_infty =216.08;       
u_infty = M*sqrt(1.4*287*T_infty);                     
rho_infty = 0.3168841;                     
nu = 4.48E-05;      
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;


disp("Flow regions that can be solved")
disp("  ")

disp(" 1 --> for Shock effected region ")

disp(" ")

n = input("Enter the region number n = ");

if n ==1


disp(" ");

% -----------------------------------------------------------------
%                       Compression Corner
% -----------------------------------------------------------------

disp("Reference number for compression region")
disp(" ")
disp("1.1 --> Compression Corner Plateau Region, Transonic Flow")
disp("1.2 --> Compression Corner Reattachement Region, Transonic Flow")

disp("  ")

i = input("Enter the reference no for the region i = ");



xs = 14.4;
Rexs = (u_infty*xs)/nu;
delta_star_s = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
omega_s = (2*pi*f*delta_star_s)/u_infty;

if i == 1.1

% ==========================================================
%     Compression Corner Plateau Region, Transonic Flow
% ==========================================================

CPT_comp = [3;6;9;12;15;100];
factor = CPT_comp;
Prms_comp_Pq = 0.025/(1 + M^2);                                                    
P1 = Prms_comp_Pq;


elseif i == 1.2

% ====================================================================
%      Compression Corner Reattachement Region, Transonic Flow
% ====================================================================


CRT_comp = [9;30];
factor = CRT_comp;
Prms_comp_Rq = (0.10)/(1 + M^2);
P1 = Prms_comp_Rq;
end


% ==========================================================
%        Power Spectral Density for Compression Corner
% ==========================================================

% term1 = 4*(P)^2;                                                                                                      
% term2 = factor*F^1.433./(1 + factor^2 *F^2.867*omega_s.^2);
% term3 = delta_star_s*q^2/u_infty;                                                                                           
% G = term1*term2*term3;

end

markers = {'o','*','o','*','o','*'};
colors = lines(length(factor));

figure(1)
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================
for idx = 1:2
    G(idx,:) = (4*(P1)^2*factor(idx)*(F^(1.433)./(1 + factor(idx)^2*F^(2.867).*omega_s.^2)))*(q^2*delta_star_s/u_infty);
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
    if i == 1.1
        title(" Compression Corner Plateau Region, Transonic Flow")
        legend("C = 3(actual)","C = 6","C = 9","C = 12","C = 15","C = 18")
    elseif i ==1.2
         title("Compression Corner Reattachement Region, Transonic Flow")
        legend("C = 9(actual)","C = 5")
    end
    pbaspect([1,1,1]);
    %daspect([50,2,1]);
    xlim([-1000 12000]);
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
    ax = gca;
    ax.XMinorGrid = "on";
    ax.YMinorGrid = "on";
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.TickLength = [0.02 0.01];  % [major minor]
    xlabel("Frequency (in Hz)")
    ylim([100 160])
    ylabel("SPL (in dB)")
       if i == 1.1
         title(" Compression Corner Plateau Region, Transonic Flow")
         legend("C = 3(actual)","C = 6","C = 9","C = 12","C = 15","C = 18")
        elseif i ==1.2
         title("Compression Corner Reattachement Region, Transonic Flow")
        legend("C = 9(actual)","C = 12","C = 15","C = 18","C = 21","C = 24")
      end
end

 hold on
load("shock_peak.txt")
semilogx(shock_peak(:,1),shock_peak(:,2),"o","Color","k")
 legend("C = 9(actual)","C = 5","Literature","location","northwest")
