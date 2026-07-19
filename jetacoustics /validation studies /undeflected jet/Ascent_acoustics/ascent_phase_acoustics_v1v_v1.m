clear ;
clc;

% =========================
%     Octave Frequency
% =========================

nf=31;                              %total number of frequencies considered
f=ones(nf,1);
f(1)=20;                            % initially it was 20
for i=2:nf
f(i)=f(i-1).*(2.^(1/3));             %1/3 octave band frequencies
end
f_signal = f;
f=f';
f_lr=f.*(2.^(-1/6));
f_up=f.*(2.^(1/6));
delta_f=f_up-f_lr;  

alpha = 0;                             % angle in degree  4
M = 0.95;                              % 0.81
u_infty = 303.610;                     % 264.43           % As per the simulation data
nu = 0.00003221;                       % 0.00001258
T_infty = 253.544;                     % 265.20
rho_infty = 0.5020;                    % 0.252
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;


%% 

% =============================
%       ATTACHED FLOW 
% =============================


x_a = 1.5;                            % 0.9    % x is the distace from the start of the boundary layer growth
Rexa = (u_infty*x_a)/nu;

delta_star_a = (0.0371/(Rexa)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*x_a;

 F = 1 + 0.13*M^2;

%Prms_q = (0.010)/F;                              % Prms = (0.010*q)/F
Prms_q = 0.006/(1 + 0.14*M^2);

% ==========================================
%  Power spectral density for attached flow
% ==========================================


ff_a = (2*pi*f*delta_star_a)/u_infty;

G=(4*(Prms_q)^2*((F)^(1.433)./(1 + F^(2.867)*(ff_a).^2)))*(delta_star_a*q^2/u_infty);                     % G=(4*(Prms)^2*((F)^(1.433)./(1 + F^(2.867)*(ff_a).^2)))*(delta_star_a/u_infty);


figure(1)
subplot(2,1,1)
plot(f,G,"Color","b","LineWidth",0.7,"Marker","o");
grid on
ax = gca;
ax.XMinorGrid = "on";
ax.YMinorGrid = "on";
ax.XMinorTick = "on";
ax.YMinorTick = "on";
ax.TickLength = [0.02 0.01];  % [major minor]
xlabel("Frequency (in Hz)")
ylabel("(PSD)")
title("Attached flow")
pbaspect([3,1,1]);


%% 

% ---------------------------------------------------------------------------------------------
%                              SEPERATED FLOW AND SHOCK WAVES
% ---------------------------------------------------------------------------------------------

%hh_data = load("hammerhead_data.txt");
hh_data = load("v1v_ascent_1.3ma.txt");
figure("Name","Mach")
plot(hh_data(:,2),hh_data(:,1),"LineWidth",0.5,"Color","b");
grid on;
%plot(hh_data(:,2),hh_data(:,1),"linewidth",0.7,"Color","r");
xlabel("X- Location");
ylabel("mach no(M)");
legend("AutoUpdate","on");
grid on
ax = gca;
ax.XMinorTick = "on";
ax.YMinorTick = "on";
ax.XMinorGrid = "on";
ax.YMinorGrid = "on";
%print("Mach","-dpng","-r300");


%% 


%[m1,indx] = max(hh_data(:,1));
m1 = 1.34;                     
xs = 3.5;                       
%indx = find(hh_data(:,2) == xs);
%rho = hh_data(indx,5);
%xs = hh_data(indx,2);                 
%T = hh_data(indx,6);
%u = hh_data(indx,7);
%rho = 0.268125;
% T = 220.261;
% u = 400.191;
Rexs = (u_infty*xs)/nu;
%% 


delta_star_s = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;

ff_s = (2*pi*f*delta_star_s)/u_infty;



% ==========================================================
%     Compression Corner Plateau Region, Transonic Flow
% ==========================================================


CPT_comp = 3;

Prms_comp_Pq = 0.025/(1 + M^2);                         %Prms_comp_P = 0.025*q/(1 + M^2);                                 



% ====================================================================
%      Compression Corner Reattachement Region, Transonic Flow
% ====================================================================


CRT_comp = 9;

Prms_comp_Rq = (0.10)/(1 + M^2);                                           %Prms_comp_R = (0.10*q)/(1 + M^2);                                



% ============================================================
%    Compression Corner Plateau Region, Supersonic FLow
% ============================================================



CPS_comp = 10;

theta = alpha + asind(1/m1);                                                % angles are in degrees



P_ratio = 1/2.4*(2.8*m1^2*sind(theta)^2 - 0.4);


 Prms_tbl_q = (0.006)/F;                                                   %Prms_tbl = (0.006*q)/F;                                                        

 Prms_PSq = Prms_tbl_q*P_ratio;                                            %Prms_PS = Prms_tbl*P_ratio;                                               


% ==========================================================================
%           Compression Corner Separation or Reattachement Shockwave
% ==========================================================================

CRS_comp = 30;

Prms_shock_tbl = -1.181 + 1.713*P_ratio + 0.468*(P_ratio)^2;                 % check this line again


Prms_shock_q = Prms_tbl_q*Prms_shock_tbl;                                    %Prms_shock = Prms_tbl*Prms_shock_tbl;                                        

%Prms_shock_comp = Prms_tbl*q;                                               % Check this line


% ==========================================================
%        Power Spectral Density for Compression Corner
% ==========================================================

term1 = 4*(Prms_comp_Pq)^2;                                                     %term1 = 4*(Prms_shock)^2;                                                           
term2 = CPT_comp*F^1.433./(1 + CPT_comp^2 *F^2.867*ff_s.^2);
term3 = delta_star_s*q^2/u_infty;                                           %term3 = delta_star_s/u_infty;                                                 % 

G_C = term1*term2*term3;

%G_C = 4*(Prms_shock)^2*CRS_comp*(F^(1.433)*delta_star./(1 + (CRS_comp)^2 * F^(2.867)*ff.^2*u_infty));

figure(1)
subplot(2,1,2)
plot(f,G_C,"Color","b","LineWidth",0.7,"Marker","o");
grid on
ax = gca;
ax.XMinorGrid = "on";
ax.YMinorGrid = "on";
ax.XMinorTick = "on";
ax.YMinorTick = "on";
ax.TickLength = [0.02 0.01];  % [major minor]
xlabel("Frequency (in Hz)")
ylabel("(PSD)")
title("Compression Corner separation region")
pbaspect([3,1,1]);
%print("Ascent_results_agnikul_veh","-dpng");



%% 

% -----------------------------------------------------------------
%                        Expansion Corner
% -----------------------------------------------------------------

xe = 2.87;
Rexe = (u_infty*xe)/nu;

delta_star_e = (0.0371/(Rexe)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xe;

ff_e = (2*pi*f*delta_star_e)/u_infty;

% ====================================================================
% Expansion Corner Plateau Region, Transonic and Supersonic Flow
% ====================================================================

C_exp = 3;

Prms_ratio_exp = 0.040/ (1 + M^2);

% =========================================================
%   Expansion corner Reattachement Region, Transonic Flow
% ==========================================================

CRT_exp = 9;


Prms_RT_exp = 0.16/(1 + M^2);


% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================


G_ET = 4*(Prms_ratio_exp)^2*C_exp*(F^(1.433)./(1 + C_exp^2*F^(2.867).*ff_e.^2))*(q^2*delta_star_e/u_infty);


subplot(3,1,3)
plot(f,G_ET,"Color","b","LineWidth",0.7,"Marker","o");
grid on
ax = gca;
ax.XMinorGrid = "on";
ax.YMinorGrid = "on";
ax.XMinorTick = "on";
ax.YMinorTick = "on";
ax.TickLength = [0.02 0.01];  % [major minor]
xlabel("Frequency (in Hz)")
ylabel("(PSD)")
title("Expansion Corner")
pbaspect([3,1,1]);

%print("cc_ec","-dpng","-r300");




%% 





% =============================
%      PSD to SPL conversion
% =============================

PSD = G_C'; 


p_ref = 20e-6; % Reference sound pressure level (20 µPa)
delta_f = delta_f';
P = PSD.*delta_f;

SPL = 10 * log10(P / p_ref^2);

% Display the results
disp('Sound Pressure Level (SPL) in dB:');
%disp(SPL);


semilogx(f,SPL);
xlabel("SPL (in dB)");
ylabel("Frequency (hz)")
grid on
ax = gca;
ax.XMinorGrid = "on";
ax.XMinorTick = "on";
ax.YMinorTick = "on";
pbaspect([3,1,1]);