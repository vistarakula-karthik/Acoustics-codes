clear ;
clc;
%% 

% =========================
%     Octave Frequency
% =========================

%de = 0.410;
%ue = 2755;               % 2755 m/s for transonic regime      % 2911.2381 (supersonic)

nf=31;                              %total number of frequencies considered
f=ones(nf,1);
f(1)=20;                            % initially it was 20
for i=2:nf
f(i)=f(i-1).*(2.^(1/3));             %1/3 octave band frequencies
end
f_signal = f;
f=f';
%stn=f.*de/ue;
f_lr=f.*(2.^(-1/6));
f_up=f.*(2.^(1/6));
delta_f=f_up-f_lr;  
%% 

%f = [20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]';       % octave centered frequencies
alpha = 0;                             % angle in degree  4
M = 0.95;                              % 0.81
Re = 1.286e6;
T_infty = 300;                     % 265.20
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
mu = 0.000016203;
rho_infty = 0.22518;  
nu = mu/rho_infty;                 % 0.00001258
q = 0.5*rho_infty*u_infty^2;
 F = 1 + 0.13*M^2;


%% 


% =============================
%       ATTACHED FLOW 
% =============================


x_a = 0.65024;                            % 0.9
Rexa = (u_infty*x_a)/nu;

delta_star_a = (0.0371/(Rexa)^(0.2)*((9/7 + 0.457*M^2)/(1 + 0.13*M^2)^0.64))*x_a;

 F = 1 + 0.13*M^2;

Prms_q = (0.010)/F;                              % Prms = (0.010*q)/F


% ==========================================
%  Power spectral density for attached flow
% ==========================================


ff_a = (2*pi*f*delta_star_a)/u_infty;

G=(4*(Prms_q)^2*((F)^(1.433)./(1 + F^(2.867)*(ff_a).^2)))*(delta_star_a*q^2/u_infty);                     % G=(4*(Prms)^2*((F)^(1.433)./(1 + F^(2.867)*(ff_a).^2)))*(delta_star_a/u_infty);
G = G';

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
hh_data = load("v1t_transonic_data_final1.txt");
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
m1 = 1.14076;                     % 1.50485
xs = 5.4056600;                       %6.75          %9.9375  5.4057
indx = find(hh_data(:,2) == xs);
rho = hh_data(indx,5);
%xs = hh_data(indx,2);                 %9.9375
T = hh_data(indx,6);
u = hh_data(indx,7);
Rexs = (u_infty*xs)/nu;
%% 


delta_star_s = (0.0371/(Rexs)^(0.2)*((9/7 + 0.457*M^2)/(1 + 0.13*M^2)^0.64))*xs;

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



P_ratio = 1/2.4*(2.8*m1^2*sin(theta)^2 - 0.4);


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



xe = 0.65024;
Rexe = (u_infty*xe)/nu;

delta_star_e = (0.0371/(Rexe)^(0.2)*((9/7 + 0.457*M^2)/(1 + 0.13*M^2)^0.64))*xe;

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


G_ET = 4*(Prms_RT_exp)^2*CRT_exp*(F^(1.433)./(1 + CRT_exp^2*F^(2.867).*ff_e.^2))*(q^2*delta_star_e/u_infty);


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
%% 

SPL = 10 * log10(P / p_ref^2);
OASPL = 10 * log10(sum(10.^(SPL./10)));
%% 

% Display the results
disp('Sound Pressure Level (SPL) in dB:');
%disp(SPL);


plot(f,SPL);
grid on
ax = gca;
ax.XMinorGrid = "on";
ax.XMinorTick = "on";
ax.YMinorTick = "on";
pbaspect([3,1,1]);