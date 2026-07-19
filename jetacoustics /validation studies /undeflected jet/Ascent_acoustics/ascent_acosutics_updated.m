clear ;
clc;

% =========================
%     Octave Frequency
% =========================

nf=31;                             
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
M = 0.95; 
Re = 1.286e6;
T_infty = 300;
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
mu = 0.000016203;
rho_infty = 0.22518;  
nu = mu/rho_infty;   
%u_infty = 303.610;                     
%nu = 0.00003221;              
%T_infty = 253.544;                     
%rho_infty = 0.5020;                 
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;
%% 

disp("Enter the region")
disp("  ")

disp(" 1 --> for attached flow region ")
disp(" 2 --> for shock effected region ")
disp(" 3 --> for Expansion fan effected region ")
disp(" ")

n = input("Enter the region number n = ");

if n ==1

% =============================
%       ATTACHED FLOW 
% =============================

x_a = 1.5;                            
Rexa = (u_infty*x_a)/nu;

delta_star_a = (0.0371/(Rexa)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*x_a;

 F = 1 + 0.13*M^2;

Prms_q = (0.010)/F;                              

% ==========================================
%  Power spectral density for attached flow
% ==========================================


ff_a = (2*pi*f*delta_star_a)/u_infty;

G=(4*(Prms_q)^2*((F)^(1.433)./(1 + F^(2.867)*(ff_a).^2)))*(delta_star_a*q^2/u_infty);                    

figure(1)
plot(f,G,"Color","b","LineWidth",0.7,"Marker","o");
grid on
ax = gca;
ax.XMinorGrid = "on";
ax.YMinorGrid = "on";
ax.XMinorTick = "on";
ax.YMinorTick = "on";
ax.TickLength = [0.02 0.01]; 
xlabel("Frequency (in Hz)")
ylabel("(PSD)")
title("Attached flow")
pbaspect([3,1,1]);



elseif n==2

disp("Reference number for compression region")
disp(" ")
disp("1.1 --> Compression Corner Plateau Region, Transonic Flow")
disp("1.2 --> Compression Corner Reattachement Region, Transonic Flow")
disp("1.3 --> Compression Corner Plateau Region, Supersonic Flow")
disp("1.4 --> Compression Corner Separation or Reattachement shock wave")

disp(" ");

i = input('Enter the reference number corresponding to the region of interest i = ');

m1 = 1.34;                     
xs = 3.5;                       
Rexs = (u_infty*xs)/nu;
delta_star_s = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
ff_s = (2*pi*f*delta_star_s)/u_infty;
theta = alpha + asind(1/m1);                                                
P_ratio = 1/2.4*(2.8*m1^2*sind(theta)^2 - 0.4);

if i == 1.1

% ==========================================================
%     Compression Corner Plateau Region, Transonic Flow
% ==========================================================

CPT_comp = 3;
factor = CPT_comp;

Prms_comp_Pq = 0.025/(1 + M^2);                                                    
P = Prms_comp_Pq;

elseif i == 1.2

% ====================================================================
%      Compression Corner Reattachement Region, Transonic Flow
% ====================================================================


CRT_comp = 9;
factor = CRT_comp;

Prms_comp_Rq = (0.10)/(1 + M^2);                                                                
P = Prms_comp_Rq;


% ============================================================
%    Compression Corner Plateau Region, Supersonic FLow
% ============================================================


elseif i == 1.3
CPS_comp = 10;
factor = CPS_comp;

 Prms_tbl_q = (0.006)/F;                                                                                                    

 Prms_PSq = Prms_tbl_q*P_ratio;                                                                                
 P = Prms_PSq;

% ==========================================================================
%           Compression Corner Separation or Reattachement Shockwave
% ==========================================================================

elseif i == 1.4
CRS_comp = 30;
factor = CRS_comp;

Prms_shock_tbl = -1.181 + 1.713*P_ratio + 0.468*(P_ratio)^2;                 


Prms_shock_q = Prms_tbl_q*Prms_shock_tbl;                                                                      
P = Prms_shock_q;
                                            
end


% ==========================================================
%        Power Spectral Density for Compression Corner
% ==========================================================

term1 = 4*(P)^2;                                                                                                      
term2 = factor*F^1.433./(1 + factor^2 *F^2.867*ff_s.^2);
term3 = delta_star_s*q^2/u_infty;                                                                                           
G = term1*term2*term3;





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
title("Compression Corner separation region")
pbaspect([3,1,1]);
%print("Ascent_results_agnikul_veh","-dpng");


elseif n==3
% -----------------------------------------------------------------
%                        Expansion Corner
% -----------------------------------------------------------------

disp("  ")
disp("Reference number for Expansion region")
disp(" ")
disp("2.1 --> Expansion Corner Plateau Region, Transonic and Supersonic Flow")
disp("2.2 -->  Expansion corner Reattachement Region, Transonic Flow")
disp("  ")

i = input("Enter the reference no for the region i = ");

xe = 2.87;
Rexe = (u_infty*xe)/nu;

delta_star_e = (0.0371/(Rexe)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xe;

ff_e = (2*pi*f*delta_star_e)/u_infty;

% ====================================================================
% Expansion Corner Plateau Region, Transonic and Supersonic Flow
% ====================================================================


if i == 2.1

C_exp = 3;
factor = C_exp;

Prms_ratio_exp = 0.040/ (1 + M^2);
P = Prms_ratio_exp;

% =========================================================
%   Expansion corner Reattachement Region, Transonic Flow
% ==========================================================

elseif i == 2.2

CRT_exp = 9;
factor = CRT_exp;

Prms_RT_exp = 0.16/(1 + M^2);
P = Prms_RT_exp;

end
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================


G = 4*(P)^2*factor*(F^(1.433)./(1 + factor^2*F^(2.867).*ff_e.^2))*(q^2*delta_star_e/u_infty);


subplot(3,1,3)
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
title("Expansion Corner")
pbaspect([3,1,1]);
end

%print("cc_ec","-dpng","-r300");



% =============================
%      PSD to SPL conversion
% =============================
% 
PSD = G'; 


p_ref = 20e-6; % Reference sound pressure level (20 µPa)
delta_f = delta_f';
P = PSD.*delta_f;
%% 

SPL = 10 * log10(P / p_ref^2);
OASPL = 10 * log10(sum(10.^(SPL./10)))
%% 

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


