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
T_infty = 253.72544;
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
rho_infty = 0.503358164;              
nu = 0.00003214958791;
q = 0.5*rho_infty*u_infty^2;
F = 1 + 0.13*M^2;



disp("Flow regions that can be solved")
disp("  ")
disp(" 1 --> for shock effected region ")
disp(" ")

n = input("Enter the region number n = ");


if n==1

disp("Reference number for compression region")
disp(" ")
disp("1.1 --> Compression Corner Plateau Region, Transonic Flow")
disp("1.2 --> Compression Corner Reattachement Region, Transonic Flow")


disp(" ");

i = input('Enter the reference number corresponding to the region of interest i = ');

 % ====================================
 %    shock upstream conditions
 % ====================================
                  

xs =3.5;
Rexs = (u_infty*xs)/nu;
delta_star_s = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
%delta_star_s = 0.1;
omega_s = (2*pi*f*delta_star_s)/u_infty;


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
%Prms_comp_Rq = 0.07;
P = Prms_comp_Rq;                                            
end


% ==========================================================
%        Power Spectral Density for Compression Corner
% ==========================================================

term1 = 4*(P)^2;                                                                                                      
term2 = factor*F^1.433./(1 + factor^2 *F^2.867*omega_s.^2);
term3 = delta_star_s*q^2/u_infty;                                                                                           
G = term1*term2*term3;

end

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
if i == 1.1
    title("Compression Corner Plateau Region, Transonic Flow")
elseif i==1.2
    title("Compression Corner Reattachement Region, Transonic Flow")
elseif i ==1.3
    title("Compression Corner Plateau Region, Supersonic Flow")
else
    title("Compression Corner Separation or Reattachement shock wave")
end
pbaspect([3,1,1]);


PSD = G'; 
p_ref = 20e-6;                       % Reference sound pressure level (20 µPa)
delta_f = delta_f';
P = PSD.*delta_f;
SPL = 10*log10(P./p_ref^2);
results = [f' SPL];
T = array2table(results, 'VariableNames', {'Frequency','SPL'});
disp(T)
semilogx(f,SPL,"Marker","o","MarkerFaceColor","b")
grid on

OASPL = 10 * log10(sum(10.^(SPL./10)))

%% 


% alpha = 0;                             
% M = 0.933;
% Re = 1.38E+07;
% T_infty = 254.129;
% u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
% mu = 1.62E-05;
% rho_infty = 0.74;  
% nu = mu/rho_infty;                   
% q = 0.5*rho_infty*u_infty^2;
% F = 1 + 0.13*M^2;
