clear ;
clc;

% =========================
%     Octave Frequency
% =========================

nf=35;                             
f=ones(nf,1);
f(1)=10;                            
for i=2:nf
f(i)=f(i-1).*(2.^(1/3));             
end
f_signal = f;
f=f';
f_lr=f.*(2.^(-1/6));
f_up=f.*(2.^(1/6));
delta_f=f_up-f_lr;

alpha = 0;                             
M = 2; 
T_infty = 166.66;
u_infty = M*sqrt(1.4*287*T_infty);                     % 264.43           % As per the simulation data
%  u_infty = 713.1940079;
rho_infty = 0.703358164;              
nu = 0.00003214958791;
q = 0.5*rho_infty*u_infty^2;
%q = 18195;
F = 1 + 0.13*M^2;


disp("Flow regions that can be solved")
disp("  ")

disp(" 1 --> for Attached flow region ")

disp(" ")


disp(" ");

% -----------------------------------------------------------------
%                      Attached Region
% -----------------------------------------------------------------

xs =1.5;
%Rexs = (u_infty*xs)/nu;
%delta_star_a = (0.0371/(Rexs)^(0.2)*((9/7 + 0.475*M^2)/(1 + 0.13*M^2)^0.64))*xs;
delta_a = 0.0152;
omega_a = (2*pi*f*delta_a)/u_infty;


C_a = 3;
factor = C_a;

Prms_q_a = 0.010/F;                                                    
P1 = Prms_q_a;


figure(1)
% =================================================
%  Power Spectral Density for Expansion Corner
% =================================================

G= (4*(P1)^2*factor*(F^(1.433)./(1 + factor^2*F^(2.867).*omega_a.^2)));
G = G*(u_infty/q^2*delta_a);
semilogx(omega_a,G,"Color","b","LineWidth",0.7,"Marker","o");
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


