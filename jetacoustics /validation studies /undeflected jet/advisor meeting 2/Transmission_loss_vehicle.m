function [TL,Com,f_c,f_natural] = Transmission_loss_vehicle(f_signal)

format Long g

rho = 286;             %kg/m^3(7750, SS)
E = 7.65e9;               %E = input('Provide the young modulas of the material( in Pa) :'); %Pa
nu = 0.14;              %nu = input('Poisson ratio of the material :');      
%L = 19.761;                 %input('Length of the sheet(in m) :');(2.5)  length of booster
r = 1.4;
h = 15.86e-3;                %input("Thickness of the material(in m)[..,..,..] :");(20e-03)
C_air = 340 ;            % Velocity of sound in the air( in m/sec)
rho_air = 1.225 ;        % density of air in kg/m^3                         
damping_ratio = 0.005 ;  % damping ratio of composite 
                                                                                         %Com = (12*(1-(nu*nu))*W^4)/(E*h);
                                                                                               %f_signal = linspace(20,10000,20000);% Acosutic Signal Frequency
f_natural = 14.7;
n=length(f_signal);
m =length(h);
TL = [];
for i =1:n
    for j =1:m
        C_m = sqrt(E/(rho*(1-(nu*nu))));                            % Velocity of sound in the sheet material( in m/sec)
        f_c = sqrt(3)*(C_air^2)/(pi*C_m*h);                                    %Coincidence frequency for the sheet in Hz
        %Com = (12*(1-(nu*nu))*W^4)/(E*h(1,j));
         Com = (r / h) * (1 - nu^2) / E;
         if f_signal(1,i)>=f_c                                                                                               
           TL(i) = (10*log10(1+((rho*h(1,j)*f_c*pi)/(C_air*rho_air))^2))+(10*log10(damping_ratio))+(33.2*log10(f_signal(1,i)/f_c))-5.7;
        elseif f_signal(1,i)<f_c && f_signal(1,i)>=f_natural 
           TL(i) = 10*log10(1+(rho*pi*h(1,j)*f_signal(1,i)/(rho_air*C_air))^2) - 5;
        else                                                        
           Tl(i) =  20*log10(1/(2*2*pi*f_signal(1,i)*C_air*rho_air*Com)) - 10*log10(0.23*(10*log10(1 + (1/(2*rho_air*C_air*2*pi*f_signal(1,i)*Com))^2)));
        end
    end

end
