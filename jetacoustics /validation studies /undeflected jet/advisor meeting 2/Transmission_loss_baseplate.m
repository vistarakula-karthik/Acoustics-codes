function [TL,Com,f_c,f_natural] = Transmission_loss_baseplate(f_signal)

format Long g

  
rho = 2780;                   %kg/m^3(7750, SS)
E = 71e9;                  %E = input('Provide the young modulas of the material( in Pa) :'); %Pa
nu = 0.3;                   %nu = input('Poisson ratio of the material :');      
D = 2.8;
h = 12e-3;                  %input("Thickness of the material(in m)[..,..,..] :");(20e-03)
C_air = 334 ;              % Velocity of sound in the air( in m/sec)
rho_air = 1.12 ;          % density of air in kg/m^3                         
damping_ratio = 0.01;      % damping ratio of aluminum
                                                                                         
                                                                                              
f_natural = 15.31;                          %6.3357 %10.42
n=length(f_signal);
m =length(h);
TL = [];
for i =1:n
    for j =1:m
        C_m = sqrt(E/(rho*(1-(nu*nu))));                            % Velocity of sound in the sheet material( in m/sec)
        f_c = sqrt(3)*(C_air^2)/(pi*C_m*h);                                    %Coincidence frequency for the sheet in Hz
        Com = (3*D^4*(1 - nu^2))/256*E*h^3;
        
        if f_signal(1,i)>=f_c                                                                                               
           TL(i) = (10*log10(1+((rho*h(1,j)*f_c*pi)/(C_air*rho_air))^2))+(10*log10(damping_ratio))+(33.2*log10(f_signal(1,i)/f_c))-5.7;
        elseif f_signal(1,i)<f_c && f_signal(1,i)>=f_natural 
           TL(i) = 10*log10(1+(rho*pi*h(1,j)*f_signal(1,i)/(rho_air*C_air))^2) - 5;
        else                                                        
           Tl(i) =  20*log10(1/(2*2*pi*f_signal(1,i)*C_air*rho_air*Com)) - 10*log10(0.23*(10*log10(1 + (1/(2*rho_air*C_air*2*pi*f_signal(1,i)*Com))^2)));
        end
    end

end
                                                        