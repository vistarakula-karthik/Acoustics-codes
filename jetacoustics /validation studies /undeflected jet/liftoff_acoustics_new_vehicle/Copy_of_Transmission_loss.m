function [TL,Com,f_c,f_natural] = Copy_of_Transmission_loss(f_signal)

format Long g


rho = 7850;             %kg/m^3(7750, SS)
E = 210e9;               %E = input('Provide the young modulas of the material( in Pa) :'); %Pa
nu = 0.275;              %nu = input('Poisson ratio of the material :');      
%L = 19.761;                 %input('Length of the sheet(in m) :');(2.5)  length of booster
%r = 1.4;
W = 13.3;                 %input('Width of the sheet(in m) :') ;(2.5)  2pir
L = 6.3;
h = 30e-3;                %input("Thickness of the material(in m)[..,..,..] :");(20e-03)
C_air = 340 ;            % Velocity of sound in the air( in m/sec)
rho_air = 1.225 ;        % density of air in kg/m^3                         
damping_ratio = 0.0013;  % damping ratio of aluminum
z = rho_air*C_air;                                                                                          
M = rho*h;                                                                                    %Com = (12*(1-(nu*nu))*W^4)/(E*h);
                                                                                               %f_signal = linspace(20,10000,20000);% Acosutic Signal Frequency
f_natural = (pi/(4*sqrt(3)))*sqrt(E/(rho*(1 - nu^2)))*h*((1/L)^2 + (1/W)^2);

n=length(f_signal);
m =length(h);
TL = [];
for i =1:n
    for j =1:m
        C_m = sqrt(E/(rho*(1-(nu*nu))));                            % Velocity of sound in the sheet material( in m/sec)
        f_c = sqrt(3)*(C_air^2)/(pi*C_m*h);                           %Coincidence frequency for the sheet in Hz
        Com = 768*(1 - nu*nu)/(pi^8*E*h^3*((1/L^2) + (1/W^2))^2);

        if f_signal(1,i)>=f_c                                                                                                % if acoustic signal frequency is greater than the coincidence frequency then attenuation in acoustic signal
           TL(i) = abs((10*log10(1+((rho*h(1,j)*f_c*pi)/(z))^2))+(10*log10(damping_ratio))+(33.2*log10(f_signal(1,i)/f_c))-5.7);
           %TL(i) = 10*log10(1 + ((pi*M*f_c)/(z))^2) + 10*log10(damping_ratio) + 33.22*log10(f_signal(1,i)/f_c) - 5.7;
        elseif f_signal(1,i)<f_c && f_signal(1,i)>=f_natural  % If acoustic signal frequent is geater than the natural frequency and less than the coincidence frequenct then attenuation in the acoustic signal
           TL(i) = 10*log10(1+(M*pi*f_signal(1,i)/(z))^2) - 5;
        else                                                         % If acoustic signal frequency is lesser than the acoustic signal frequency
           TL(i) = 10*log10(1+(1/(2*z*2*pi*f_signal(1,i)*Com))^2);
        end
    end

end