function [TL,Com,f_c,f_natural] = Transmission_loss_baseplate(f_signal)

format Long g


rho = 2840;             %kg/m^3(7750, SS)
E = 73.8e9;               %E = input('Provide the young modulas of the material( in Pa) :'); %Pa
nu = 0.33;              %nu = input('Poisson ratio of the material :');      
%L = 19.761;                 %input('Length of the sheet(in m) :');(2.5)  length of booster
r = 1.4;
%W = 2*pi*r;                 %input('Width of the sheet(in m) :') ;(2.5)  2pir
h = 4e-3;                %input("Thickness of the material(in m)[..,..,..] :");(20e-03)
C_air = 340 ;            % Velocity of sound in the air( in m/sec)
rho_air = 1.225 ;        % density of air in kg/m^3                         
damping_ratio = 0.04;  % damping ratio of aluminum
                                                                                         %Com = (12*(1-(nu*nu))*W^4)/(E*h);
                                                                                               %f_signal = linspace(20,10000,20000);% Acosutic Signal Frequency
f_natural = 6.370;
%f_natural = ((pi*pi*(1+(W/L)^2)))*(((E*(h^3))/(12*(1-(nu*nu))))*(L/(rho*L*h*(W)^4)))^0.5;        % Natural frequency of the sheet
n=length(f_signal);
m =length(h);
TL = [];
for i =1:n
    for j =1:m
        C_m = sqrt(E/(rho*(1-(nu*nu))));                            % Velocity of sound in the sheet material( in m/sec)
        f_c = sqrt(3)*(C_air^2)/(pi*C_m*h(1,j));                             %Coincidence frequency for the sheet in Hz
        %Com = (12*(1-(nu*nu))*W^4)/(E*h(1,j));
         Com = (r / h) * (1 - nu^2) / E;
        
        %f_natural = ((pi*pi*(1+(W/L).^2)))*(((E*(h(1,j).^3))/(12*(1-(nu*nu))))*(L/(rho*L*h(1,j)*(W).^4))).^0.5;               % Natural frequency of the sheet
        if f_signal(1,i)>=f_c                                                                                                % if acoustic signal frequency is greater than the coincidence frequency then attenuation in acoustic signal
           TL(i) = (10*log10(1+((rho*h(1,j)*f_c*pi)/(C_air*rho_air))^2))+(10*log10(damping_ratio))+(33.2*log10(f_signal(1,i)/f_c))-5.7;
        elseif f_signal(1,i)<f_c && f_signal(1,i)>=f_natural % If acoustic signal frequent is geater than the natural frequency and less than the coincidence frequenct then attenuation in the acoustic signal
           TL(i) = 10*log10(1+(rho*2*pi*h(1,j)*f_signal(1,i)/(2*rho_air*C_air))^2);
        else                                                         % If acoustic signal frequency is lesser than the acoustic signal frequency
           %TL(i) = 10*log10(1+(1/(2*rho_air*C_air*2*pi*f_signal(1,i)*Com))^2);
           Tl(i) =  20*log10(1/(2*2*pi*f_signal(1,i)*C_air*rho_air*Com)) - 10*log10(0.23*(10*log10(1 + (1/(2*rho_air*C_air*2*pi*f_signal*Com))^2)));
        end
    end

end
%x=smooth()
% hold on
% plot(TL(1:n,1),TL(1:n,2),'LineWidth',2)
% plot(TL(1:n,1),TL(1:n,3),'LineWidth',2)
% plot(TL(1:n,1),TL(1:n,4),'LineWidth',2)
% ax = gca;
% ax.FontSize = 16;
% lgd=legend(["10mm","8mm"," 5mm"]);
% fontsize(lgd,18,"points")
% xlabel("Frequency of Acoustic Signal( in Hz)","FontSize",18)
% ylabel("Transmission Loss( in dB)","FontSize",18)
% hold off
%TL;
%end