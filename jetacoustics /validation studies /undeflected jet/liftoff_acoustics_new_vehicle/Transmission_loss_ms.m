function [TL,Com,f_c,f_natural] = Transmission_loss_ms(f_signal)

format Long g

% Material Properties
%Material_name=stringtype(input('provide the material name :'));
%rho = input("Provide the material density in Kg/m^3 :");


rho = 2400;             %kg/m^3(7750, SS)
E = 20.7e9;               %E = input('Provide the young modulas of the material( in Pa) :'); %Pa
nu = 0.13;              %nu = input('Poisson ratio of the material :');      
L = 17.5;                 %input('Length of the sheet(in m) :');(2.5)
W = 17.5;                 %input('Width of the sheet(in m) :') ;(2.5)
h = 150e-3;                %input("Thickness of the material(in m)[..,..,..] :");(20e-03)
C_air = 340 ;            % Velocity of sound in the air( in m/sec)
rho_air = 1.225 ;        % density of air in kg/m^3                         
damping_ratio = 0.020 ;  % damping ratio of steel 
                                                                                              %Com = (12*(1-(nu*nu))*W^4)/(E*h);
                                                                                               %f_signal = linspace(20,10000,20000);% Acosutic Signal Frequency
f_natural = ((pi*pi*(1+(W/L)^2)))*(((E*(h^3))/(12*(1-(nu*nu))))*(L/(rho*L*h*(W)^4)))^0.5;      % Natural frequency of the sheet
n=length(f_signal);
m =length(h);
TL = [];
for i =1:n
    for j =1:m
        C_m = sqrt(E/(rho*(1-(nu*nu))));                            % Velocity of sound in the sheet material( in m/sec)
        f_c = sqrt(3)*(C_air^2)/(pi*C_m*h(1,j));                             %Coincidence frequency for the sheet in Hz
        Com = (12*(1-(nu*nu))*W^4)/(E*h(1,j));
        
        f_natural = ((pi*pi*(1+(W/L).^2)))*(((E*(h(1,j).^3))/(12*(1-(nu*nu))))*(L/(rho*L*h(1,j)*(W).^4))).^0.5;               % Natural frequency of the sheet
        if f_signal(1,i)>=f_c                                                                                                % if acoustic signal frequency is greater than the coincidence frequency then attenuation in acoustic signal
           TL(i) = abs((10*log10(1+((rho*h(1,j)*f_c*pi)/(C_air*rho_air))^2))+(10*log10(damping_ratio))+(33.2*log10(f_signal(1,i)/f_c))-5.7);
        elseif f_signal(1,i)<f_c && f_signal(1,i)>=f_natural % If acoustic signal frequent is geater than the natural frequency and less than the coincidence frequenct then attenuation in the acoustic signal
           TL(i) = 10*log10(1+(rho*2*pi*h(1,j)*f_signal(1,i)/(2*rho_air*C_air))^2);
        else                                                         % If acoustic signal frequency is lesser than the acoustic signal frequency
           TL(i) = 10*log10(1+(1/(2*rho_air*C_air*2*pi*f_signal(1,i)*Com))^2);
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