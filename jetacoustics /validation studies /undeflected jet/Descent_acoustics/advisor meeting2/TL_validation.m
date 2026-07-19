

nf = 25;                      % Total number of bands
f = zeros(nf,1);
n = 3;                        % For 1/3 octave bands, n = 3
f(1) = 31.5;                    % Start frequency (could also use standardized frequencies)

% Generate 1/3 octave center frequencies
for i = 2:nf
    f(i) = f(i-1) * 2^(1/n);
end

f = f';                       % Convert to row vector (optional)
f_signal = f;

% Strouhal number calculation (if needed)

% Correct lower and upper frequency bounds
f_lr = f .* 2.^(-1/(2*n));
f_up = f .* 2.^(1/(2*n));
delta_f = f_up - f_lr;


format Long g


rho = 2700;             %kg/m^3(7750, SS)
E = 70e9;               %E = input('Provide the young modulas of the material( in Pa) :'); %Pa
nu = 0.3;              %nu = input('Poisson ratio of the material :');      
L = 2.5;                 %input('Length of the sheet(in m) :');(2.5)
W = 1.5;                 %input('Width of the sheet(in m) :') ;(2.5)
h = 8e-3;                %input("Thickness of the material(in m)[..,..,..] :");(20e-03)
C_air = 334;            % Velocity of sound in the air( in m/sec)
rho_air = 1.12;         % density of air in kg/m^3                         
damping_ratio = 0.001;  % damping ratio of steel 
z = rho_air*C_air;                                                                                          
 M = 21.6;

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
        elseif f_signal(1,i)<f_c  % If acoustic signal frequent is geater than the natural frequency and less than the coincidence frequenct then attenuation in the acoustic signal
           TL(i) = 10*log10(1+(M*pi*f_signal(1,i)/(z))^2) - 5;
        else                                                         % If acoustic signal frequency is lesser than the acoustic signal frequency
           TL(i) = 10*log10(1+(1/(2*z*2*pi*f_signal(1,i)*Com))^2);
        end
    end

end


semilogx(f,TL,"o","MarkerSize",8,"MarkerFaceColor","b");
grid on
data = load("al_tl.txt");
hold on
semilogx(data(:,1),data(:,2),"o","MarkerSize",8,"MarkerFaceColor","r");
xlabel("Frequency (in Hz)")
ylabel("TL (in dB)")
legend("current","literature")
%% 

data = load("ascent_validation .txt");
hold on
semilogx(data(:,1),data(:,2),"*","MarkerSize",10);
grid on

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