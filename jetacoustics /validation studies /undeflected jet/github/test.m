% Example parameters
Uinf = 298;          % m/s
M = 0.933;
x =0.04053;
q = 32887;           % Pa
Rex= 5.520e+05;
 ratio = 0.0371 * (Rex).^(-0.2) .*((9/7 + 0.475*M.^2) ./ (1 + 0.13*M.^2).^0.64);
    delta_star = ratio .* x;
%delta_star = 1.6955e-04;  % m
L = 0.28;
prms =    1.7595e+03;           % Pa
Cexp = 3;          % correction factor
F =1.1132;               % nondimensional frequency
f = linspace(40,480,12); % frequency range
Omega = 2*pi*L*f/Uinf;
%Omega = linspace(1,10,200);
%Omega = linspace(0.1,10,200); % frequency range

% Compute G(f)
Gvals = Gf_of_Omega(Omega, Uinf, q, delta_star, prms, Cexp, F);

% Plot
semilogx(f, Gvals, 'Marker',"o")
xlabel('\Omega','FontSize',14)
ylabel('G(f)','FontSize',14)
grid on


function G = Gf_of_Omega(Omega, Uinf, q, delta_star, prms, Cexp, F)
% Gf_of_Omega - Computes G(f) as a function of Omega
%
% Inputs:
%   Omega      - nondimensional frequency (vector or scalar)
%   Uinf       - freestream velocity
%   q          - dynamic pressure
%   delta_star - displacement thickness (δ*)
%   prms       - rms pressure (experimental)
%   Cexp       - experimental correction factor
%   F          - frequency parameter
%
% Output:
%   G          - spectral gain G(f) as function of Omega

    % Term from prms/q
    term1 = 4 * (prms/q)^2;

    % Frequency-dependent part
    term2 = Cexp .* (F.^1.433) ./ (1 + (Cexp.^2) .* (F.^2.867) .* (Omega.^2));

    % Scaling
    term3 = (q^2 * delta_star) / Uinf;

    % Final expression
    G = term1 .* term2 .* term3;
end
