function Agr = ground_atten(f, dp, hs, hr, Gs, Gm, Gr)


f = f(:);   % force column vector
N = length(f);
Agr = zeros(N,1);

%% -------------------------------
% q-factor
%% -------------------------------
if dp <= 30*(hs + hr)
    q = 0;
else
    q = 1 - (30*(hs + hr)/dp);
end

%% -------------------------------
% Coefficient functions
%% -------------------------------
a_p = @(h) 1.5 ...
    + 3.0 * exp(-0.12*(h-5).^2) .* (1 - exp(-dp/50)) ...
    + 5.7 * exp(-0.09*h.^2)     .* (1 - exp(-2.8e-6*dp.^2));

b_p = @(h) 1.5 ...
    + 8.6 * exp(-0.09*h.^2)     .* (1 - exp(-dp/50));

c_p = @(h) 1.5 ...
    + 14.0 * exp(-0.46*h.^2)    .* (1 - exp(-dp/50));

d_p = @(h) 1.5 ...
    + 5.0 * exp(-0.9*h.^2)      .* (1 - exp(-dp/50));

%% -------------------------------
% Loop over frequencies
%% -------------------------------
for i = 1:N
    switch f(i)
        case 63
            As = -1.5;
            Ar = -1.5;

        case 125
            As = -1.5 + Gs * a_p(hs);
            Ar = -1.5 + Gr * a_p(hr);

        case 250
            As = -1.5 + Gs * b_p(hs);
            Ar = -1.5 + Gr * b_p(hr);

        case 500
            As = -1.5 + Gs * c_p(hs);
            Ar = -1.5 + Gr * c_p(hr);

        case 1000
            As = -1.5 + Gs * d_p(hs);
            Ar = -1.5 + Gr * d_p(hr);

        case {2000, 4000, 8000}
            As = -1.5 * (1 - Gs);
            Ar = -1.5 * (1 - Gr);

        otherwise
            error('Frequency must be standard octave band (63–8000 Hz)');
    end

    % Middle region
    if f(i) == 63
        Am = -3 * q;
    else
        Am = -3 * q * (1 - Gm);
    end

    Agr(i) = As + Am + Ar;
end

end
