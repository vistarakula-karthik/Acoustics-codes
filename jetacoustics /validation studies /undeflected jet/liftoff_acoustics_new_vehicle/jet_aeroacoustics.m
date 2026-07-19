% 1. The parsed Y-values (Amplitude/Power Data)
y_vals = [1.33E+06, 8.98E+05, 5.91E+05, 3.83E+05, 2.46E+05, 1.57E+05, ...
    9.94E+04, 6.29E+04, 3.98E+04, 2.51E+04, 1.58E+04, 9.98E+03, ...
    6.29E+03, 3.96E+03, 2.50E+03, 1.57E+03, 9.91E+02, 6.25E+02, ...
    3.93E+02, 2.48E+02, 1.56E+02, 9.84E+01, 6.20E+01, 3.90E+01, ...
    2.46E+01, 1.55E+01, 9.76E+00, 6.15E+00];

% 2. Define the frequency limits
fl = 20;
fu = 10000;

% 3. Generate the 28 frequency points
% Since this is acoustics data, it uses logarithmic spacing (1/3 octave bands).
f_vals = logspace(log10(fl), log10(fu), length(y_vals));

% 4. Perform the Numerical Integration using the Trapezoidal Rule
total_area = trapz(f_vals, y_vals);

% 5. Display the results
fprintf('Numerical Integration from %d Hz to %d Hz\n', fl, fu);
fprintf('Total Integrated Value: %.2f\n', total_area);