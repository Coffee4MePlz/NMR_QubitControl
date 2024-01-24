function [m2, b2] = linefit(x,yn)
% Fitting a best-fit line to data, both noisy and non-noisy
% Determine coefficients for noisy line yn=m2*x+b2
Xcolv = x(:); % Make X a column vector
Yncolv = yn(:); % Make Yn a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term
NoisyCoeffs = [Xcolv Const]\Yncolv; % Find the coefficients
m2 = NoisyCoeffs(1);
b2 = NoisyCoeffs(2);