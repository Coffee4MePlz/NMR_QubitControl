function [XCoeff, YCoeff, CCoeff] = planefit(auxx,auxy,auxz)
% Fitting a best-fit line to data, both noisy and non-noisy
% Plane Equation
%
% z = XCoeff*X + YCoeff*Y + CCoeff
%
Xcolv = auxx(:); % Make X a column vector
Ycolv = auxy(:); % Make Y a column vector
Zcolv = auxz(:); % Make Z a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term
Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
XCoeff = Coefficients(1); % X coefficient
YCoeff = Coefficients(2); % X coefficient
CCoeff = Coefficients(3); % constant term
