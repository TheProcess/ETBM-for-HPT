function xs = xsegregate(n)

% Written by Yajun Wei c2004
% calculating the Sb segregation value x(n) in the nth As layer in InAs layers

% adjustable parameters:
xi = 0.124; % initial excessive antimony on the finished GaSb layer surface
x0 = 0.010; % excessive antimony from the background Sb pressure in the growth chamber
R = 0.67; % the incorporation ratio of the previous Sb into the next atomic layer of InAs

% the resulting ratio of Sb atoms in the nth As layer in the SL InAs layer
xs = xi * R^(n-1) * (1-R) + x0 * (1-R^n);
