d=1;
s2e=1;
s2q=1e-4;


H=[0 0 1; d^2 d 1; 4*d^2 2*d 1];
H2=[0 0 1; 0.5 0.5 1; 2 1 1];
T=diag([1/(2*d^2) 1/(2*d) 1]);

% Normal
(s2e+s2q)*inv(H'*H)

% Quant
s2q*T*inv(H2'*H2)*T'+s2e*(T*T')

s2q*T*inv(H2'*H2)*T'
s2e*(T*T')