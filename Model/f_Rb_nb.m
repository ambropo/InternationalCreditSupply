function F = f_Rb_nb(Rb)
global f beta lambda yH1 yF1 yH2 yF2 eta
% This is the demand function after substituting for s2
F = -Rb + (1/beta).*...
    (((lambda.*yH1)/(lambda.*yF1+(1-lambda).*(1+eta).*f)).^(1-lambda))./...      % this is s1
    (((lambda.*yH2)/(lambda.*yF2-(1-lambda).*(1+eta).*f.*Rb)).^(1-lambda));     % this is s2