function [Cest, yorg_predict, nmse] = function_LET1Q2_coeffestimate(y, VV, L, Q, Mdiscard)

y = y(:); N = length(y); 
y = y(Mdiscard+1:N); 

VVorg = VV;
VV = VV(Mdiscard+1:N,:);

% estimate the kernel coefficients in least-squares sense
XtX = VV' * VV;
[EV, ED] = eig(XtX);

TH_sv = 1e-6; 
if min(abs(diag(ED)))>TH_sv,
    ALL_C = inv(XtX) * VV' * y;
else
    [Usvd, Ssvd, Vsvd] = svd(VV); SV = abs(diag(Ssvd));
    Idx_sv = find(SV>=TH_sv, 1, 'last'); 
    Pseudo_invVV = Vsvd(:,1:Idx_sv) * diag(1./SV(1:Idx_sv)) * Usvd(:,1:Idx_sv)';
    ALL_C = Pseudo_invVV * y;
end

yorg_predict = VVorg * ALL_C;
y_predict = VV * ALL_C;
nmse = mean((y - y_predict).^2) / mean(y.^2);

% Laguerre coefficients
tempC = ALL_C; 
Cest.c0 = tempC(1); tempC = tempC(2:end);

c1 = tempC(1:L); tempC = tempC(L+1:end);
Cest.c1 = c1;

if Q==2, 
    c2 = zeros(L, L);
    mmm = 1;
    for kk=1:L
        c2(kk,kk:L) = tempC(mmm:mmm+L-kk)';
        mmm = mmm + L - kk + 1; 
    end
    c2 = (c2 + c2') / 2;
    Cest.c2 = c2; 
end

return



