function [Cest, yorg_predict, nmse] = function_LET2Q2_coeffestimate(y, VV, L1, L2, Q, Mdiscard)

y = y(:); N = length(y); 
y = y(Mdiscard+1:N); 

N1 = L1; 
if Q==2, N1_2nd = L1*(L1+1)/2; end
N2 = L2; 
if Q==2, N2_2nd = L2*(L2+1)/2; end
if Q==2, N12_2nd = L1 * L2; end

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
%Cest.all = ALL_C; 

% Laguerre coefficients
tempC = ALL_C; 
Cest.c0 = tempC(1); tempC = tempC(2:end);

c10 = tempC(1:N1); tempC = tempC(N1+1:end);
c01 = tempC(1:N2); tempC = tempC(N2+1:end);
Cest.c10 = c10; Cest.c01 = c01;

if Q==2, 
    buf_c20 = tempC(1:N1_2nd); tempC = tempC(N1_2nd+1:end);
    c20 = zeros(L1, L1);
    mmm = 1;
    for kk=1:L1
        c20(kk,kk:L1) = buf_c20(mmm:mmm+L1-kk)';
        mmm = mmm + L1 - kk + 1; 
    end
    c20 = (c20 + c20') / 2;
    Cest.c20 = c20;

    buf_c02 = tempC(1:N2_2nd); tempC = tempC(N2_2nd+1:end);
    c02 = zeros(L2, L2);
    mmm = 1;
    for kk=1:L2
        c02(kk,kk:L2) = buf_c02(mmm:mmm+L2-kk)';
        mmm = mmm + L2 - kk + 1; 
    end
    c02 = (c02 + c02') / 2;
    Cest.c02 = c02; 

    if ~isempty(tempC), 
        buf_c11(1:N12_2nd) = tempC(1:end);
        c11 = zeros(L1, L2);
        for kk=1:L1,
            c11(kk,:) = buf_c11((kk-1)*L2+1:kk*L2);
        end
        Cest.c11 = c11; 
    end
end

return











