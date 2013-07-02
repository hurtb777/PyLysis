function [yest, nmse, ANFs, uu] = function_estimate_3rd_order_ANFs(x, y, PDMs)

x = x(:); y = y(:); N = length(y);
Mdiscard = size(PDMs, 1); 

uu = function_LET1PDMoutputs(x, PDMs);
    
% nonlinear ANFs
ANFs_order = 3; 
Npdms = size(PDMs, 2); 

ANFs_coeff = zeros(ANFs_order * Npdms, 1); 
Nunknowns1 = size(ANFs_coeff, 1);

Mtx = zeros(N, Nunknowns1); 
cnt = 0; 
for mmm=1:Npdms,
    for kk=1:ANFs_order,
        cnt = cnt + 1;
        Mtx(:,cnt) = uu(:,mmm).^kk;
    end
end

Mtx = [ones(N, 1), Mtx];

Mtx_org = Mtx;
Mtx = Mtx(Mdiscard+1:N,:);
y = y(Mdiscard+1:N);

[EV, ED] = eig(Mtx' * Mtx);
min(abs(diag(ED)));
TH_sv = 1e-2; 
if min(abs(diag(ED)))>TH_sv,
    ANFs_coeff = inv(Mtx' * Mtx) * Mtx' * y;
else
    [Usvd, Ssvd, Vsvd] = svd(Mtx); 
    SV = diag(Ssvd);  
    Idx_sv = find(SV>=TH_sv, 1, 'last'); 
    Pseudo_invMTX = Vsvd(:,1:Idx_sv) * diag(1./SV(1:Idx_sv)) * Usvd(:,1:Idx_sv)';
    ANFs_coeff = Pseudo_invMTX * y;
end

%ANFcoeff_org = ANFs_coeff;
yest.all = Mtx_org * ANFs_coeff;
yest.no_transition = Mtx * ANFs_coeff;
nmse = mean((y - yest.no_transition).^2) / mean(y.^2);

%ANFs.all = ANFs_coeff;
ANFs.const = ANFs_coeff(1); ANFs_coeff = ANFs_coeff(2:end);
for mmm=1:Npdms,
    eval(['ANFs.pdm' num2str(mmm) ' = ANFs_coeff(1:ANFs_order);']); 
    ANFs_coeff = ANFs_coeff(ANFs_order+1:end);
end

return


        
        
        
        
        
        
        