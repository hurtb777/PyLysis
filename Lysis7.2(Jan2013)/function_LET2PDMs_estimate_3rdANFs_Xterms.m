function [ANFs, yorg_predict, nmse, ANFs_org] = ...
    function_LET2PDMs_estimate_3rdANFs_Xterms(PDMs1, x1, PDMs2, x2, y, ...
                        THvalue, ANFs1_order, ANFs2_order, Mdiscard)

x1 = x1(:); x2 = x2(:); y = y(:); N = length(y);

[uu1, uu2, uuX, index_pdm1, index_pdm2] = ...
        function_LET2PDMoutputs_self_and_Xterms(PDMs1, x1, PDMs2, x2, y, THvalue);
    
ANFs_org.uu1 = uu1; ANFs_org.uu2 = uu2; ANFs_org.uuX = uuX; 
ANFs_org.index_pdm1 = index_pdm1; ANFs_org.index_pdm2 = index_pdm2; 

    
% ANFs for channel 1
Npdms1 = size(PDMs1, 2); 

ANFs_coeff1 = zeros(ANFs1_order * Npdms1, 1); 
Nunknowns1 = size(ANFs_coeff1, 1);

Mtx1 = zeros(N, Nunknowns1); 
cnt = 0; 
for mmm=1:Npdms1,
    for kk=1:ANFs1_order,
        cnt = cnt + 1;
        Mtx1(:,cnt) = uu1(:,mmm).^kk;
    end
end

% ANFs for channel 2
Npdms2 = size(PDMs2, 2); 

ANFs_coeff2 = zeros(ANFs2_order * Npdms2, 1); 
Nunknowns2 = size(ANFs_coeff2, 1);

Mtx2 = zeros(N, Nunknowns2); 
cnt = 0; 
for mmm=1:Npdms2,
    for kk=1:ANFs2_order,
        cnt = cnt + 1;
        Mtx2(:,cnt) = uu2(:,mmm).^kk;
    end
end

Mtx = [ones(N, 1), Mtx1, Mtx2, uuX]; 

Mtx_org = Mtx;
Mtx = Mtx(Mdiscard+1:N,:);
y = y(Mdiscard+1:N);

[EV, ED] = eig(Mtx' * Mtx);
min(abs(diag(ED)));
TH_sv = 1e-3; 
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
yorg_predict.all = Mtx_org * ANFs_coeff;
yorg_predict.input1 = Mtx_org(:,2:size(Mtx1,2)+1) * ANFs_coeff(2:size(Mtx1,2)+1);
yorg_predict.input2 = Mtx_org(:,size(Mtx1,2)+2:size(Mtx1,2)+size(Mtx2,2)+1) ...
                                * ANFs_coeff(size(Mtx1,2)+2:size(Mtx1,2)+size(Mtx2,2)+1);
yorg_predict.xterms = Mtx_org(:,size(Mtx1,2)+size(Mtx2,2)+2:end) ...
                        * ANFs_coeff(size(Mtx1,2)+size(Mtx2,2)+2:end);
y_predict = Mtx * ANFs_coeff;
nmse = mean((y - y_predict).^2) / mean(y.^2);

ANFs_org.all = ANFs_coeff;
ANFs_org.const = ANFs_coeff(1); ANFs.const = ANFs_coeff(1); 
ANFs_coeff = ANFs_coeff(2:end);
for mmm=1:Npdms1,
    eval(['ANFs_org.i1pdm' num2str(mmm) ' = ANFs_coeff(1:ANFs1_order);']); 
    eval(['ANFs.i1pdm' num2str(mmm) ' = ANFs_coeff(1:ANFs1_order);']); 
    ANFs_coeff = ANFs_coeff(ANFs1_order+1:end);
end
for mmm=1:Npdms2,
    eval(['ANFs_org.i2pdm' num2str(mmm) ' = ANFs_coeff(1:ANFs2_order);']); 
    eval(['ANFs.i2pdm' num2str(mmm) ' = ANFs_coeff(1:ANFs2_order);']);
    ANFs_coeff = ANFs_coeff(ANFs2_order+1:end);
end

ANFs_org.cross_terms = ANFs_coeff;
ANFs.cross_terms = ANFs_coeff;

return


        
        
        
        
        
        
        