function Kest = function_LET2Q2_kernelestimate(Cest, Laguerre1, Laguerre2, L1, L2, Q)    

M1 = max(size(Laguerre1)); M2 = max(size(Laguerre2)); 

% calculate the first, second, and third order kernels
Kest.k0 = Cest.c0; 
k10 = Laguerre1 * Cest.c10; Kest.k10 = k10;
k01 = Laguerre2 * Cest.c01; Kest.k01 = k01;

if Q==2, 
    c20 = Cest.c20; 
    temp = zeros(M1, M1); 
    for k1=1:L1
        for k2=1:L1,
            temp = temp +c20(k1,k2) * Laguerre1(:,k1) * Laguerre1(:,k2)';
        end
    end    
    Kest.k20 = temp; clear temp

    c02 = Cest.c02; 
    temp = zeros(M2, M2); 
    for k1=1:L2
        for k2=1:L2,
            temp = temp +c02(k1,k2) * Laguerre2(:,k1) * Laguerre2(:,k2)';
        end
    end    
    Kest.k02 = temp; clear temp

    fieldnamesCest = fieldnames(Cest);
    Xterms = 0; 
    for kk=1:length(fieldnamesCest),
        if strcmp(fieldnamesCest(kk), 'c11'), Xterms = 1; end
    end
    if Xterms==1, 
        c11 = Cest.c11;
        temp = zeros(M1, M2);  
        for m1=1:M1,
            for m2=1:M2,
                for j1=1:L1,
                    for j2=1:L2,
                        temp(m1,m2) = temp(m1,m2) + c11(j1,j2) * Laguerre1(m1,j1)*Laguerre2(m2,j2);
                    end
                end
            end
        end
        Kest.k11 = temp; clear temp
    end
end

return











