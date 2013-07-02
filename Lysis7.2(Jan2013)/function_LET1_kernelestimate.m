function Kest = function_LET1_kernelestimate(Cest, Laguerre, L, Q)    

M = max(size(Laguerre));

% calculate the first kernel
Kest.k0 = Cest.c0; 
k1 = Laguerre * Cest.c1; Kest.k1 = k1;

if Q==2, 
    % calculate the second kernel
    c2 = Cest.c2; 
    temp = zeros(M, M); 
    for k1=1:L
        for k2=1:L,
            temp = temp +c2(k1,k2) * Laguerre(:,k1) * Laguerre(:,k2)';
        end
    end    
    Kest.k2 = temp; clear temp
end

return











