function Laguerre = function_generate_laguerre(alpha, L, M)
% Generate the Laguerre functions: Laguerre = generate_laguerre(alpha, L, M);
%                  order of Laguerre function: 0, 1, 2, ..., L
%                  sample length of each Lagureer function: 0, 1, 2, ..., M-1

mm = (0:M-1);
beta = 1 - alpha;   rootAlpha = sqrt(alpha);
L_buf = sqrt((alpha .^ mm.') * beta);

mmm = zeros(M,1); 
Laguerre(:,1) = L_buf;

for n1=1:L-1,
    for n2=mm+1,
        mmm(n2) = L_buf(n2);
        if n2==1, L_buf(n2) = rootAlpha * L_buf(n2);
        else      L_buf(n2) = rootAlpha * (L_buf(n2-1) + mmm(n2)) - mmm(n2-1);
        end
    end
    Laguerre(:,n1+1) = L_buf;
end

return