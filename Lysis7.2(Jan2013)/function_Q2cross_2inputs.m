function Vxx = function_Q2cross_2inputs(V1, V2)

N = size(V1, 1); L1 = size(V1, 2); L2 = size(V2, 2);

% cross 2nd-order terms
Nxx = L1 *L2; 
Vxx = zeros(N, Nxx);

cnt = 0; 
for k1=1:L1,
    for k2=1:L2,
        cnt = cnt + 1; 
        Vxx(:,cnt) = V1(:,k1) .* V2(:,k2);
    end
end

return











