function Vxx = function_Q2self_each_input(V)

N = size(V, 1); L = size(V, 2); 

% self 2nd-order terms of one input
Nxx = L * (L + 1 )/ 2; 
Vxx = zeros(N, Nxx);

cnt = 0; 
for k1=1:L,
    for k2=k1:L,
        cnt = cnt + 1; 
        Vxx(:,cnt) = V(:,k1) .* V(:,k2);
    end
end

return











