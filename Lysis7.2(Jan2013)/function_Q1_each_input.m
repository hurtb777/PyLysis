function Vx = function_Q1_each_input(x, Laguerre)

N = length(x); L = size(Laguerre, 2); 

% first-order of one input
Vx = zeros(N, L);
for k=1:L,
    v_column = conv(Laguerre(:,k), x(:));   
    Vx(:,k) = v_column(1:N);
end

return











