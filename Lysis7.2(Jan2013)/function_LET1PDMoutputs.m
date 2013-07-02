function uu = function_LET1PDMoutputs(x, PDMs)
                                                        
x = x(:);
N = length(x); 

Npdms = size(PDMs, 2);                                                      
uu = zeros(N, Npdms);
for kk=1:Npdms,
    temp = conv(x, PDMs(:,kk));
    uu(:,kk) = temp(1:N);
end

return
    

        


