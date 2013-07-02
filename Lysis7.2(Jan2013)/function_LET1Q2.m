function Kest = function_LET1Q2(x, y, alpha, L, M, Nfig)

x = x(:); y = y(:);
N = length(x);

 % LET parameters
Mdiscard = M; Q = 2; 

LETparameters.alpha = alpha; 
LETparameters.Q = Q;
LETparameters.L = L;
LETparameters.M = M;
LETparameters.Mdiscard = Mdiscard;

[NMSE, Cest, Kest, prediction_LET1Q2] = function_LET1Q2_LFparameters(x, y, LETparameters);             
            
figure(Nfig), clf; subplot(211), plot((0:M-1), Kest.k1, 'linewidth', 2); 
        title(['LET1 analysis: \alpha = ', num2str(alpha), ', L = ', num2str(L), ' and Q = 2']); 
        grid; ylabel('1^{st}-order kernel'); 
    subplot(212), mesh((0:M-1), (0:M-1), Kest.k2); colorbar; 
        zlabel('2^{nd}-order kernel'); set(gca, 'xlim', [0, M]); set(gca, 'ylim', [0, M]); drawnow  

figure(Nfig+1), clf; plot((0:N-1), y, 'b', (0:N-1), prediction_LET1Q2, 'r', 'linewidth', 2); 
        grid; legend('output','estimated output'); set(gca, 'xlim', [Mdiscard+1 N-1]); 
        title(['LET1 analysis: \alpha = ', num2str(alpha), ...
            ', L = ', num2str(L), ' and Q = 2']);
        xlabel(['NMSE = ', num2str(NMSE)]); drawnow 

return

    


    