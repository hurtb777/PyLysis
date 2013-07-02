function [Npdms, PDMs, ANFs, pred, NMSE] = PDM_1(x, y, alpha, L, Nfig)
% [Npdms, PDMs, ANFs, pred, NMSE] = PDM_1(x, y, alpha, L, Nfig)
% input variables:
%   x: input-data vector
%   y: output-data vector
%   alpha: alpha value of Laguerre functions (default value: 0.5)
%   L: number of Laguerre functions (default value: 5)
%   Nfig: an integer value for plotting windows
%
% output variables:
%   Npdms: number of PDMs
%   PDMs: estimates of Principal Dynamic Modes
%   ANFs: cubic polynomial function (with zero constant term) for each PDM
%           ANFs.const - constant term of ANFs
%           ANFs.pdmX - ANFs for PDM #X (first, second, and third order coefficients)
%   pred: PDM-based model prediction
%   NMSE: normalized mean-squre error of PDM-based model prediction

x = x(:); y = y(:);
if length(x)~= length(y), 
    error('lengths of input and output are different'); 
end

N = length(x); 

if alpha<=0 || alpha>=1, 
    disp('alpha should be between 1 and 9, default vlaue alpha = 0.5'); 
    alpha = 0.5; 
end
M = (-30 - log(1-alpha)) / log(alpha); 
M = ceil(M);

if L<=0 || L>9, 
    disp('L should be between 1 and 9, default vlaue L = 5'); 
    L = 5; 
end

if isempty(Nfig), 
    Nfig = 1; 
end

%% LET1Q2: estimate first-order and second-order kernels from "alpha" and "M"        
Kest = function_LET1Q2(x, y, alpha, L, M, Nfig); 

%% estimate PDMs from the first-order and second-order kernel estimates 
[PDMs, Npdms] = function_LET1Q2_PDMs(x, Kest, L, Nfig+1, []);     

[pred, NMSE, ANFs, uu] = function_estimate_3rd_order_ANFs(x, y, PDMs); 
function_plot_nlANFs(uu, PDMs, ANFs, Nfig+3); 

figure(Nfig+4), clf; subplot(211), plot((0:N-1), y, 'b', (0:N-1), pred.all, 'r', 'linewidth', 2); 
        grid; legend('output signal','model prediction'); set(gca, 'xlim', [M+1, N-1]); 
        title(['PDM and ANF analysis: ' num2str(Npdms) ' PDMs']); drawnow 
    subplot(212), plot((0:N-1), y-pred.all, 'linewidth', 2); 
        grid; set(gca, 'xlim', [M+1, N-1]); ylabel('residual'); 
        xlabel(['nmse = ', num2str(NMSE)]); drawnow 

pred = pred.all; 

return



