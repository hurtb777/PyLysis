function [Cest, Kest, pred, NMSE] = LET_1(x, y, alpha, L, Q, Nfig)
% [Cest, Kest, pred, NMSE] = LET_1(x, y, L, M, Q, Nfig)
% input variables:
%   x: input-data vector
%   y: output-data vector
%   alpha: alpha parameter of Laguerre functions (from 0 to 1, default value: 0.5)
%   L: number of Laguerre functions used (from 1 to 9, default value: 5)
%   Q: order of the model : 1 for first-order or 2 for second-order (default value: 2)
%   Nfig: an integer value for the first plotting window
%
% output variables:
%   Cest: Laguerre coefficient estimates (structure array)
%           Cest.c0 - constant
%           Cest.c1 - first-order Laguerre coefficient estimates
%           Cest.c2 - second-order Laguerre coefficient estimates
%   Kest: kernel estimates (structure array)
%           Kest.k0 - constant
%           Kest.k1 - first-order kernel estimate
%           Kest.k2 - second-order kernel estimate
%   pred: model prediction using the kernel estimates
%   NMSE: normalized mean-square error of model prediction

x = x(:); y = y(:);
if length(x)~= length(y), 
    error('lengths of input and output are different'); 
end

N = length(x); 

if alpha<=0 || alpha>=1, 
    disp('alpha must be between 0 and 1'); 
    alpha = 0.5; 
end
M = (-30 - log(1-alpha)) / log(alpha); 
M = ceil(M);

if L<=0 || L>9, 
    disp('L should be between 1 and 9, default vlaue L = 5'); 
    L = 5; 
end
L = round(L); 

if Q~=1 && Q~=2, 
    Q = 2; 
end
if isempty(Nfig), 
    Nfig = 1; 
end
Mdiscard = M;

LETparameters.alpha = alpha; 
LETparameters.M = M;
LETparameters.Mdiscard = Mdiscard;
LETparameters.Q = Q;
LETparameters.L = L;

%% LET1 Q1 or Q2
% first-order and second-order kernels for the input: Kest.k1 and Kest.k2
[Cest, Kest, pred, NMSE] = function_LET1_LFparameters(x, y, LETparameters);

%%
if ~isempty(Nfig), 
    figure(Nfig), clf; subplot(211), plot((0:M-1), Kest.k1, 'linewidth', 2); 
        title(['LET1: \alpha = ', num2str(alpha), ', L = ', num2str(L), ...
            ' and Q = ', num2str(Q)]); grid; ylabel('k_1'); 
        if Q==2, 
            subplot(212), mesh((0:M-1), (0:M-1), Kest.k2); colorbar
            zlabel('k_2'); drawnow  
        end
    figure(Nfig+1), clf; subplot(211), 
        plot((0:N-1), y, 'b', (0:N-1), pred, 'r', 'linewidth', 2); 
            grid; ylabel('y (blue) and model prediction (red)'); set(gca, 'xlim', [0 N]);  
            title(['LET1 analysis: \alpha = ', num2str(alpha), ...
                            ', L = ', num2str(L), ' and Q = ', num2str(Q)]);
        subplot(212), plot((0:N-1), y-pred, 'b', 'linewidth', 2); 
            grid; ylabel('errors'); set(gca, 'xlim', [0 N]); 
            xlabel(['NMSE = ', num2str(NMSE)]); drawnow  
end   


    