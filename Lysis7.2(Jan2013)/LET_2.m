function [Cest, Kest, pred, NMSE] = LET_2(x1, x2, y, alpha1, alpha2, L1, L2, Q, Nfig)
% [Cest, Kest, pred, NMSE] = LET_2(x1, x2, y, alpha1, alpha2, L1, L2, Q, Nfig)
% input variables:
%   x1: first input-data vector
%   x2: second input-data vector
%   y: output-data vector
%   alpha1: parameter alpha of Laguerre functions for x1 (default value: 0.5)
%   alpha2: parameter alpha of Laguerre functions for x2 (default value: 0.5)
%   L1: number of Laguerre functions for x1 (default value: 5)
%   L2: number of Laguerre functions for x2 (default value: 5)
%   Q: order of the model (Q = 1 or 2, default value: 2)
%   Nfig: an integer value for the first plotting window
%
% output variables:
%   Cest: Laguerre coefficient estimates (structure array)
%           Cest.c0 - Laguerre constant
%           Cest.c10 - first-order Laguerre coefficient estimate for x1
%           Cest.c20 - second-order Laguerre coefficient estimate for x1
%           Cest.c01 - first-order Laguerre coefficient estimate for x2
%           Cest.c02 - second-order Laguerre coefficient estimate for x2
%           Cest.c11 - second-order cross-Laguerre coefficient estimate for x1 and x2
%   Kest: kernel estimates (structure array)
%           Kest.k0 - constant
%           Kest.k10 - first-order kernel estimate for x1
%           Kest.k20 - second-order kernel estimates for x1
%           Kest.k01 - first-order kernel estimate for x2
%           Kest.k02 - second-order kernel estimates for x2
%           Kest.k11 - second-order cross-kernel estimate for x1 and x2
%   pred: model prediction
%   NMSE: normalized mean-square error of model prediction

x1 = x1(:); x2 = x2(:); y = y(:);
if length(x1)~= length(x2) || length(x1)~= length(y) || length(x2)~= length(y), 
    error('lengths of inputs and output are different'); 
end

N = length(x1); 

if alpha1<=0 || alpha1>=1, 
    disp('alpha1 must be between 0 and 1'); 
    alpha1 = 0.5; 
end

M1 = (-30 - log(1-alpha1)) / log(alpha1); 
M1 = ceil(M1);

if alpha2<=0 || alpha2>=1, 
    disp('alpha2 must be between 0 and 1'); 
    alpha2 = 0.5; 
end
M2 = (-30 - log(1-alpha2)) / log(alpha2); 
M2 = ceil(M2);

if L1<=0 || L1>9, 
    disp('L1 should be between 1 and 9, default vlaue L1 = 5'); 
    L1 = 5; 
end
L1 = round(L1); 
if L2<=0 || L2>9, 
    disp('L2 should be between 1 and 9, default vlaue L2 = 5'); 
    L2 = 5; 
end
L2 = round(L2); 

if Q~=1 && Q~=2, 
    Q = 2; 
end
if isempty(Nfig), 
    Nfig = 1; 
end
Mdiscard = max(M1, M2);

LETparameters.alpha1 = alpha1; LETparameters.alpha2 = alpha2; 
LETparameters.L1 = L1; LETparameters.L2 = L2; 
LETparameters.M1 = M1; LETparameters.M2 = M2;
LETparameters.Q = Q;
LETparameters.Mdiscard = Mdiscard;

[Cest, Kest, pred, NMSE] = function_LET2Q2_LFparameters(x1, x2, y, LETparameters);

figure(Nfig), clf; subplot(211), plot((0:M1-1), Kest.k10, 'linewidth', 2); 
    title(['LET for first input: \alpha1 = ', num2str(alpha1), ...
              ', L_1 = ' num2str(L1) ' and Q = ', num2str(Q)]); grid; ylabel('k_1');
    set(gca, 'xlim', [0 M1]); drawnow
    if Q==2, 
        subplot(212), mesh((0:M1-1), (0:M1-1), Kest.k20); colorbar
        zlabel('k_2'); set(gca, 'xlim', [0 M1]); set(gca, 'ylim', [0 M1]); drawnow  
    end
figure(Nfig+1), clf; subplot(211), plot((0:M2-1), Kest.k01, 'linewidth', 2); 
    title(['LET for second input: \alpha2 = ', num2str(alpha2), ...
              ', L_2 = ' num2str(L2) ' and Q = ', num2str(Q)]); grid; ylabel('k_1'); 
    set(gca, 'xlim', [0 M2]); drawnow      
    if Q==2, 
        subplot(212), mesh((0:M2-1), (0:M2-1), Kest.k02); colorbar
        zlabel('k_2'); set(gca, 'xlim', [0 M2]); set(gca, 'ylim', [0 M2]); drawnow  
    end
figure(Nfig+2), clf; 
if Q==2, 
    figure(Nfig+2), clf; mesh((0:M1-1), (0:M2-1), Kest.k11.'); colorbar
        ylabel('k_{11}'); set(gca, 'xlim', [0 M1]); set(gca, 'ylim', [0 M2]);  
        title(['LET analysis: \alpha1 = ', num2str(alpha1), ' & \alpha2 = ', ...
            num2str(alpha2) '; L1 = ', num2str(L1), ' & L2 = ', num2str(L2) ...
            '; and Q = ', num2str(Q)]); drawnow 
end    
figure(Nfig+3), clf; subplot(211),
    plot((0:N-1), y, 'b', (0:N-1), pred, 'r', 'linewidth', 2); 
        grid; ylabel('y (blue) and model prediction (red)'); set(gca, 'xlim', [Mdiscard N]);  
        title(['LET2 analysis: \alpha1 = ', num2str(alpha1), ' and \alpha2 = ', ...
            num2str(alpha2) '; L1 = ', num2str(L1), ' & L2 = ', num2str(L2) ...
            '; and Q = ', num2str(Q)]);
    subplot(212), plot((0:N-1), y-pred, 'b', 'linewidth', 2); 
        grid; ylabel('errors'); set(gca, 'xlim', [Mdiscard N]); 
        xlabel(['NMSE = ', num2str(NMSE)]); drawnow

    