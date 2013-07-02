function [Npdms1, PDMs1, Npdms2, PDMs2, ANFs, pred, NMSE] = ...
                                            PDM_2(x1, x2, y, alpha1, alpha2, L1, L2, Nfig)
% [Npdms1, PDMs1, Npdms2, PDMs2, ANFs, pred, NMSE] = PDM_2(x1, x2, y, alpha1, alpha2, L1, L2, Nfig)
% input variables:
%   x1: first input-data vector
%   x2: second input-data vector
%   y: output-data vector
%   alpha1: parameter of Laguerre functions for the first input  (default value: 0.5)
%   alpha2: parameter of Laguerre functions for the second input   (default value: 0.5)
%   L1: number of Laguerre functions for the first input (default value: 5)
%   L2: number of Laguerre functions for the second input (default value: 5)
%   Nfig: an integer for the first plotting window
%
% output variables:
%   Npdms1: number of PDMs for the first input
%   PDMs1: PDM estimates for the first input
%   Npdms2: number of PDMs for the second input
%   PDMs2: PDM estimates for the second input
%   ANFs: cubic polynomials for each PDM estimate(structure array)
%           ANFs.const - constant term of ANFs
%           ANFs.iXpdmY - cubic ANF coefficients for input X and PDM Y (zero constant)
%           ANFs.cross_terms - coefficients of the significant cross-terms, 
%                              which are pair-products of PDM-outputs selected 
%                              for statistical significance by use of the w-statistic
%   pred: PDM-based model prediction
%   NMSE: normalized mean-squre error of PDM-based model prediction

x1 = x1(:); x2 = x2(:); y = y(:);
if length(x1)~= length(y) || length(x2)~= length(y), 
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

if isempty(Nfig), 
    Nfig = 1; 
end
Mdiscard = max(M1, M2);

%% LET1Q2: estimate first-order and second-order kernels from "alpha" and "M"        
LETparameters.alpha1 = alpha1; LETparameters.alpha2 = alpha2; 
LETparameters.L1 = L1; LETparameters.L2 = L2; 
LETparameters.M1 = M1; LETparameters.M2 = M2;
LETparameters.Q = 2;
LETparameters.Mdiscard = Mdiscard;
[Cest, Kest] = function_LET2Q2_LFparameters(x1, x2, y, LETparameters);

%% estimate PDMs from the first-order and second-order kernel estimates 
Kest_input1.k1 = Kest.k10; Kest_input1.k2 = Kest.k20; 
[PDMs1, Npdms1] = function_LET1Q2_PDMs(x1, Kest_input1, L1, Nfig, 'for first input');     

Kest_input2.k1 = Kest.k01; Kest_input2.k2 = Kest.k02; 
[PDMs2, Npdms2] = function_LET1Q2_PDMs(x2, Kest_input2, L2, Nfig+2, 'for second input');    

%% estimate the ANFs and cross-terms
THcrossterms = (exp(3/sqrt(N))-1) / (exp(3/sqrt(N))+1); 
ANFs1_order = 3; ANFs2_order = 3; 
[ANFs, yprt, NMSE, ANFs_org] = function_LET2PDMs_estimate_3rdANFs_Xterms(PDMs1, x1, ...
                            PDMs2, x2, y, THcrossterms, ANFs1_order, ANFs2_order, Mdiscard); 

pred = yprt.all; 
function_plot_2nlANFs(1, 'first input', ANFs_org.uu1, PDMs1, ANFs1_order, ANFs_org, 1, Nfig+4);
function_plot_2nlANFs(2, 'second input', ANFs_org.uu2, PDMs2, ANFs2_order, ANFs_org, 1, Nfig+5);

figure(Nfig+6), clf; subplot(211), plot((0:N-1), y, 'b', (0:N-1), pred, 'r', 'linewidth', 2); 
        grid; ylabel('y (blue) and model prediction (red)');set(gca, 'xlim', [Mdiscard+1, N-1]); 
        title(['PDM analysis: ' num2str(Npdms1) ' PDMs for first input and ', ...
                    num2str(Npdms2), ' PDMs for second input']);
    subplot(212), plot((0:N-1), y-pred, 'linewidth', 2); set(gca, 'xlim', [Mdiscard+1, N-1]); 
        grid; ylabel('residual'); xlabel(['nmse = ', num2str(NMSE)]); drawnow 
        
return


