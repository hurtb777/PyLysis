function [nmse, Cest, Kest, y_predict] = function_LET1Q2_LFparameters(x, y, LETparameters)

% parameters for LET function
alpha = LETparameters.alpha; L = LETparameters.L; M = LETparameters.M;
Q = LETparameters.Q; Mdiscard = LETparameters.Mdiscard; 

x = x(:); y = y(:);

% estimation of Laguerre coefficients and Volterra kernels 
Laguerre = function_generate_laguerre(alpha, L, M);

% estimate the Laguerre coefficients
V0 = ones(length(x), 1); 
V1 = function_Q1_each_input(x, Laguerre); 
V2 = function_Q2self_each_input(V1);

VV_est = [V0 V1 V2];

% estimate the Laguerre coefficients
[Cest, y_predict, nmse] = function_LET1Q2_coeffestimate(y, VV_est, L, Q, Mdiscard);
 
% estimate the Volterra kernel from the Laguerre coefficient and functions
Kest = function_LET1Q2_kernelestimate(Cest, Laguerre, L, Q);

return


        
        
        
        
        
        
        