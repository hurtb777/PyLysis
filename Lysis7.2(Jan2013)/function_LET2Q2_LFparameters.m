function [Cest, Kest, y_predict, nmse] = function_LET2Q2_LFparameters(x1, x2, y, LETparameters)

% parameters for LET function
A1 = LETparameters.alpha1; L1 = LETparameters.L1; M1 = LETparameters.M1;
A2 = LETparameters.alpha2; L2 = LETparameters.L2; M2 = LETparameters.M2;
Mdiscard = LETparameters.Mdiscard; 
Q = LETparameters.Q;

% estimation of Laguerre coefficients and Volterra kernels 
Laguerre1 = function_generate_laguerre(A1, L1, M1);
Laguerre2 = function_generate_laguerre(A2, L2, M2);

% estimate the Laguerre coefficients
V0 = ones(length(x1), 1); 
V1 = function_Q1_each_input(x1, Laguerre1); 
V2 = function_Q1_each_input(x2, Laguerre2);
VV_est = [V0 V1 V2]; 

if Q==2, 
    V1_2nd = function_Q2self_each_input(V1);
    V2_2nd = function_Q2self_each_input(V2);
    V12_2nd = function_Q2cross_2inputs(V1, V2); 

    VV_est = [VV_est V1_2nd V2_2nd V12_2nd];
end

% estimate the Laguerre coefficients
[Cest, y_predict, nmse] = function_LET2Q2_coeffestimate(y, VV_est, L1, L2, Q, Mdiscard);
 
% estimate the Volterra kernel from the Laguerre coefficient and functions
Kest = function_LET2Q2_kernelestimate(Cest, Laguerre1, Laguerre2, L1, L2, Q);

return


    