clear; clc; close all;

% length of data samples
N = 1024;

% single input
x = randn(N, 1); 

% first-order kernel for testing
DLF = function_generate_laguerre(0.5, 3, 40);
h = -0.5 * DLF(:,1) + 1 * DLF(:,2) - 1.5 * DLF(:,3);

% output
v = conv(x, h); y = v(1:N)+ v(1:N).^2;

% Laguerre parameters
alpha = 0.5;

% length of Laguerre function (M) and the length of transition phase (Mdiscard)
L = 3; Q = 2; Nfig = 1;

[Cest, Kest, Pred, NMSE] = LET_1(x, y, alpha, L, Q, Nfig);



    