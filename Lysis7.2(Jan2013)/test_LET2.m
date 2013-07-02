clear; clc; close all;

% length of data samples
N = 4096;
M = 50;

% two inputs
x1 = randn(N, 1); x2 = randn(N, 1);
DLF1 = function_generate_laguerre(0.4, 3, M);
DLF2 = function_generate_laguerre(0.6, 3, M);
    h1 = 1 * DLF1(:,1) + 0.5 * DLF1(:,2) - 0.5 * DLF1(:,3);
    h2 = 0.5 * DLF2(:,1) + 0.5 * DLF2(:,2) - 1.5 * DLF2(:,3);
    v1 = conv(x1, h1); v1 = v1(1:N);
    v2 = conv(x2, h2); v2 = v2(1:N);
    y=v1+v2.^2;
    x1 = x1 - mean(x1); x2 = x2 - mean(x2); y = y - mean(y);
    figure(1), clf; subplot(311), plot(y, 'linewidth', 2); grid; ylabel('output');
        set(gca, 'xlim', [0 N]);
    subplot(312), plot(h1, 'linewidth', 2); grid; ylabel('h1');
    subplot(313), plot(h2, 'linewidth', 2); grid; ylabel('h2');
        set(gca, 'xlim', [0 N]); drawnow

% Laguerre parameters
L1 = 5; L2 = 5; 
alpha1 = 0.5; alpha2 = 0.5; 

% length of Laguerre function (M1 and M2) and the length of transition phase (Mdiscard)
Q = 2; Nfig = 1; 

[Cest, Kest, Pred, NMSE] = LET_2(x1, x2, y, alpha1, alpha2, L1, L2, Q, Nfig);
    