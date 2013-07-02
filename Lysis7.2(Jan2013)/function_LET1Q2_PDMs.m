function [PDMs, Npdms] = function_LET1Q2_PDMs(x, Kest, L, Nfig, inputID)

rms.x1 = mean((x-mean(x)).^2); 
Mtx_kernels = [Kest.k1 rms.x1*Kest.k2];
[Uk, Sk, Vk] = svd(Mtx_kernels);
for kk=1:size(Uk,2);
    temp = Uk(:,kk);
    if sum(temp)<0, Uk(:,kk) = -Uk(:,kk); end
end
Sk = diag(Sk); 
figure(Nfig), clf; subplot(211), stem(Sk, '.', 'linewidth', 2); grid;
    if isempty(inputID), 
        title('singular values (S_k) of the kernel matrix'); 
    else
        title(['singular values (S_k) of the kernel matrix ', inputID]); 
    end
    ylabel('S_k'); set(gca, 'Xlim', [0.5 L+1+0.5]); set(gca, 'Xtick', (1:L+1)); drawnow
        
MM = length(Sk); 
sum_of_Sk = cumsum(Sk)/sum(Sk);
subplot(212), plot((1:MM), sum_of_Sk, 'b', [1 MM], [0.9, 0.9], 'r--', 'linewidth', 3); grid; 
    ylabel('cumulative S_k'); 
    axis([0.5 L+0.5 0.7 1.1]); set(gca, 'Xtick', (1:L+1)); drawnow

% threshold for number of PDMs
Npdms = input(['Enter a number of PDMs ' inputID ' ( <= ', num2str(L) '): ']); 
if isempty(Npdms), 
    Npdms = 1; 
end
if Npdms>L, 
    Npdms = L; 
end

PDMs = Uk(:,1:Npdms); 
figure(Nfig+1), clf; 
for kk=1:Npdms,
    subplot(Npdms,1,kk), plot((0:length(Kest.k1)-1), PDMs(:,kk), 'linewidth', 2); grid;
        if kk==1, title(['PDMs ', inputID]); end 
        ylabel(['PDM #', num2str(kk)]); set(gca, 'Xlim', [0 length(Kest.k1)-1]); drawnow
end

return

    