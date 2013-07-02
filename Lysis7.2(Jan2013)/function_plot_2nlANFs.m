function function_plot_2nlANFs(inputID, channelID, uu, PDMs, ANFs_order, nlANF_selfXterms, wgt_domain, Nfig) 

Npdms = size(PDMs, 2); M = size(PDMs, 1); 

for mmm=1:Npdms,
    eval(['NLcoeff = nlANF_selfXterms.i' num2str(inputID) 'pdm' num2str(mmm), ';']);
    uu_bound = ceil(wgt_domain*std(uu(:,mmm)));
    min_domain = -uu_bound; max_domain = uu_bound;
    in_domain = (min_domain:0.01:max_domain)+mean(uu(:,mmm));
    out_range = zeros(size(in_domain));
    for kk=1:ANFs_order,
        out_range = out_range + NLcoeff(kk) * in_domain.^kk;
    end
    eval(['domain.pdm' num2str(mmm) ' = in_domain;']);
    eval(['range.pdm' num2str(mmm) ' = out_range;']);
end    

xmin = zeros(Npdms, 1); xmax = zeros(Npdms, 1); 
figure(Nfig), clf;
for mmm=1:Npdms,
    eval(['x_domain = domain.pdm' num2str(mmm) ';']);
    xmin(mmm) = min(x_domain); xmax(mmm) = max(x_domain);
    eval(['y_range = range.pdm' num2str(mmm) ';']);
    if mmm==1,     plot(x_domain, y_range, 'b', 'linewidth', 3); 
    elseif mmm==2, plot(x_domain, y_range, 'g', 'linewidth', 3); 
    elseif mmm==3, plot(x_domain, y_range, 'r', 'linewidth', 3); 
    elseif mmm==4, plot(x_domain, y_range, 'k', 'linewidth', 3); 
    elseif mmm==5, plot(x_domain, y_range, 'm', 'linewidth', 3); 
    elseif mmm==6, plot(x_domain, y_range, 'c', 'linewidth', 3); 
    elseif mmm==7, plot(x_domain, y_range, 'y', 'linewidth', 3); 
    elseif mmm==8, plot(x_domain, y_range, 'b--', 'linewidth', 3);
    elseif mmm==9, plot(x_domain, y_range, 'g--', 'linewidth', 3);
    end
    hold on;
end
hold off; grid; title(['cubic ANFs of PDMs for ' channelID]) 
set(gca, 'Xlim', [min(xmin) max(xmax)]);  
warning off;
legend('PDM #1', 'PDM#2', 'PDM #3', 'PDM #4', 'PDM #5', 'PDM #6', 'PDM #7', 'PDM #8', 'PDM #9'); 
xlabel('u'); ylabel('z'); drawnow
warning on;
            
return

