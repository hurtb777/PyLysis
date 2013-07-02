function [uu1, uu2, uuX, index_pdm1, index_pdm2] = ...
        function_LET2PDMoutputs_self_and_Xterms(PDMs1, in_sig1, PDMs2, in_sig2, out_sig, TH)
                                                        
in_sig1 = in_sig1(:); in_sig2 = in_sig2(:); out_sig = out_sig(:); 

N = length(out_sig); 
Npdms1 = size(PDMs1, 2); Npdms2 = size(PDMs2, 2);                                                        
uu1 = zeros(N, Npdms1);
for kk=1:Npdms1,
    temp1 = conv(in_sig1, PDMs1(:,kk));
    uu1(:,kk) = temp1(1:N);
end

uu2 = zeros(N, Npdms2);
for kk=1:Npdms2,
    temp2 = conv(in_sig2, PDMs2(:,kk));
    uu2(:,kk) = temp2(1:N);
end

uuX = zeros(N, Npdms1*Npdms2); 
XcorrcoeffV = zeros(Npdms1*Npdms2,1); 
XcorrcoeffM = zeros(Npdms1, Npdms2); 
cnt = 0; 
for kk1=1:Npdms1,
    for kk2=1:Npdms2,
        temp = uu1(:,kk1) .* uu2(:,kk2);
        cnt = cnt + 1; 
        uuX(:,cnt) = temp;
        XcorrcoeffV(cnt) = xcorr(temp-mean(temp), out_sig-mean(out_sig), 0, 'biased');
        XcorrcoeffV(cnt) = XcorrcoeffV(cnt) / (std(temp) * std(out_sig));
        XcorrcoeffM(kk1,kk2) = XcorrcoeffV(cnt);
    end
end
clear temp
 
[XcorrcoeffV_order, index_order] = sort(abs(XcorrcoeffV), 'descend');
XcorrcoeffTH = find(XcorrcoeffV_order>=TH);
if isempty(XcorrcoeffTH), 
    index_pdm1 = []; index_pdm2 = []; uuX = [];
else
    Nxterms = length(XcorrcoeffTH);   
    index_order = index_order(1:Nxterms); 
    index_pdm1 = zeros(Nxterms, 1); index_pdm2 = zeros(Nxterms, 1);  
    for kk=1:Nxterms,
        index_pdm1(kk) = ceil(index_order(kk)/Npdms2);
        index_pdm2(kk) = rem(index_order(kk), Npdms2);
        if index_pdm2(kk)==0, index_pdm2(kk) = Npdms2; end
    end

    uuX = zeros(N, Nxterms);
    for kk=1:Nxterms,
        temp = uu1(:,index_pdm1(kk)) .* uu2(:,index_pdm2(kk));
        uuX(:,kk) = temp;
        XcorrcoeffV = xcorr(temp-mean(temp), out_sig-mean(out_sig), 0, 'biased');
        XcorrcoeffV = XcorrcoeffV / (std(temp) * std(out_sig));
        if XcorrcoeffV_order(kk)~=abs(XcorrcoeffV), error('?? wrong selection ??'); end
    end 
end

return
    

        


