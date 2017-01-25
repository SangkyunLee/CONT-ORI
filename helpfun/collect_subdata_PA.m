function [data, evt,sepinx]=collect_subdata_PA(data, evt, sevts, eORI, PAthr)
% function [data, evt]=collect_subdata_PA(data, evt, sevts, eORI, PAthr)
% select data based on population activity with PA thresholds

     
    inxsample = TP.select_subdata(evt,sevts);
    eORI = eORI(inxsample);
    D = data(inxsample,:);

    % similar number of samples between two conditions,not exceed 10% more samples               
    [~, inxsep] = TP.select_samples_evtratio(eORI,sevts{2}, 1.1);

    P=cell(1,length(inxsep));
    for ix = 1 : length(inxsep)                
        P{ix} = mean(D(inxsep{ix},:),2);
    end    
    
    sepinx0 = cellfun(@sepdata,P,repmat({PAthr},size(P)),'UniformOutput', false);    
    if length(inxsep) ==1
        sepinx0 = inxsep{1}(sepinx0{1}{1})'; 
    elseif length(inxsep) ==2
        sepinx0 = [inxsep{1}(sepinx0{1}{1})' inxsep{2}(sepinx0{2}{1})'] ;
    end
    


    data = D(sepinx0,:)';
    evt = eORI(sepinx0);
    sepinx = inxsample(sepinx0);
    
    %----------------------------------
end

function y = sepdata(x,thr)
% function y = sepdata(x,thr)
% x : 1-row/comlumn vector
% thr: 2xnseg matrix
y = cell(1,size(thr,2));
[~,inx]=sort(x,'ascend');

n = length(x);

for i = 1 : size(thr,2)
    i1 = floor(thr(1,i)*n)+1;    
    i2 = floor(thr(2,i)*n);
    y{i}=inx(i1:i2);
end
end
    
    