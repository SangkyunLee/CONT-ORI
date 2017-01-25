function [X,mX,stdX,oris,cons] = sortdata(dtype,ises,exp_type)
% function [X,oris,cons] = sortdata(dtype,ises,exp_type)
% sort two separate session data by event 
% output X: cell variable of 2 x events

% exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE'};

fnpf1='P1';
fnpf2='P2';
ctm=0.6;
cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type,DATA_thr_str);

fndata1 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf1);
fndata2 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf2);
[D1, D2]= loadData(data_path,fndata1,fndata2);
subinx = intersect(D1.cellinx_sel,D2.cellinx_sel);

[D1,D2] = subdata(subinx,D1,D2,{dtype});




inxs_valtrial = D1.events(:)>0 & ~isinf(D1.events(:));
unique_evt = unique(D1.events(inxs_valtrial)');
unique_evt = setdiff(unique_evt,0);
nevt = length(unique_evt);
E1 = D1.events(:);
E2 = D2.events(:);
inxsample = cell(2,nevt);
cons = zeros(1,nevt);
oris = zeros(1,nevt);
X = cell(size(inxsample));
mX = zeros(nevt, 2,length(D1.cellinx_sel));
stdX = zeros(nevt,2, length(D1.cellinx_sel));
for ievt0 = 1:nevt
    ievt = unique_evt(ievt0);
    inxsample{1,ievt0} = TP.select_subdata(E1,{ievt});
    inxsample{2,ievt0} = TP.select_subdata(E2,{ievt});
    cons(ievt0)=D1.events_cont(inxsample{1,ievt0}(1));
    oris(ievt0)=D1.events_ORI(inxsample{1,ievt0}(1));
    
    X{1,ievt0}=D1.(dtype)(inxsample{1,ievt0},:);
    X{2,ievt0}=D2.(dtype)(inxsample{2,ievt0},:);    
    mX(ievt,1,:) = mean(X{1,ievt0},1);
    mX(ievt,2,:) = mean(X{2,ievt0},1);
    stdX(ievt,1,:) = std(X{1,ievt0},0,1);
    stdX(ievt,2,:) = std(X{2,ievt0},0,1);
end


