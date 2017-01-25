clear all
close all
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE'};
datainxstr = {'AN1-16','AN17-22','','','AW23-40'};


iexp_type=5
ises=29



fnpf1='P1'
fnpf2='P2'
ctm=0.6
cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);

fndata1 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf1);
fndata2 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf2);
[D1, D2]= loadData(data_path,fndata1,fndata2);
subinx = intersect(D1.cellinx_sel,D2.cellinx_sel);
[D1,D2] = subdata(subinx,D1,D2,{'Xsel'});

switch iexp_type
    case {1}
        P1 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P1-AN1-16_DECSCR_CELLGRP_ctm0.60.mat')
        P1 = P1.S(ises);
        P2 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P2-AN1-16_DECSCR_CELLGRP_ctm0.60.mat')
        P2 = P2.S(ises);
    case{2}
        P1 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P1-AN17-22_DECSCR_CELLGRP_ctm0.60.mat')
        P1 = P1.S(ises);
        P2 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P2-AN17-22_DECSCR_CELLGRP_ctm0.60.mat')
        P2 = P2.S(ises);
    case{4}
    case{5}
        P1 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AWAKE\thr5\SUMMARY_P1-AW25-40_DECSCR_CELLGRP_ctm0.60.mat')
        P1 = P1.S(ises);
        P2 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AWAKE\thr5\SUMMARY_P2-AW25-40_DECSCR_CELLGRP_ctm0.60.mat')
        P2 = P2.S(ises);    
end

[contrasts, ORI_list, ORI_compindexset, ~, seslist] =get_expinfo(iexp_type);




evt1 = [D1.events_cont(:) D1.events_ORI(:)];
inxT1 = ~isinf(evt1(:,1));
evt1 = evt1(inxT1,:);

evt2 = [D2.events_cont(:) D2.events_ORI(:)];
inxT2 = ~isinf(evt2(:,1));
evt2 = evt2(inxT2,:);
evt =[evt1; evt2];

X = [D1.Xsel(inxT1,:); D2.Xsel(inxT2,:)];


X = bsxfun(@rdivide, X, sqrt(sum(X.^2,1)));


icont = 1;
icomp = 3;
Ncv=6;
NC=length(D1.cellinx_sel);
ORI_condset=combnk(ORI_list,2)';

[I11 J11]=sort(P1.C(:,3,icomp,1),'descend');
[I12 J12]=sort(P1.C(:,3,icomp,2),'descend');

[I21 J21]=sort(P2.C(:,3,icomp,1),'descend');
[I22 J22]=sort(P2.C(:,3,icomp,2),'descend');


I = [I11 I21 I12 I22];
J = [J11 J21 J12 J22];

Ws = zeros(NC+1,Ncv+1,2);
for icont = 1: 2
sevts{1} = contrasts(icont);
sevts{2} = ORI_condset(:,icomp);
inxsample = TP.select_subdata(evt,sevts);

inxs_cv = get_Kfoldcvinxs(evt(inxsample,2),Ncv,false);

lam1=0.1;
lam2=0.0;

for icv = 1 : Ncv+1
    if icv == Ncv+1
        inxsm = inxsample(:);
    else
        inxsm = inxsample(inxs_cv.inxs_test{icv});
    end
    trdat = X(inxsm,:);
    trlab = evt(inxsm,2);
    conf.W=[];conf.btrain = true;
    if icv==Ncv+1
        conf.lambda1=lam1*Ncv;
        conf.lambda2=lam2*Ncv;
    else
        conf.lambda1=lam1;
        conf.lambda2=lam2;
    end

    [W, out_tr] = sl_smlr(trdat', trlab,conf);
    Ws(:,icv,icont)=W(:,1);
end

figure; hold on;
plot(Ws(:,1:Ncv,icont));
% figure; hold on;
% plot(mean(Ws(:,1:Ncv/2,icont),2),'b--','LineWidth',2);
% plot(mean(Ws(:,Ncv/2+1:Ncv,icont),2),'r--','LineWidth',2);
% plot(Ws(:,Ncv+1,icont),'LineWidth',2,'Color','k');


end




figure; 
subplot(211);
hold on;
icont=1;
plot(mean(Ws(:,1:Ncv/2,icont),2),'r','LineWidth',2);
plot(mean(Ws(:,Ncv/2+1:Ncv,icont),2),'r--','LineWidth',2);

icont=2;
plot(mean(Ws(:,1:Ncv/2,icont),2),'b','LineWidth',2);
plot(mean(Ws(:,Ncv/2+1:Ncv,icont),2),'b--','LineWidth',2);
plot([1 NC],[0 0])

subplot(212); hold on;
plot(Ws(:,Ncv+1,1),'r.-','LineWidth',2);
plot(Ws(:,Ncv+1,2),'b.-','LineWidth',2);
plot([1 NC],[0 0])

plot([1 NC],[-1 -1])
plot([1 NC],[1 1])




%--------


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
    
    X{1,ievt0}=D1.Xsel(inxsample{1,ievt0},:);
    X{2,ievt0}=D2.Xsel(inxsample{2,ievt0},:);    
    mX(ievt,1,:) = mean(X{1,ievt0},1);
    mX(ievt,2,:) = mean(X{2,ievt0},1);
    stdX(ievt,1,:) = std(X{1,ievt0},0,1);
    stdX(ievt,2,:) = std(X{2,ievt0},0,1);
end
m100= mX(cons==100,:,:);
m40= mX(cons==40,:,:);
s100 = stdX(cons==100,:,:);
s40 = stdX(cons==40,:,:);


cl=[5 10 19 20 33 38 ]%1
cl=[8 9 17 18 25 35 37]

cl=[1 8 10 16 21 27 34 44 47 58]
cl=[1  21 27  29 42 46 57]

n1=80; n2=60
figure;
for k= 1: length(cl)
    subplot(3,3,k);
 hold on;
 icell = cl(k);
errorbar(unique(oris)',m100(:,1,icell),s100(:,1,icell)/sqrt(n1),'r');
errorbar(unique(oris)',m100(:,2,icell),s100(:,2,icell)/sqrt(n2),'r--');
errorbar(unique(oris)',m40(:,1,icell),s40(:,1,icell)/sqrt(n1),'b');
errorbar(unique(oris)',m40(:,2,icell),s40(:,2,icell)/sqrt(n2),'b--');
title(num2str(icell))

end
