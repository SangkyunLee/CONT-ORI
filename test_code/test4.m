clear all
close all
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE'};
datainxstr = {'AN1-16','AN17-22','','','AW23-40'};


iexp_type=1
ises=4



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


if iexp_type==2
    P1 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P1-AN17-22_DECSCR_CELLGRP_ctm0.60.mat')
    P1 = P1.S(ises);
    P2 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P2-AN17-22_DECSCR_CELLGRP_ctm0.60.mat')
    P2 = P2.S(ises);
elseif iexp_type==1
    P1 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P1-AN1-16_DECSCR_CELLGRP_ctm0.60.mat')
    P1 = P1.S(ises);
    P2 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_P2-AN1-16_DECSCR_CELLGRP_ctm0.60.mat')
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


% X = bsxfun(@rdivide, X, sqrt(sum(X.^2,1)));


icont = 1;
icomp = 3;
Ncv=10;
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
    else
        conf.lambda1=lam1;
    end
    conf.lambda2=0;
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

figure; hold on;
icont=1;
plot(mean(Ws(:,1:Ncv/2,icont),2),'b','LineWidth',2);
plot(mean(Ws(:,Ncv/2+1:Ncv,icont),2),'b--','LineWidth',2);
% plot(Ws(:,Ncv+1,icont),'k','LineWidth',2);
icont=2;
plot(mean(Ws(:,1:Ncv/2,icont),2),'r','LineWidth',2);
plot(mean(Ws(:,Ncv/2+1:Ncv,icont),2),'r--','LineWidth',2);
% plot(Ws(:,Ncv+1,icont),'k--','LineWidth',2);

plot([1 NC],[0 0])
plot([1 NC],[-1 -1])
plot([1 NC],[1 1])
%------



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


cl=[5 8 12 17 29 35 36 43 45]%1
cl=[4 6 18 28 38 40 42 53 56 60 64 68]%17
cl=[11 13 20 21 26 49 53 ]%18
cl=[7 11 12 18 23 33 36 37 39 41 43 44 46]%19
cl=[3 7 16 20 21 37 43 45 49 59 68 70 76 79 80]%20
n=1
figure;
for k= 1: length(cl)
    subplot(3,3,k);
 hold on;
 icell = cl(k);
errorbar(unique(oris)',m100(:,1,icell),s100(:,1,icell)/sqrt(n),'r');
errorbar(unique(oris)',m100(:,2,icell),s100(:,2,icell)/sqrt(n),'r--');
errorbar(unique(oris)',m40(:,1,icell),s40(:,1,icell)/sqrt(n),'b');
errorbar(unique(oris)',m40(:,2,icell),s40(:,2,icell)/sqrt(n),'b--');
title(num2str(icell))

end

c100= find(cons==100);

a=squeeze(mean(m100,2));
b=squeeze(mean(m40,2));
c=[a; b];
[~,pj]= max(c,[],1);
pj(pj>4)=pj(pj>4)-4;


FF =[];
Var=[];
nevt = length(cons)
ncont = length(unique(cons));
for i=1:2
    for j=1:nevt
        FF(i,j,:)=var(100*X{i,j},0,1)./mean(100*X{i,j},1);
        Var(i,j,:)=var(100*X{i,j},0,1);
    end
end

O=unique(oris);
C=unique(cons);
nc = length(cl);
aa1=[];
aa2=[];
for ic0 = 1 : nc
    ic = cl(ic0);
    inx1 = find(oris==O(pj(ic)));
    inx2 = setdiff(1:length(oris),inx1);
    aa1(:,:,ic0)= FF(:,inx1,ic);
    aa2(:,:,:,ic0) = reshape(FF(:,inx2,ic),[2, length(C), length(O)-1]) ;
end

aa3 = squeeze(mean(aa2,3));

clr={'r', 'b', 'g'}
m={'*','+'}
figure; hold on;
for ic=1:length(C)
    for id=1:2
        x=aa1(id,ic,:);
        y=aa3(id,ic,:);
        plot(x(:),y(:),m{id},'Color',clr{ic});
    end
end
plot([0.01 20],[0.01 20],'k')
xlabel('FF-pref.ORI')
ylabel('FF-other ORIs')

c100= squeeze(aa1(:,1,:));
c40= squeeze(aa1(:,2,:));
figure; hold on;
plot(c100(:),c40(:),'.')
plot([0.01 20],[0.01 20],'k')

%%
FF=zeros(NC,length(oris),2);
for i=1:length(oris);
    for j=1:2
        FF(:,i,j) = var(X{j,i},0,1)./mean(X{j,i},1);
    end
end

FF20 = mean(FF(:,1:3:end,:),2);
FF40 = mean(FF(:,2:3:end,:),2);
FF100 = mean(FF(:,3:3:end,:),2);

FF20 = FF(:,1:3:end,:);
FF40 = FF(:,2:3:end,:);
FF100 =FF(:,3:3:end,:);


%%


P1=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\P1-P2_AN1-16CELLSEL_NCELL_CRSCON_SMLR_L2_ctm0.60.mat')
P2=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\P2-P1_AN1-16CELLSEL_NCELL_CRSCON_SMLR_L2_ctm0.60.mat')
P3=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\P1-P2_AN17-22CELLSEL_NCELL_CRSCON_SMLR_L2_ctm0.60.mat')
P4=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\P2-P1_AN17-22CELLSEL_NCELL_CRSCON_SMLR_L2_ctm0.60.mat')

clear DEC
ses=[(1:8)   (17:22)]%[(1:8) 11 12 15 16]% (17:22) ]    
icomp=3
for icont1 = 1 : 2
    for icont2 = 1 : 2
        for i0= 1:length(ses)
            i = ses(i0);
            if i>16
                DEC(i0,icont2,icont1)=0.5*(P3.DEC_SELCELL{icomp,i}(1,icont2,icont1)+P4.DEC_SELCELL{icomp,i}(1,icont2,icont1));
            else
                DEC(i0,icont2,icont1)=0.5*(P1.DEC_SELCELL{icomp,i}(1,icont2,icont1)+P2.DEC_SELCELL{icomp,i}(1,icont2,icont1));
            end
        end
    end
end

% [a b]=ttest(DEC(:,1,1)-DEC(:,1,2))
% 
% [a b]=ttest(DEC(:,2,2)-DEC(:,2,1))

% anova1([DEC(:,1:2,1) DEC(:,1:2,2)])

[p,t,st] = anova1([DEC(:,1:2,1) DEC(:,1:2,2)]);
[c,m,h,nms] = multcompare(st);

x=DEC(:,1,1);
y=DEC(:,1,2);
label={'100->100','40->100'};

x=DEC(:,2,2);
y=DEC(:,2,1);
label={'40->40','100->40'};
p = ranksum(x,y)

mx=mean([x y]);
sx=std([x y],0,1)/sqrt(length(ses));
figure; hold on;
bar(mx,0.4); errorbar(mx,sx,'.','LineWidth',2);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',label);
ylabel('Dec Acc','FontSize',20)
set(gca,'FontSize',20)

%%


