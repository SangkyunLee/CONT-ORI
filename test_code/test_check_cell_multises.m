exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE'};

iexp_type=2
ises=17
% P1=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\P1-DEC_CELLGRP_ctm0.60_ses18.mat');
% P2=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\P2-DEC_CELLGRP_ctm0.60_ses18.mat');

datainxstr = {'AN1-16','AN17-22','','','AW23-40'};



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



%-
Mfn1 = sprintf('SUMMARY_%s-%s_DECSCR_CELLGRP_ctm%0.2f.mat',fnpf1,datainxstr{iexp_type},ctm);
Mfn1 = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, Mfn1);
Mfn2 = sprintf('SUMMARY_%s-%s_DECSCR_CELLGRP_ctm%0.2f.mat',fnpf2,datainxstr{iexp_type},ctm);
Mfn2 = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, Mfn2);

M1 = load(Mfn1);
M2 = load(Mfn2);

icomp = 3;
seslist = 17:22;
ORI_list=[-10 0 30 90]
ORI_condset=combnk(ORI_list,2)';


figure; 
for ises0 = 1 : 6
    ises = seslist(ises0);
    Mx1=M1.S(ises);
    Mx2=M2.S(ises);
    C1 = squeeze(Mx1.C(:,end,:,:));
    Ci1 = C1(:,icomp,1);
%     Ci1 = (Ci1-min(Ci1))/(max(Ci1)-min(Ci1));
    C2 = squeeze(Mx2.C(:,end,:,:));
    Ci2 = C2(:,icomp,1);
%     Ci2 = (Ci2-min(Ci2))/(max(Ci2)-min(Ci2));
    subplot(3,2,ises0);        
    plot([Ci1 Ci2])
end

ises = 17
Mx1=M1.S(ises);
Mx2=M2.S(ises);
C1 = squeeze(Mx1.C(:,end,:,:));
Ci1 = C1(:,icomp,1);
% Ci1 = (Ci1-min(Ci1))/(max(Ci1)-min(Ci1));
C2 = squeeze(Mx2.C(:,end,:,:));
Ci2 = C2(:,icomp,1);
% Ci2 = (Ci2-min(Ci2))/(max(Ci2)-min(Ci2));

figure; plot([Ci1 Ci2])
% [I1 J1]=sort(Ci1,'descend');
% [I2 J2]=sort(Ci2,'descend');
% thr=0.5
% inx1 = find(Ci1>=thr);
% inx2 = find(Ci2>=thr);
% 
% iA = union(inx1, inx2);

[I J]=sort((Ci1-Ci2),'descend');
% [I J]=sort(abs(Ci1-Ci2),'descend');


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
m100= mX(2:2:end,:,:);
m40= mX(1:2:end,:,:);

NC=length(D1.cellinx_sel);
for ic = 1 : NC
    if mod(ic,35)==1,
        figure; 
        k=1;
    end
    subplot(5,7,k); hold on;
    plot(m100(:,:,J(ic)));
    plot(m40(:,:,J(ic)),'-.');
    title(num2str(J(ic)));
    k = k+1;
end

orisel = ORI_condset(:,icomp)';

in1 = find(oris==orisel(1) & cons==100);
X11 = D1.Xsel(inxsample{1,in1},J);
X21 = D2.Xsel(inxsample{2,in1},J);

in2 = find(oris==orisel(2) & cons==100);
X12 = D1.Xsel(inxsample{1,in2},J);
X22 = D2.Xsel(inxsample{2,in2},J);

figure;
subplot(2,2,1);
imagesc(corr(X11)); caxis([0 0.5]);
subplot(2,2,2);
imagesc(corr(X21)); caxis([0 0.5]);
subplot(2,2,3);
imagesc(corr(X12)); caxis([0 0.5]);
subplot(2,2,4);
imagesc(corr(X22)); caxis([0 0.5]);

figure;
subplot(411);
plot([mean(X11,1)' mean(X21,1)'])
subplot(412);
plot([mean(X12,1)' mean(X22,1)'])
subplot(413)
plot([var(X11,0,1)'  var(X21,0,1)'])
subplot(414)
plot([var(X12,0,1)'  var(X22,0,1)'])





%%



icomp=3

for ises0 = 1 : 6

    ises = seslist(ises0);
    Mx1=M1.S(ises);
    Mx2=M2.S(ises);
    C1 = squeeze(Mx1.C(:,end,:,:));
    Ci1 = C1(:,icomp,1);
    C2 = squeeze(Mx2.C(:,end,:,:));
    Ci2 = C2(:,icomp,1);

end

[I J]=sort((Ci1-Ci2),'descend');


y1 = quantile(Ci1,[ .25 .50 .75]); 
y2 = quantile(Ci2,[ .25 .50 .75]); 

i1 = find(Ci1<y1(1));
i1h = find(Ci1>y1(3));
i2 = find(Ci2<y2(1));
i2h = find(Ci2>y2(3));

Hd = [intersect(i1,i2h); intersect(i1h,i2)];

iL = union(i1,i2);
iH = union(i1h,i2h);
[iL Ci1(iL)-Ci2(iL)]
Ci1(iH)-Ci2(iH)



