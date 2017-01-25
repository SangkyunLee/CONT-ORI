dtype='Xsel'
iexp_type=1;
exp_type='AN'
clear dF dV dM
for ises = [ (1:8) 15 16 (17:22)]

[X,mX,stdX,oris,cons] = sortdata(dtype,ises,exp_type);

c40 = find(cons==40);
c100 = find(cons==100);
%fanofactor    
FF40 = stdX(c40,:,:).^2./mX(c40,:,:);
FF100 = stdX(c100,:,:).^2./mX(c100,:,:);    

%var
var40 = stdX(c40,:,:).^2;
var100 = stdX(c100,:,:).^2;    
%mean
M40 = mX(c40,:,:);
M100 = mX(c100,:,:);    


dF{ises,1} = (FF100(:)-FF40(:));
dV{ises,1} = (var100(:)-var40(:));
dM{ises,1} = (M100(:)-M40(:));

end
dF1 = cell2mat(dF);
dV1 = cell2mat(dV);
dM1 = cell2mat(dM);

[mean(dM1) mean(dV1) mean(dF1)]

%--------------------------------
% for ises=15:16
ises=1
iexp_type=1
[contrasts, ORI_list, ORI_compindexset, ~, seslist] =get_expinfo(iexp_type);
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


THR=0.6
[X,mX,stdX,oris,cons] = sortdata(dtype,ises,exp_type);

NC=size(X{1},2);
nC = length(unique(cons));

ORI_condset=combnk(ORI_list,2)';

C1 = squeeze(P1.C(:,end,:,:));
figure; hist(C1(:),100);
% title(num2str(median(C1(:))));
% end
C2 = squeeze(P2.C(:,end,:,:));
Clist =[];
for ic = 1 : nC
    for io = 1 : size(ORI_condset,2)
        Clist = union(Clist,find(C1(:,io,ic)>THR)');
        Clist = union(Clist,find(C2(:,io,ic)>THR)');
    end
end

%%
dtype='Xsel'
iexp_type=1;
exp_type='AN'

ises=1
[X,mX,stdX,oris,cons] = sortdata(dtype,ises,exp_type);

cl=[5 8 12 17 29 35 36 43 45]%1
figure;
for ie = 1:12
    subplot(4,3,ie)
    imagesc(corr(X{1,ie}(:,cl)));
    caxis([0 0.5])
end


