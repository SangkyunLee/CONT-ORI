clear all
% close all

exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
iexp_type=1;
ises=15;


if iexp_type == 3,
    inx_selcomp=[4 1]; % 4: 0vs 90, 1: 0 vs 30
    contrasts = [100 40];
    sel_compset = [1 2 3 4 7 12 16 19 21];  
    ncellgrp_plot=[3 4 5];
    %ORI_list=[0 30 35 60 90 120 150];
else
    inx_selcomp=[2 3 6]; % 2: 0vs 90, 3: 0 vs 30 4: -15 vs 0
    contrasts = [100 40];%[100 40 20];
    sel_compset = 1:6;
    ncellgrp_plot=[3 4 5];
end


ctm=0.6;
DATA_thr_str = 'thr5';
rThr =0.95;
mThr =0.6;
%------load 
fntmp = {'AN1-16','AN17-22','','',''};
fnsave = sprintf('%sSUBCELL-%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type},'Xsel',ctm); %L2norm -regularization
RLR_fn= fullfile(fileparts(pwd),'NEW_DECODING_NOBIAS_ZMEAN',exp_type{iexp_type},DATA_thr_str, fnsave); 
R=load(RLR_fn);

fndata = sprintf('DEC_CELLGRP_ctm%0.2f_ses%d.mat',ctm,ises);
fullfndata = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fndata);
fullfndata1 = strrep(fullfndata,'DEC_CELLGRP','CELLGRP');

load(fullfndata);
load(fullfndata1);


figure('Position',[100 100 1000 500]);

for icomp0 = 1: length(inx_selcomp)
        icomp =sel_compset(inx_selcomp(icomp0)); 
    for icont = 1 : length(contrasts)    
        
        nvals = zeros(61,length(ncellgrp_plot));
        for incell0 = 1:length(ncellgrp_plot)
            incell = ncellgrp_plot(incell0);

            ACC = single(deaccN{incell,icomp,icont})/100;

            a = hist(ACC,[40:100]/100); 
            nvals(:,incell0) = a./sum(a);              
        end

        iplot = icomp0 + 3*(icont-1);
        subplot(2,3,iplot); 
        plot([40:100]/100,nvals,'LineWidth',2);
        hold on; plot(deaccALL(icomp,icont)*ones(1,2), [0 0.1],'k','LineWidth',2);
        
%         Rd = R.DEC_M(icont,icomp,ises)*100
%         plot(Rd*ones(1,2), [0 0.1],'k--','LineWidth',2);
        set(gca,'XTick',[0.5 0.75 1]);

        xlim([0.5 1]);
        set(gca,'FontSize',22)         
    end    
end


%
load('/home/slee/data/data_2photon/matlab_2ndLev/NEW_DECODING/AN/thr5/SUMMARY_AN1-16_DEC_CELLGRP_ctm0.60.mat')
M=load('/home/slee/data/data_2photon/matlab_2ndLev/NEW_DECODING/AN/thr5/AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
figure; plot(S(ises).C(:,:,iori,icont))

[~,inx] = sort(S(ises).C(:,5,iori,icont),'descend');
A=M.CELLSEL_INX{icont,iori,ises};
cinx = cell2mat(A(4,:));
cuinx = unique(cinx);
[a b]=hist(cinx,cuinx);
[~,inx2]=sort(a,'descend');
b(inx2)


