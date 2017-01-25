% function CELLCONT=plot_hist_grpdec
%% plot histogram of group decoding acc data

exp_type={'AN','AN_0TO150','AWAKE'};
iexp_type=1;
ises=15;


if iexp_type==2,
    inx_selcomp=[4 1]; % 4: 0vs 90, 1: 0 vs 30
    contrasts = [100 40];
    sel_compset = [1 2 3 4 7 12 16 19 21];  
    ncellgrp_plot=[3 4 5];
    %ORI_list=[0 30 35 60 90 120 150];
else
    inx_selcomp=[2 3]; % 2: 0vs 90, 3: 0 vs 30
    contrasts = [100 40 20];
    sel_compset = 1:6;
    ncellgrp_plot=[3 4 5];
end


ctm=0.6;
DATA_thr_str = 'thr5';
rThr =0.95;
mThr =0.6;


fnsave = sprintf('SUMMARY_DEC_CELLGRP_ctm%0.2f.mat',ctm);
fullfnsave = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);
summary = load(fullfnsave);
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

            ACC = single(deaccN{incell,icomp,icont});

            a = hist(ACC,40:100); 
            nvals(:,incell0) = a./sum(a); 

            len_HACC = length(find(ACC>=deaccALL(icomp,icont)*100));
            if len_HACC>50
                %cellgrp1 = find(ACC>=deaccALL(icomp,icont)*100 & ACC>(mThr*100));
                if deaccALL(icomp,icont)>mThr
                    cut = deaccALL(icomp,icont)*100;
                else
                    cut =mThr*100;
                end
            else
                %cellgrp1 = find(ACC >= rThr*deaccALL(icomp,icont)*100  & ACC>(mThr*100));
                if rThr * deaccALL(icomp,icont)*100>(mThr*100)
                    cut = deaccALL(icomp,icont)*100*rThr;
                else
                    cut =mThr*100;
                end

            end           
        end

        iplot = icont + 3*(icomp0-1);
        subplot(2,3,iplot); 
        plot(40:100,nvals,'LineWidth',2);
        hold on; plot(deaccALL(icomp,icont)*100*ones(1,2), [0 0.1],'k','LineWidth',2);
        plot(cut*ones(1,2), [0 0.1],'k--','LineWidth',2);
        xlim([50 100]);
        set(gca,'FontSize',22)         
    end    
end


%---------- get CELLCONT
ngcell=3;
ncont=length(contrasts);
CELLCONT = zeros(20,ngcell*ncont*2);
 
icomp=sel_compset(inx_selcomp(2))
kk=1;
maxlen=0;
for icont=1:length(contrasts)    
    for inx_NCG= 1:length(ncellgrp_plot)
        if ~isempty(summary.CELL_CBT{icomp,inx_NCG,icont,ises})
            inxcell = summary.CELL_CBT{icomp,inx_NCG,icont,ises}(:,2) > 0.02;
            M = summary.CELL_CBT{icomp,inx_NCG,icont,ises};
            tmp =[summary.CELLINX{ises}(M(inxcell,1))' M(inxcell,2)*100 ];
            CELLCONT(1:size(tmp,1),kk:kk+1)=tmp;
            maxlen = max(maxlen, size(tmp,1));
        end
        kk = kk + 2;        
    end
    

end
CELLCONT = CELLCONT(1:maxlen,:);

