
%% analysis of group decoding acc data
clear all
close all

bpart=true;
if bpart
    fnpf='P1'
else
    fnpf=''
end


exp_type={'AN','AN_0TO150','AWAKE'};



%--- AN
datainxstr='AN17-22'
ORI_list=[-10 0 30 90];
contrasts = [ 100 40];
sel_compset = 1:6;

iexp_type=1;
seslist=17:22;



ctm=0.6;
DATA_thr_str = 'thr5';


if bpart
    fnsave = sprintf('SUMMARY_%s-%s_DECSCR_CELLGRP_ctm%0.2f.mat',fnpf,datainxstr,ctm);
else
    fnsave = sprintf('SUMMARY_DEC_CELLGRP_ctm%0.2f.mat',ctm);
end
fullfnsave = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);


S.C=[];
S.Ncell =0;
S.deaccALL=[];
S.cellinx_sel=[];
S.ncellgrp = [];
S = repmat(S,[1 max(seslist)]);

for ises0 = 1: length(seslist)
    ises = seslist(ises0);
    if bpart
        fndata = sprintf('%s-DEC_CELLGRP_ctm%0.2f_ses%d.mat',fnpf,ctm,ises);
    else
        fndata = sprintf('DEC_CELLGRP_ctm%0.2f_ses%d.mat',ctm,ises);
    end
    fullfndata = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fndata);
    fullfndata1 = strrep(fullfndata,'DEC_CELLGRP','CELLGRP');
    
    if ~(exist(fullfndata,'file') && exist(fullfndata1,'file'))
        continue;
    end
    
    D1=load(fullfndata);
    D2=load(fullfndata1);

    nc = D1.Ncell; 
    deaccN = D1.deaccN;
    deaccALL = D1.deaccALL;
    grplistALL = D2.grplistALL;
    L = 1:nc;
    C = zeros(nc,length(ncellgrp_plot),length(sel_compset),length(contrasts));
    for icont = 1 : length(contrasts)        
        for icomp2 = 1 : length(sel_compset) 
            icomp =sel_compset(icomp2);            
            
            
            for incell0 = 1  : length(ncellgrp_plot)
                incell = ncellgrp_plot(incell0);
                
                ACC = single(deaccN{incell,icomp,icont});
                M = round(deaccALL(icomp,icont)*100);
                A = grplistALL{incell};
                
                
                s = (ACC-50)/(M-50);
                s(s<0)=0;
                s = repmat(s, [1 incell]);                
                
                J = zeros(nc,1);
                K = A(:);               
                
                for ii = L(:)'
                    J(ii)=mean(s(K==ii));
                end
                
                C(:,incell0,icomp2,icont) = J;
            end

        end
    end
    S(ises)=struct('C',C,'Ncell',nc,'deaccALL',deaccALL,'cellinx_sel',D1.cellinx_sel,'ncellgrp',sort(D1.ncellgrp));
    
end



