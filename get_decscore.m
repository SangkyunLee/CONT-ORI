
%% analysis of group decoding acc data
% clear all
% close all
% 
% bpart=true;
% iexp_type=1;

function get_decscore(iexp_type,part,thr)

if exist('part','var') && ~isempty(part)
    fnpf=part;    
    bpart = true;
else
    fnpf='';
    bpart = false;
end



exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
datainxstr = {'AN1-16','AN17-22','','','AW23-40'};
ncellgrp_plot=3:5;
ctm=0.6;
DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'};


[contrasts, ORI_list, sel_compset, nses, seslist] =get_expinfo(iexp_type);






if bpart
    fnsave = sprintf('SUMMARY_%s-%s_DECSCR_CELLGRP_ctm%0.2f.mat',fnpf,datainxstr{iexp_type},ctm);
else
    fnsave = sprintf('SUMMARY_%s_DEC_CELLGRP_ctm%0.2f.mat',datainxstr{iexp_type},ctm);
end
fullfnsave = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str{iexp_type}, fnsave);


S.C=[];
S.Ncell =0;
S.deaccALL=[];
S.cellinx_sel=[];
S.ncellgrp = [];
S.thr=[];
S = repmat(S,[1 nses]);

for ises0 = 1: length(seslist)
    ises = seslist(ises0);
    fprintf('ises: %d\n',ises);
    if bpart
        fndata = sprintf('%s-DEC_CELLGRP_ctm%0.2f_ses%d.mat',fnpf,ctm,ises);
    else
        fndata = sprintf('DEC_CELLGRP_ctm%0.2f_ses%d.mat',ctm,ises);
    end
    fullfndata = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str{iexp_type}, fndata);
    fullfndata1 = strrep(fullfndata,'DEC_CELLGRP','CELLGRP');
    
    
    
    D1=load(fullfndata);
    D2=load(fullfndata1);

    nc = D1.Ncell; 
    deaccN = D1.deaccN;
    deaccALL = D1.deaccALL;
    grplistALL = D2.grplistALL;
    L = 1:nc;
    C = zeros(nc,length(thr),length(ncellgrp_plot),length(sel_compset),length(contrasts));

    if exist(fullfnsave,'file')
        load(fullfnsave);
        C = S(ises).C;
    end
    if ~(exist(fullfndata,'file') && exist(fullfndata1,'file')) 
        continue;
    end
    for icont = 1 : length(contrasts)        
        for icomp2 = 1 : length(sel_compset) 
            icomp =sel_compset(icomp2);                      
     
            
            for incell0 = 1  : length(ncellgrp_plot)          
                
                
                if ~any(sum(C(:,:,incell0,icomp,icont),1)==0),
                    continue;
                end
                fprintf('icont:%d, icomp:%d, incell0%d\n', icont,icomp2, incell0);                
                
                incell = ncellgrp_plot(incell0);
                
                
                ACC = single(deaccN{incell,icomp,icont});
                M = round(deaccALL(icomp,icont)*100);
                A = grplistALL{incell};
                K = A(:);
                
                for ithr = 1 : length(thr)
                    s = (ACC-50)/(M-50);
                    s(s<thr(ithr))=NaN;
                    s = repmat(s, [1 incell]);    
                    J = zeros(nc,1);
                                   

                    for ii = L(:)'
                        J(ii)=nansum(s(K==ii));
                    end

                    C(:,ithr,incell0,icomp2,icont) = J;
                end
            end

        end
    end
    S(ises)=struct('C',C,'Ncell',nc,'deaccALL',deaccALL,'cellinx_sel',D1.cellinx_sel,'ncellgrp',ncellgrp_plot,'thr',thr);
    save(fullfnsave,'S');
end



