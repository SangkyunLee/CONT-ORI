% function get_ORI_conscale_PADEP(calfunstr)
calfunstr='fit_ab4'

dpath{1} = '../GRP_data/AN/thr5/';
dpath{2} = '../GRP_data/AN/thr5/';
dpath{3} = '../GRP_data/AWAKE_EYE/thr5_eyethr_xy1_p1/';
D(1) = load(fullfile(dpath{1},'PADEPDATA_AN1-16-Xsel_ctm0.60.mat'));
D(2) = load(fullfile(dpath{2},'PADEPDATA_AN17-22-Xsel_ctm0.60.mat'));
D(3) = load(fullfile(dpath{3},'PADEPDATA_-Xsel_ctm0.60.mat'));
explist =[1 2 5];
datainxstr = {'AN1-16','AN17-22','AW23-40'};
ctm =0.6

% parpoolid = set_env(true);

cx = [[1 2]' [3 4]' [1 3]' [2 4]' [1 4]' [2 3]'];
for iexp = 1:3
    %fnsave = sprintf('%s_PADEP_ORIsc_ctm%0.2f_%s.mat',datainxstr{iexp},ctm,calfunstr); 
    fnsave = sprintf('L_%s_PADEP_ORIsc_ctm%0.2f_%s.mat',datainxstr{iexp},ctm,calfunstr); 
    iexp0 =explist(iexp);
    [~, ORI_list, ~, nsub, seslist] =get_expinfo(iexp0);
    
    comcont = D(iexp).comcont;
    NO = length(ORI_list);
    ORIsc = cell(size(cx,2),nsub);
    for isub = seslist
        
        for ix = 1 : size(cx,2)
            
            X = cell(1,NO*2);
            cons = zeros(1,NO*2);
            oris = zeros(1,NO*2);
            
            for j = 1 : 2
                X1 = D(iexp).DPA(cx(j,ix),:,isub);
                
                X(1,(j-1)*NO+1 : j*NO) = X1;

                ct = strsplit(comcont{cx(j,ix)},',');
                if strcmp(ct{2},'L')
                    cons1 = -1*str2double(ct{1})*ones(size(X1));
                else
                    cons1 = str2double(ct{1})*ones(size(X1));
                end
                oris((j-1)*NO+1:j*NO) = ORI_list;
                cons((j-1)*NO+1:j*NO) = cons1;
                
                
            end
            conref = cons(1);
            concom = cons(1 +NO);
            len =@(x) size(x,2);
            ns1 = cellfun(len,X);
            inxsample = cell(1, length(ns1));
            k=0;
            for is = 1 : length(ns1)
                inxsample{is}=k+(1:ns1(is));
                k = sum(ns1(1:is));
            end
            X0= cell2mat(X)';
        
        
            evt_cond.inxsample = inxsample;
            evt_cond.cons = cons;
            evt_cond.oris = oris;
            evt_cond.conref = conref;
            evt_cond.concom = concom;
            ORIsc{ix,isub} = cal_wori(X0,evt_cond,calfunstr);
        end
    end
    Out.ORIsc = ORIsc;
    Out.cx=cx;
    Out.comcont = comcont;
    Out.ORI_list = ORI_list;
    Out.PATHR = D(iexp).PATHR;
    Out.sc = mfilename('fullpath');
    
    
    
    fullfnsav = fullfile(dpath{iexp},fnsave);
    save(fullfnsav,'Out','-v7.3') ;

end

if ~isempty(parpoolid),
    delete(parpoolid)
    pause(3);
end
