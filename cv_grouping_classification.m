function [out, opts] = cv_grouping_classification(dataset,label, opts,grplist)
% function out = cv_grouping_classification(dataset,label, opts, grplist)
% 
% INPUT:
%     dataset: 3D matrix [nCell(no channel) X frames X trials]
%     label: 1D vector [1 X trials]
%     grplist: gropulist [Ngrp X Ncomponent]
%     opts:
%         - sel_cl: selected classes
%         - out_ch: excluded channel
%         - Nch: number of channel (cell)
%         - lambda1s: regularization for L2 norm
%         - lambda2s: regularization for L1 norm
%         - mode: classification mode
%                 mode == 1, avg timesamples(one trial --> one prediction)
%         - inxsample: index of samples to be used for classification.
%         - cvmode: cross-validation method ('random','kfold')        
%         - Ncv: number of cross-validation
%         - pertest: percentage of testing data in cross-validation
%         - perval: percentage of validation data in cross-validation
%         - acc_precision: to determine accuracy precision with the decimal places
%           to interprete the accuracy value acc/(10^acc_precision)
%         - classifier: LDA, LDA_opt1 (for fast calculation, this mode
%         should be used)

        
% 2015-07-23, written by Sangkyun Lee
% 2015-10-19, modified by Sangkyun Lee to save memory
%             add acc_precision to determine the decimal places of the accuracy

def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...    
    'mode',1,...
    'inxsample',1,...
    'cvmode','random','Ncv',10,'pertest',0.1,'perval',0.1,'acc_precision',1,'bget_trainacc',false);

if nargin < 3,
	opts = def_opts;
else
	fnms = fieldnames(def_opts);
	for i=1:length(fnms),
		if ~isfield(opts,fnms{i}),
			opts.(fnms{i}) = def_opts.(fnms{i});
		end;
	end;
end;




sel_cl=opts.sel_cl;
out_ch=opts.out_ch;
Nch=opts.Nch;
in_ch=setdiff((1:Nch),out_ch);


Ncv = opts.Ncv;
pertest = opts.pertest;
perval = opts.perval;
inxsample =opts.inxsample;
acc_precision = opts.acc_precision;

label =label(:)';
sel_labinx = false(size(label));
for inx=1:length(sel_cl)
    sel_labinx(label==sel_cl(inx))=true;        
end
label=label(sel_labinx);
dataset=dataset(in_ch,:,sel_labinx);


%%
%

if opts.mode  == 1
    dataset=reshape(mean(dataset(:,inxsample,:),2),[size(dataset,1) size(dataset,3)]);
elseif opts.mode==4    
    dataset=reshape(dataset(:,inxsample,:),[size(dataset,1)*length(inxsample) size(dataset,3)]);
end

%%

X = dataset;
X = X(~isnan(sum(X,2)),:);

if isempty(X)
    out.acc_train = [];
    out.acc_test = [];
    return;
end
%% classifier setting
if strcmp(opts.classifier,'LDA')
    perval=0; %set validation set as [];
end



%% CV mode setting
if strcmp(opts.cvmode,'random')    
    inxs_cv = get_cvinxs_rand(label,Ncv,pertest,perval);    
elseif strcmp(opts.cvmode,'kfold')    
    if perval>0
        bval = true;
        opts.perval = 1/Ncv;
        opts.pertest = 1/Ncv;
    else
        bval = false;
        opts.perval = 0;
        opts.pertest = 1/Ncv;
    end
    inxs_cv = get_Kfoldcvinxs(label,Ncv,bval);
else
    error('Only random or kfold cross-validation mode allowed');
end

bget_trainacc = opts.bget_trainacc;

    
if acc_precision == 0,
    acc_test = zeros(size(grplist,1), Ncv, 'uint8');
    if bget_trainacc
        acc_train = zeros(size(grplist,1), Ncv, 'uint8');
    end
else
    acc_test = zeros(size(grplist,1), Ncv, 'uint16');
    if bget_trainacc
        acc_train = zeros(size(grplist,1), Ncv, 'uint16');
    end
end

for cvinx = 1:Ncv
    
    inx_te=inxs_cv.inxs_test{cvinx};      
    inx_tr=inxs_cv.inxs_train{cvinx}; 

    
    
    switch opts.mode
        case {1, 4}
            trdat = single(X(~isnan(sum(X,2)),inx_tr));
            trlab = single(label(inx_tr));
            tedat = single(X(~isnan(sum(X,2)),inx_te));
            telab = single(label(inx_te));

            if strcmp(opts.classifier,'LDA_opt1')
                
                y0 = trlab(:);                
                y1 = telab(:);
                cond = unique(y0);
                y0b = zeros(size(y0));                
                y0b(y0==cond(1))=1;
                y0b(y0==cond(2))=-1;
                y1b  = zeros(size(y1));                
                y1b(y1==cond(1))=1;
                y1b(y1==cond(2))=-1;
                
                Ncell=size(trdat,1);  
                sc = mean(trdat(:));
                trdatA = [trdat; sc*ones(1, size(trdat,2))];
                tedatA = [tedat; sc*ones(1, size(tedat,2))];

                Ntr = size(trdatA,2);
                Nte = size(tedatA,2);
                Ngrp=size(grplist,1);  

    

                %This code is more than 10 times faster than the following code with the builtin function.    
                %tic
                if acc_precision == 0,
                    if bget_trainacc                        
                        acc_train1=zeros(1,Ngrp,'uint8');
                    end
                    acc_test1=zeros(1,Ngrp,'uint8');
                else
                    if bget_trainacc
                        acc_train1=zeros(1,Ngrp,'uint16');
                    end
                    acc_test1=zeros(1,Ngrp,'uint16');
                end
                  
                M1= trdatA*trdatA';
                Xy1 = trdatA*y0b; 
                    
                parfor igrp = 1 : Ngrp
                    %---- temporal variable for parfor
                    M1_tmp = M1;
                    Xy1_tmp = Xy1;     
                    trdatA_tmp =trdatA;
                    tedatA_tmp =tedatA;
                    %--------------
                    
                    inxcells = [grplist(igrp,:) Ncell+1];
                    M = M1_tmp(inxcells,inxcells);
                    w = M\Xy1_tmp(inxcells);       
                    if bget_trainacc
                        ytr_pre=sign(trdatA_tmp(inxcells,:)'*w);                        
                        res = abs(ytr_pre-y0b);
                        err_tr1 = sum(res,1)/2/Ntr;
                        if acc_precision == 0,
                            acc_train1(igrp)= uint8((1-err_tr1)*100);
                        else
                            acc_train1(igrp)= uint16((1-err_tr1)*100*(10^acc_precision));
                        end
                    end


                    yte_pre=sign(tedatA_tmp(inxcells,:)'*w);                        
                    res = abs(bsxfun(@minus,yte_pre,y1b));
                    err_te1 = sum(res,1)/2/Nte;
                    if acc_precision == 0,
                        acc_test1(igrp)=uint8((1-err_te1)*100);
                    else
                        acc_test1(igrp)=uint16((1-err_te1)*100*(10^acc_precision));
                    end
                end
                
            else
                error('not specified classifier');
            end
        
        otherwise
            fprintf('not implemented yet');
    end
    if bget_trainacc
        acc_train(:,cvinx) = acc_train1;
    end
    acc_test(:,cvinx) =  acc_test1;         
end
if bget_trainacc
    out.acc_train = acc_train;
end
out.acc_test = acc_test;

 
       
        
        



