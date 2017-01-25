% iexp_type = 1;
% metric ='var'
% selcont = [1 2];
% %
% 
% D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\BASIC_SUMMARY_AN1-16_ctm0.60.mat');
% D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\BASIC_SUMMARY_AN17-22_ctm0.60.mat');
% D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\BASIC_SUMMARY_AW23-40_ctm0.60.mat');
% V1=cell(3,1);
% for iexp = 1 : 3
%     n = length(D(iexp).S);
%     v = cell(n,1);    
%     for ises= 1: n        
%         switch metric
%             case 'snr'
%                 V = D(iexp).S(ises).m./sqrt(D(iexp).S(ises).V);
%             case 'var'
%                 V = D(iexp).S(ises).V;
%             case 'mean'
%                 V = D(iexp).S(ises).m;
%             case 'fano'
%                 V = D(iexp).S(ises).V./D(iexp).S(ises).m;
%             case 'nr'
%                 V = D(iexp).S(ises).NR;                
%             case 'sr'
%                 V = D(iexp).S(ises).SR;
%         end
%         
%         if ~isempty(V)
%             V = mean(V,1);
%             K = squeeze(mean(V(:,:,selcont),2)); %mean across orientations
%             v{ises}=K';
%             
% 
%             
% %             K = squeeze(mean(V(:,:,selcont),2)); %mean across orientations
% %             v{ises}=K;
%             
%         end
%     end
%     v = cell2mat(v);    
%     V1{iexp} = v;
% end
% % V1= cell2mat(V1);
% % p = ranksum(V1(:,1),V1(:,2))  
% % anova1(V1)
% 
% 
% %----------------------------------
% K1(1) =load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_AN1-16_DEC_CELLGRP_ctm0.60.mat');
% K1(2) =load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUMMARY_AN17-22_DEC_CELLGRP_ctm0.60.mat');
% K1(3) =load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AWAKE_EYE\thr5_eyethr_xy1_p1\SUMMARY_AW23-40_DEC_CELLGRP_ctm0.60.mat');
% 
% 
% K0(1) =load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% K0(2) =load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% K0(3) =load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% 


method ='smlr'
% method ='descore'
% incell = 3; ithr=2
ncs =10;
icont1 = 1; % 100%
icont2 = 2; % 40% contrast
icont3 = 3; % dummy for randomly selected cellinx
stattype ='V'

M = cell(3,6);

for icomp = 1 : 6
    
    for iexp = 1:3
        ORIcomp = combnk(D(iexp).ORI_list,2)';
        
        icomp1 =find(D(iexp).ORI_list==ORIcomp(1,icomp));
        icomp2 =find(D(iexp).ORI_list==ORIcomp(2,icomp));
    
        if strcmp(method,'smlr')
            K = K0(iexp);
            nses = size(K.WEIGHT,3);
        elseif strcmp(method,'descore')
            K = K1(iexp);
            nses = length(K.S);
        end
        
        M0 = cell(nses,1);
        for ises = 1: nses
            x =D(iexp).S(ises).(stattype);
            if isempty(x)
                continue;
            end
            fprintf('icomp:%d, iexp:%d, ises:%d\n',icomp, iexp, ises);
            %inx3: distinctive cells
            switch method
                case 'smlr'
                    mean2 = @(x)mean(x,2);
                    C = cellfun(mean2,K.WEIGHT([icont1 icont2],icomp,ises),'UniformOutput', false);
                    C = cell2mat(C');
                case 'descore'                   
                    C = squeeze(K.S(ises).C(:,ithr,incell,icomp,[icont1 icont2]));                    
            end
            if isempty(C)||sum(C(:))==0
                fprintf('\nempty C: iexp%d, icomp%d, ises%d\n',iexp,icomp,ises);
            end
            
%             C0 = bsxfun(@rdivide,C,max(abs(C),[],1))
%             [~,inx1] =sort(abs((C0(:,1)-C0(:,2))),'descend');
%             [~,inx2] =sort(abs(C0(:,1)+C0(:,2)),'descend');
%             
%             %%--------------
%             C0 = bsxfun(@rdivide,C,max(abs(C),[],1));
%             [~,inx1] =sort(abs((C0(:,1)-C0(:,2))),'descend');
%        
%             [~,inx21] =sort(abs(C(:,1)),'descend');
%             [~,inx22] =sort(abs(C(:,2)),'descend');
%             
%             inx2 =intersect(inx21(1:ncs),inx22(1:ncs));
%             
%             inx3 = setdiff(inx1(1:ncs),inx2);
%             inx4 = setdiff(inx2,inx1(1:ncs));
%             [~,inx5] = sort(abs(C(:,1)+C(:,2)),'ascend');
%             inx5 = inx5(1:ncs);
%             inx3 = setdiff(inx3,inx5);
            %---------------------------------
            C0 = bsxfun(@rdivide,C,max(abs(C),[],1));
            [~,inx1] =sort(abs(C0(:,1)),'descend');
            [~,inx2] =sort(abs(C0(:,2)),'descend');
            
            inx0 = unique([inx1(1:ncs)' inx2(1:ncs)']);
            inx4 = intersect(inx1(1:ncs)', inx2(1:ncs)');
            inx3 = setdiff(inx0,inx4);
 

            [~,inx5] = sort(abs(C0(:,1)+C0(:,2)),'ascend');
            inx5 = inx5(1:ncs);
            inx3 = setdiff(inx3,inx5);
%             nx =min(length(inx3), length(inx4));
%             inx3 = inx3(1:nx);
%             inx4 = inx4(1:nx);
%             inx5 = inx5(1:nx);

            M1=[]; M2=[]; M3=[];
            switch stattype
                case 'm'
                    M1 = x(inx3,icomp1,[icont1 icont2])- x(inx3,icomp2,[icont1 icont2]);
                    M1= squeeze(abs(M1));
                    if size(M1,2)==1,
                        M1= M1';                        
                    end
                    M2 = x(inx4,icomp1,[icont1 icont2])- x(inx4,icomp2,[icont1 icont2]);
                    M2= squeeze(abs(M2));
                    if size(M2,2)==1,
                        M2= M2';                        
                    end
                    M3 = x(inx5,icomp1,[icont1 icont2])- x(inx5,icomp2,[icont1 icont2]);
                    M3= squeeze(abs(M3));
                    if size(M3,2)==1,
                        M3= M3';                        
                    end

                    M1 =abs(M1*[1 -1]');%./M1*[1 1]';
                    M2 =abs(M2*[1 -1]');%./M2*[1 1]';
                    M3 =abs(M3*[1 -1]');%./M3*[1 1]';
                    
                case 'V' %variance
                    M1 = x(inx3,icomp1,[icont1 icont2])+ x(inx3,icomp2,[icont1 icont2]);
                    M1= squeeze((M1/2));
                    if size(M1,2)==1,
                        M1= M1';                        
                    end
                    M2 = x(inx4,icomp1,[icont1 icont2])+ x(inx4,icomp2,[icont1 icont2]);
                    M2= squeeze((M2/2));
                    if size(M2,2)==1,
                        M2= M2';                        
                    end
                    M3 = x(inx5,icomp1,[icont1 icont2])+ x(inx5,icomp2,[icont1 icont2]);
                    M3= squeeze((M3/2));
                    if size(M3,2)==1,
                        M3= M3';                        
                    end

                   M1 =abs(M1*[1 -1]');%./M1*[1 1]';
                    M2 =abs(M2*[1 -1]');%./M2*[1 1]';
                    M3 =abs(M3*[1 -1]');%./M3*[1 1]';
            end
                    
            %M0{ises}=[M1 M2];
            
            if ~isempty(M1) && ~isempty(M2) && ~isempty(M3) 
                M0{ises}= [mean(M1) mean(M2) mean(M3,1)];
            else
                M0{ises}= NaN*ones(1,3);
                fprintf('\nempty M2: iexp%d, icomp%d, ises%d\n',iexp,icomp,ises);
                
            end

        end

        M{iexp,icomp}=cell2mat(M0);    
        
    end
    
end


Y = cell(6,1);

for icomp = 1:6
    tmp2 = cell2mat(M(:,icomp));
    Y{icomp} = tmp2;
    if icomp==1
        am =  isnan(tmp2(:,1));        
    else        
        am = am | isnan(tmp2(:,1));
    end
end
am = repmat(am,[6 1]);
Y1 = cell2mat(Y);
nrep =(length(Y1)-sum(am))/6;
[p,table,stats] =anova2(Y1(~am,:),nrep);

 [c,m,h,nms] = multcompare(stats);

 

[P,ANOVATAB,stats] = kruskalwallis(Y{3});
 [c,m,h,nms] = multcompare(stats)
 
 [c,m,h,nms] = multcompare(stats)

