clear all; close all;
% comcont={'100,H','100,L','40,H','40,L'};
% ctm=0.6
% DATA_thr_str ={'thr0','thr0','','','thr0_eyethr_xy1_p1'}
% for iexp_type= 5%[1 2]
%     get_datawithPA(iexp_type, DATA_thr_str{iexp_type}, comcont, ctm)
% end
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN1-16-Xsel_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN17-22-Xsel_ctm0.60.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\PADEPDATA_-Xsel_ctm0.60.mat');
iexp =1; isub=2;
A = D(iexp).DPA(:,:,isub);

nC = size(A{1,1},1);
k = cell2mat(cellfun(@size,A, 'UniformOutput', false));
T = min(min((k(:,2:2:end))));


mR = max(cell2mat(A(:)));
figure;
for icont = 1: 4
    for icomp = 1: 4
        subplot(4,4,(icont-1)*4+icomp);
        imagesc(A{icont,icomp});
        caxis([0 max(mR)]);
    end
end

figure;
for icont = 1: 4
    for icomp = 1: 4
        subplot(4,4,(icont-1)*4+icomp);
        imagesc(A{icont,icomp});        
    end
end


mm = zeros(nC,16);
for icont = 1: 4
    for icomp = 1: 4
        mm(:,(icont-1)*4+icomp) = mean(A{icont,icomp},2);
    end
end
figure;
imagesc(mm);


mm = zeros(nC*4,4*T);
for icont = 1: 4
    for icomp = 1: 4
        mm((icont-1)*nC+(1:nC),(icomp-1)*T+(1:T)) =A{icont,icomp}(:,1:T);
    end
end
figure;hold on;
imagesc(mm);
for icont = 1: 3    
    plot([1 T*4],[nC*icont nC*icont],'r');
end
for icomp = 1: 3
    plot([T*icomp T*icomp],[1 nC*4],'r');
end


%%

iexp =1;
sublist = [1:8 11 12 15 16]

iexp =2;
sublist = [17:22]


iexp =3;
sublist = [23 25 26 27 29 30 32 33 36 40]


for isub = sublist

    A = D(iexp).DPA(:,:,isub);
    P = cellfun(@mean,A, 'UniformOutput', false);
    nC = size(A{1,1},1);
%     k = cell2mat(cellfun(@size,A, 'UniformOutput', false));
%     T = min(min((k(:,2:2:end))));


    %figure(isub)
    c = zeros(nC,4,4);
    for icont = 1 : 4
        for icomp = 1 : 4
            Ai = A{icont,icomp};
            Pi = P{icont,icomp};
            
            for ic = 1: nC
                c(ic,icont,icomp)=corr(Ai(ic,:)', Pi');
            end                
            
        end
    end
    
end


cinx =unique([unique(cell2mat(W1(2,:))) unique(cell2mat(W2(2,:)))])

M=zeros(2,length(cinx));
for ii=1:length(cinx)
    M(1,ii) = length(find(cell2mat(W1(2,:))==cinx(ii)));
    M(2,ii) = length(find(cell2mat(W2(2,:))==cinx(ii)));
end

figure; imagesc(corr(mean(c,3)))

%% PCA
% iexp =1;
% sublist = [1:8 11 12 15 16]
% 
% iexp =2;
% sublist = [17:22]
% 
% 
% iexp =3;
% sublist = [23 25 26 27 29 30 32 33 36 40]
% 
% 
% for isub = sublist
% 
%     A = D(iexp).DPA(:,:,isub);
%     % P = cell2mat(cellfun(@mean,A, 'UniformOutput', false));
%     nC = size(A{1,1},1);
%     k = cell2mat(cellfun(@size,A, 'UniformOutput', false));
%     T = min(min((k(:,2:2:end))));
% 
% 
%     figure(isub)
%     for icont = 1:4
%         for icomp = 1:4
%             mA = bsxfun(@minus, A{icont,icomp},mean(A{icont,icomp},2));
%             d=flipud(eig(mA*mA'/T));
%             subplot(4,4, (icont-1)*4+icomp);
%             plot(d/max(d),'.-')
%         end
%     end
% end

