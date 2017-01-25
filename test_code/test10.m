clear all; close all;
% comcont={'100,H','100,L','40,H','40,L'};
% ctm=0.6
% DATA_thr_str ={'thr0','thr0','','','thr5_eyethr_xy1_p1'}
% for iexp_type= [1 2]
%     get_datawithPA(iexp_type, DATA_thr_str{iexp_type}, comcont, ctm)
% end
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr0\PADEPDATA_AN1-16-Xsel_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr0\PADEPDATA_AN17-22-Xsel_ctm0.60.mat');
iexp =1; isub=1;
A = D(iexp).DPA(:,:,isub);
nC = size(A{1,1},1);
mR = max(cell2mat(A(:)));
figure;
for icont = 1: 4
    for icomp = 1: 4
        subplot(4,4,(icont-1)*4+icomp);
        imagesc(A{icont,icomp});
        caxis([0 max(mR)]);
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
