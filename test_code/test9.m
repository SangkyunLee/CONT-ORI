clear all

D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\WSWC_AN1-16_ORIsc_ctm0.60_fit_ab4.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\WSWC_AN17-22_ORIsc_ctm0.60_fit_ab4.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\WSWC_AW23-40_ORIsc_ctm0.60_fit_ab4.mat');

M(1)=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN1-16_ORIsc_ctm0.60_fit_ab4.mat');
M(2)=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN17-22_ORIsc_ctm0.60_fit_ab4.mat');
M(3)=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\AW23-40_ORIsc_ctm0.60_fit_ab4.mat');





thr = 0.5
A = cell(30,2);
B = cell(30,1);
varA = NaN*ones(30,2);
varB = NaN*ones(30,1);
for icont = 1 : 2
    j = 1;
    for iexp = 1 : 3
        nsub = size(D(iexp).ORIsc,1);
        for isub = 1 : nsub

            X0 = D(iexp).ORIsc{isub,icont};
            if icont==1,
                X1 = M(iexp).ORIsc{isub};
            end
            if isempty(X0)
                continue;
            end

            inx = find(X0.ev>thr);
            p = squeeze(X0.as(1,:,inx));
            A{j,icont} = p;
            varA(j,icont) = var(p);
            
            if icont==1,
                inx = find(X1.ev>thr);
                p = squeeze(X1.as(1,:,inx));
                
                B{j} = p;
                varB(j) = var(p);
            end
            j = j+1;
        end
    end
end
        

varA(varA(:)==0)=NaN;
varB(varB(:)==0)=NaN;





