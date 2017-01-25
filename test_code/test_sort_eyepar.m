
mp_os ={'/media/sdb_WD4T/data_2photon','Y:/data_2photon'}
mainpath = 'Y:\data_2photon\AWAKE'
info.session_path = fullfile(mainpath,'022216BL','05132016_022216BL_2');
info.animal_state = 'awake';
ctm = 0.6
fndata=sprintf('DISK_Cont-ORI_ctm%0.2f_HPF100_tau0.85_F0-10.0sigma.mat',ctm) ;   
datafullfn=fullfile(info.session_path,'/matlab/data',fndata);


% generate A OBJECT    
CO = cont_ori_session(datafullfn,info);
    
    

    
% adding path and filename for eyetrack parameters estimated
CO = get_eyeparfn(CO,scanlist,mp_os)

CO = TP.get_eyeparfn(CO);
eyepar.lmovthr=2;
eyepar.dxythr=1;
eyepar.puppar.ls = 1;
eyepar.puppar.us = 1;
CO = CO.set_eyepar(eyepar);



F = CO.scans(1).Params.files;
mp = F.mainpath;
mp = get_OSpath(mp,mp_os);




XY = cell(CO.no_scan,1);
R = cell(CO.no_scan,1);
dfx = cell(CO.no_scan,1);
vfinx = cell(CO.no_scan,1);
% pooling all samples
XY1 = cell(CO.no_scan,1); 
R1 = cell(CO.no_scan,1); 
PP1= cell(CO.no_scan,1); 

for iscan = 1 : CO.no_scan
    S = CO.scans(iscan);
    S.Params.files.mainpath = mp;
    [XYi, Ri, dPi, PP, vfinxi] = get_eyepar_trial(S.Params,S.timeinfo);
    XY{iscan} = XYi;    
    R{iscan} = Ri;    
    dfX{iscan} = dPi;
    vfinx{iscan} = vfinxi;
    PP1{iscan}=PP;
end


% detect large bias
close all

[XY2, R2] = detect_lmov(XY,R,dfX, 2);
[XY3, R3] = detect_xybias(XY2,R2,1);
nT=zeros(6,2);
for i=1:6

inx0 = ~(cellfun(@isempty,XY{i}));
inx2=~(cellfun(@isempty,XY2{i}));
inx3=~(cellfun(@isempty,XY3{i}));

nT(i,:)=[length(find(inx0)) length(find(inx3))];

vfinx1 = vfinx{i};

inx1 = cell2mat(vfinx1(inx0));
inx2 = cell2mat(vfinx1(inx2));
inx3 = cell2mat(vfinx1(inx3));
figure; plot(PP1{i}(:,1:2));
hold on;
plot(inx1,0.95*mean(PP1{i}(inx1,1))*sign(PP1{i}(inx1,1)),'b.');
plot(inx2,0.95*mean(PP1{i}(inx2,1))*sign(PP1{i}(inx2,1)),'r.');
plot(inx3,0.95*mean(PP1{i}(inx3,1))*sign(PP1{i}(inx3,1)),'k.');
% plot(inx2,(PP1{i}(inx2,1)),'k.');
end
