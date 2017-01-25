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
scanlist=1:7;
CO = get_eyeparfn(CO,scanlist,mp_os);

CO = TP.get_eyeparfn(CO);


%------ apply all processes used in data analysis
eyepar.lmovthr=2;
eyepar.dxythr=1;
eyepar.puppar.ls = 1;
eyepar.puppar.us = 1;
CO = CO.set_eyepar(eyepar);


motionthr.twin = 1;
motionthr.IMG_motion_thr =2;
motionthr.dist_twin_thr = 1;




        
%------------------ load rotaroad
iscan =1;
channel=[3 5];
files = CO.scans(iscan).Params.files;
DAQpath = fullfile(CO.session_path,files.subpath_xml);
FMOTpath = DAQpath;
frame_period = CO.scans(iscan).Params.msperframe/1000;
pixel_resolution = CO.scans(iscan).ROI{1}(1).pixelres;
params=struct('DAQpath',DAQpath,...
    'FMOTpath',FMOTpath,...
    'frame_period',frame_period,...
    'pixel_resolution',pixel_resolution,...
    'twin',twin,...
    'channel',channel);
[~, dist_twin, ~, ~]= TP.TPSession.extract_motion_combo1(params);
tinfo =CO.scans(iscan).timeinfo;
DAQinx_1ststimon = find(tinfo.stimtime>0,true,'first');
% rotaroad time aligned to 1st stimonset
t_rota = ((1:length(dist_twin))-DAQinx_1ststimon)/CO.scans(iscan).Params.samplingfreq_NI;

%--------------------------------------


eye1 = CO.eye(1);
xy = eye1.PP(:,1:2);
P = eye1.PP(:,3);
inxn0 = xy(:,1)>0;
xy0 = bsxfun(@minus, xy,mean(xy(inxn0,:),1));
% eyemovie time aligned to 1st stimonset
t=((1:length(xy0))-eye1.vf_trial{1}(1))/30;



figure('Position', [680 678 500 300]);
hold on;
u = 3.4/410; % I assume the eye horizontal 3.4mm and the pixel : 410
ts=1500:2900;
inx_rota = find(t_rota>=t(ts(1)) & t_rota<=t(ts(end)));




%------- plot visual stim selected----------

nt = length(eye1.vf_trial);
S = zeros(nt,1);
for i = 1 : nt
    if ~isempty(eye1.vf_trial{i}) && (eye1.vf_trial{i}(1)>=ts(1) && eye1.vf_trial{i}(end)<=ts(end))
        S(i)=1;
    end
end
fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.99);
k=eye1.vf_trial(S==1);

% figure('Position', [680 678 500 300]); hold on;
ma1=-0.5; ma2=1.5;
for i = 1 : length(k)
    ts2 = k{i};
    % eyetracking sample time is prestim(100ms) + stim(500ms) + poststim(100ms)
    %therefore, here I present only for stim time
    t2 = t(ts2(4:end-3));
    
    fill_rec(t2([1 end]),[ma1 ma1],[ma2 ma2],[0.9 0.9 0.9] );
end

rota0 =0.7;
tbase = t(ts(1)-100:ts(end)+100);
plot(tbase,rota0*ones(1,length(tbase)),'k');
plot(tbase,zeros(1,length(tbase)),'k');
plot(t(ts),xy0(ts,1)*u,'r','LineWidth',1.5); %eye x,y
plot(t(ts),xy0(ts,2)*u,'b','LineWidth',1.5); %eye x,y
plot(t(ts),P(ts)*u,'k','LineWidth',1.5); %pupil radius
plot(t_rota(inx_rota),dist_twin(inx_rota)/40+rota0,'y','LineWidth',1.5) % rotaroad


%--- plot sample point of eyeframe
si =[1715 1886 1947 2215]
plot(t(si),xy0(si,1)*u,'LineStyle','none','Marker','.','MarkerSize',20,'Color',[0.5 0.5 0.5]);
plot(t(si),P(si,:)*u,'LineStyle','none','Marker','.','MarkerSize',20,'Color',[0.5 0.5 0.5]);
set(gca,'FontSize',16);
xlim([25 75])



%---------- figure- summary statistics
%--- apply only sccade detection
eyepar.lmovthr=2;
eyepar.dxythr=100;
eyepar.puppar.ls = 100;
eyepar.puppar.us = 100;
CO = CO.set_eyepar(eyepar);
eyexy = CO.eye;

i=2
Params = CO.scans(i).Params;
timeinfo = CO.scans(i).timeinfo;
Params.files.mainpath = CO.session_path;
[XYi, Ri, dPi, PP, vfinxi] = TP.get_eyepar_trial(Params,timeinfo);
eye(i).xy = XYi;    
eye(i).r = Ri;    
eye(i).PP = PP;
eye(i).vf_trial = vfinxi;
dfX = dPi;                

eye1 = TP.detect_lmov(eye(i),dfX, 2);

PP = eye(2).PP;
t = (1:length(PP))/30;
inxt = 4000:12000;
XY0 = bsxfun(@minus,PP(:,1:2), median(PP(:,1:2),1));
figure; hold on;
plot(t(inxt),XY0(inxt,1)*u,'r');
plot(t(inxt),XY0(inxt,2)*u,'b');
