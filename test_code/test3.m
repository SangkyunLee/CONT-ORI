% clear all
% eye_dir = 'Y:\data_2photon\AWAKE\0822BRLF\11082015_0822BRLF\eye11082015'
% % eye_dir='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015';


% eye_dir = 'Y:\data_2photon\AWAKE\0822BRLF\11082015_0822BRLF_2\eye11082015_2'
% recordnum=4;
% 
% 
% [daqdata,vidtime,Info] = load_eyedaq(eye_dir,recordnum);
% [stimstart, saminx] = identify_stimstart(daqdata(:,1), daqdata(:,2));

% TPdaqfr = data.Params.samplingfreq_NI;

clear all


pre_stimtime_ms=500;
motionthr.twin = 1;
motionthr.IMG_motion_thr =2;
motionthr.dist_twin_thr = 1;

ctm=0.6
fndata =sprintf('RIM_Cont-ORI_ctm%0.2f_HPF100_tau0.85_F0-5.0sigma.mat',ctm) ;
datafullfn = fullfile(locpath,fndata)
info.animal_state='AW';
% info.session_path = path_list{inx_session};
% generate A OBJECT    
CO = cont_ori_session(datafullfn,info);            
scanlist=1:CO.no_scan;


logfn = dir(fullfile(CO.session_path,'doc','*.xls'))
logfn =fullfile(CO.session_path,'doc',logfn.name)
sheet=1
finfo = load_logfile(logfn, sheet);

for iscan = scanlist
    Params = CO.scans(iscan).Params;
    [~,Imgloc] = fileparts(Params.files.subpath_xml);
    inx=strfind(Imgloc,'TSeries-');
    if isempty(inx)
        error('incorrect path');
    end
    Imgloc = Imgloc(inx:end);
    for iscan2 =1 : length(finfo)
        imdir = finfo(iscan2).Image_directory;
        if ~strcmp(Imgloc,imdir)
            continue;
        end
        eyefn = finfo(iscan2).eye;        
        eye(iscan).eyefn = eyefn;
    end
end

%------ calculate_motion_rotaroadNimageframe   
if isstruct(motionthr) && length(motionthr)==1
    twin = motionthr.twin;
    CO = CO.set_motion_rotaroad_imageframe(twin);
end

%%
iscan=7
motion= CO.motion(iscan);
timeinfo = CO.scans(iscan).timeinfo;
Params = CO.scans(iscan).Params;
frametime = find(timeinfo.frame_start==10)/Params.samplingfreq_NI;
stimstart_IMG = find(timeinfo.stimtime>0,1)/Params.samplingfreq_NI;



recordnum = str2num(eye(iscan).eyefn)
dlist = dir(fullfile(Params.files.mainpath,'eye*'));
eye_dir = fullfile(Params.files.mainpath,dlist.name)


locpath= fileparts(eye_dir);

locpath = fullfile(fileparts(eye_dir),'matlab','eyepar');

sfn = dir(fullfile(locpath, sprintf('%04d_smallfitinfo*.mat',recordnum)));
if ~isempty(sfn)
    sfn =sfn(end).name;
    M = load(fullfile(locpath,sfn));
    [inxbadfit, pars] =get_badfit(M.fitInfo,[25 80 15 5]);
    figure; plot(pars(3,:),'.'); hold on; plot(inxbadfit, pars(3,inxbadfit),'k.')
    P= pars;
%     P(:,inxbadfit)=NaN;
end

[daqdata,vidtime,Info] = load_eyedaq(eye_dir,recordnum);
[stimstart_eyedaq, saminx] = identify_stimstart(daqdata(:,1), daqdata(:,2));


timeoffset = stimstart_eyedaq - stimstart_IMG;
rTime_eye = vidtime-timeoffset;


Nframe = length(frametime)
Pframe = zeros(3,Nframe);
for ifr = 1 : Nframe-1
    tst = frametime(ifr);
    ted = frametime(ifr+1);
    vix = find(rTime_eye>=tst & rTime_eye<ted);
    Pframe(:,ifr) = nanmean(P(:,vix),2);
end
Pframe(:,Nframe)=P

M1 =CO.motion(iscan).IMG_motion;
M2 =CO.motion(iscan).dist_twin_frame;
% 
% M =[M1';M2';Pframe];
% K = nanmean(M,2) + 3*nanstd(M,0,2);
% nM = bsxfun(@rdivide,M,K);
figure; 
subplot(311);
plot(Pframe(1:2,:)')
subplot(312);
plot(Pframe(3,:)')
subplot(313);
plot([M1 M2])
%%
% 
% onset_tstamp_vector = zeros(1, length(daqtime));
% stim_onset_tstampinNI(1,1) = start;
% onset_tstamp_vector(1,start) = 10;
% 
% Nrep = stimparam.repetitions;
% blank_sample = stimparam.blank_samplesinNI/TPdaqfr*daqfr;
% stim_sample = stimparam.stim_samplesinNI/TPdaqfr*daqfr;
% 
% j = start;
% for i = 1 : Nrep
%     if blank_sample>0
%         j = j + round(stim_sample*0.9) ;
%         if j<length(flatPDsig)
%             while flatPDsig(j)<thresh_PD
%                 if j<length(flatPDsig)
%                     j=j+1;
%                 else
%                     break;
%                 end
%             end
%         end
%         if j<length(flatPDsig)           
%             onset_tstamp_vector(1,j) = -10;
%         end
%     end
%     if i<stimparam.repetitions
%         if blank_sample>0
%             j = j + round(blank_sample*0.9) ;
%         else
%             j = j + round(stim_sample*0.9) ;
%         end
%         if j<length(flatPDsig)
%             while flatPDsig(j)<thresh_PD
%                 if j<length(flatPDsig)
%                     j=j+1;
%                 else
%                     break;
%                 end
%             end        
%         end
%         if j<length(flatPDsig)
%             stim_onset_tstampinNI(1,i+1)=j;
%             onset_tstamp_vector(1,j) = 10;
%         end
%     end    
% end
% 
% 
% 
% len=length(onset_tstamp_vector);
% len2 = length(flatPDsig);
% onset_tstamp_vector2 =[onset_tstamp_vector'; zeros(len2-len,1)];
% len2 = min(len2, length(onset_tstamp_vector2));
% hfig=figure; plot(daqtime(1:len2),[onset_tstamp_vector2(1:len2)/5*(max(flatPDsig(1:len2))) flatPDsig(1:len2)])