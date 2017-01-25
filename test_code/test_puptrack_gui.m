clear all;close all;


mainpath = 'Y:\data_2photon\AWAKE\thy1_0129B3\05102016_thy1_0129B3_1'


logfn='Y:\data_2photon\AWAKE\thy1_0129B3\05102016_thy1_0129B3_1\doc\05102016_thy1_0129B3_1.xls'
sheet=1
finfo = load_logfile(logfn, sheet);

exptypestr =   'Cont-ORI' % 'ORI' %
list_selT=[];
for ifinfo=1:length(finfo)
    if length(finfo(ifinfo).Experiment)>=length(exptypestr)...
            & strcmp(finfo(ifinfo).Experiment(1:length(exptypestr)),exptypestr)...
            & (strcmp(finfo(ifinfo).use,'1') | finfo(ifinfo).use==1) 
       list_selT = [list_selT ifinfo];
    end
    
end


eye_dir = fullfile(mainpath,'eye0129B3_1')
i=list_selT(1)

motion_dir = fullfile(mainpath, 'data',finfo(i).Image_directory)
motionfn = fullfile(motion_dir,'motionpar*.mat')
mtn =dir(motionfn);

[wheel, time_sec]=load_wheel(motion_dir);



videofn = fullfile(eye_dir, [finfo(i).eye '.mj2']);
parfn = fullfile(mainpath,'matlab/eyepar',[finfo(i).eye '_fitinfo.mat']);

% chk_pupiltrack(videofn,parfn,5001:10000)
chk_pupiltrack(videofn,parfn)


% VR= filt_video(videofn);
% VR =read_subsample(VR, 1:1000);


% load('Y:\data_2photon\AWAKE\thy1_0129B3\05102016_thy1_0129B3_1\matlab\eyepar\0003_fitinfo.mat')
% par.out=fitout.out;
% par.fit = fitout.fitInfo;
% par.searchInfo = searchInfo
% parfn= 'Y:\data_2photon\AWAKE\thy1_0129B3\05102016_thy1_0129B3_1\matlab\eyepar\0002_fitinfo.mat'