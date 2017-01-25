clear all
% eye_dir = 'Y:\data_2photon\AWAKE\0822BRLF\11082015_0822BRLF\eye11082015'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015';
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF_2/eye11082015_2';
% eye_dir ='Y:\data_2photon\AWAKE\0705BRL\10112015_0705BRL_Nat_Grating\eye10112015'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0828BLF/11082015_0828BLF/eye11082015_PM_1'
% eye_dir ='Y:\data_2photon\AWAKE\0828BLF\11082015_0828BLF_2\eye11082015_PM_2'
% eye_dir = '/media/sdb_WD4T/data_2photon/AWAKE/thy1_0822BL/11042015_thy1_0822BL_2/eye11042015_2'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/thy1_0822BN/11052015_thy1_0822BN/eye11052015'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/thy1_0822BN/11052015_thy1_0822BN_2/eye11052015_2'

eye_dir = '/media/sdb_WD4T/data_2photon/AWAKE/thy1_0129B3/05102016_thy1_0129B3_1/eye0129B3_1'

recordnum = 2
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));

self= filt_video(videofn);
self =read_subsample(self, 1:10: self.NumberOfFrames);
heq = @(self) histeq_(self);
H=fspecial('gaussian',5,5);self.Xr
sF = @(self) filt(self,H);
self = addfiltstep(self,heq);
self = addfiltstep(self,sF);
applyfilt(self);  



%% pick up the threshold to find the pupil in the selected frame
baseimg = single(self.data(:,:,2));
par.parname={'Mean'};
par.op1={'>p'};
par.pmaps = single(self.data(:,:,2));
par.fnsave='THRtemp.mat';
validate_thresholds(baseimg,par);


%---------------
opt1.winsize=10;
opt1.thr = 5;
opt1.bdisp =false;

SEL.thr = opt1.thr ;
SEL.op ='>';
InitPar = gen_manual_pupilmask(self.data,1, SEL);
out1 = track_rawdata(self.data,opt1,InitPar);

% inx=150
% CM=out1.CM(:,inx);
% A = self.data(:,:,inx);
% pupil_radius=50
% center=[round(CM(2)) round(CM(1))];
% [pixel_list,~,bound] = disk_roi(255-A,center,pupil_radius,[1 1],0);
% msk = A;
% msk(bound)=200;
% figure; imagesc(msk)





%% second pass filter
Vw = self.Width;
Vh = self.Height;
CM=out1.CM;
Nframe = self.NumberOfFrames;
hSize=100;
Xr=floor([min(CM(1,:))-hSize max(CM(1,:))+hSize]);
Yr=floor([min(CM(2,:))-hSize max(CM(2,:))+hSize]);

Xr=self.Xr;
Yr=self.Yr
% - loading data
self =read_subsample(self, [],Xr,Yr);
% heq = @(self) histeq_(self);
% H=fspecial('gaussian',5,2);
% sF = @(self) filt(self,H);
% self = addfiltstep(self,heq);
% self = addfiltstep(self,sF);
applyfilt(self);  



iimg = 180;
tmpImg = self.data(:,:,iimg);
baseimg = single(tmpImg);
par.parname={'Mean'};
par.op1={'>p'};
par.pmaps = single(tmpImg);
par.fnsave='THRtemp.mat';
validate_thresholds(baseimg,par);

imgseg.hSize = hSize;
imgseg.Xr = self.Xr;
imgseg.Yr = self.Yr;


opt2.winsize=10;
opt2.thr =5;
opt2.bdisp =false;


out2.CM=zeros(2,Nframe);
out2.sPIX=cell(Nframe,1);

mCM = mean(out1.CM,2)-[Xr(1);Yr(1)];

iimg=1;
SEL.thr = opt2.thr;
SEL.op ='>';
InitPar = gen_manual_pupilmask(self.data(:,:,iimg),1, SEL);

chunksize =1000;
hf = @(data,selfr)track_rawdata(data,opt3,InitPar,CM0,2,selfr);
out2 = apply_fun_chunk (hf, self.data, chunksize, inxbadfit,fitout.out2);





% save fitInfo fitInfo -v7.3
fitout.opt1=opt1;
fitout.out1=out1;
fitout.imgseg = imgseg;
fitout.opt2=opt3;
fitout.out2=out3;
fitout.fitInfo = fitInfo;
fitout.error_frames = find(fail==1);

P=zeros(3,Nframe);
for ii=1:length(fitInfo)
    if isempty(fitInfo(ii).par)
        continue;
    end
    P(:,ii)=fitInfo(ii).par;
end


[locpath,sfn]= fileparts(videofn);
sfn = [sfn '_fitinfo.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');
mkdir(locpath)
save(fullfile(locpath,sfn),'fitout','-v7.3');



%% re-est on failed data with data loading



clear all
addpath(genpath('/home/slee/data/codes2P'));
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015';
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF_2/eye11082015_2';
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/thy1_0822BN/11052015_thy1_0822BN_2/eye11052015_2'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0828BLF/11082015_0828BLF/eye11082015_PM_1';
%eye_dir='/media/sdb_WD4T/data_2photon/AWAKE/0828BLF/11082015_0828BLF_2/eye11082015_PM_2'

eye_dir = '/media/sdb_WD4T/data_2photon/AWAKE/thy1_0822BL/11042015_thy1_0822BL_2/eye11042015_2'
recordnum = 3;
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));

[locpath,sfn]= fileparts(videofn);
sfn = [sfn '_fitinfo.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');

Modelfn = fullfile(locpath,sfn);
if exist(Modelfn,'file')
    load(Modelfn)
    fitInfo =fitout.fitInfo;
end

self= filt_video(videofn);


%----- loading data

Nframe = self.NumberOfFrames;

if exist('fitout','var')
    Xr = fitout.imgseg.Xr;
    Yr = fitout.imgseg.Yr;
else
    Xr = [1 self.Width];
    Yr = [1 self.Height];
end
self =read_subsample(self, [] ,Xr,Yr);



heq = @(self) histeq_(self);
H=fspecial('gaussian',5,2);
sF = @(self) filt(self,H);
self = addfiltstep(self,heq);
self = addfiltstep(self,sF);

applyfilt(self);  

% -----------
% out =fitout.out2;
% [inxbadfit, pars] =get_badfit(fitout.fitInfo,[25 80 15 5]);
% inxbadCM = unique(get_unstableCM(out.CM,1));
% inxfail0 = union(inxbadCM,inxbadfit);
% inxfail = inxfail0;
% inxgoodfit = setdiff(1:Nframe, inxfail0);
% CM0 = mean(out.CM(:,inxgoodfit),2);
% CM0 = repmat(CM0,[1 Nframe]);
% 
% figure; plot([out.CM(1,:)' pars(1,:)'],'.'); hold on; plot(inxfail0, out.CM(1,inxfail0),'k.')
% ylim([0 500])
% 


clear CM0
% out =fitout.out2;
% figure; plot(out.CM(1,:))

mImg = mean(self.data(151:350,101:550,:),3);
figure; imagesc(mImg)
[~,inx] = min(mImg(:));
[y0, x0]=ind2sub(size(mImg),inx);
x0 = x0 + 100;
y0 = y0 + 150;
CM0=[x0 y0]';

CM0 = repmat(CM0,[1 Nframe]);
%------------------------
inxfail = 1:Nframe;
PL=true;
thrs =[3 5 10 15 30 60];
[out, inxfails,outtmps] = reest_pupcms(out,self.data,CM0,inxfail,thrs,3,PL);
[out, inxfails,outtmps] = reest_pupcms(out,self.data,out.CM,inxfail,thrs,1,PL);

chunksize=1000;
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,150,6,selfr,1);
fitInfo1 = apply_fun_chunk (hf2, chunksize, fitInfo, self.data, inxfail,out.CM,out.sPIX);
[inxfail1, pars] =get_badfit(fitInfo1,[20 100 15 5]);
[~,inxfail1] = detect_spk(pars(3,:),5,10);
figure; plot([out.CM(1,:)' pars(1,:)'],'.'); hold on; plot(inxfail1, out.CM(1,inxfail1),'k.')
figure; plot([pars(3,10001:end)'],'.'); hold on; plot(inxfail1, pars(3,inxfail1),'k.')



%-this is particularly large pupil and largemove--------------------
[inxfail1, pars] =get_badfit(fitInfo1,[25 80 15 5]);
[~,inxfailx] = detect_spk(pars(3,:),5,10);
inxfail1 = union(inxfail1, inxfailx);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail1, pars(3,inxfail1),'k.')

thrs =[10 20 30 60];
[out, inxfails,outtmps2] = reest_pupcms(out,self.data,out.CM,inxfail1,thrs,1,PL);
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,100,3,selfr,3);
fitInfo2 = apply_fun_chunk (hf2, chunksize, fitInfo1, self.data, inxfail1,out.CM, out.sPIX);
[inxfail2, pars] =get_badfit(fitInfo2,[20 140 15 5]);
% figure; plot([out.CM(1,:)' pars(1,:)'],'.'); hold on; plot(inxfail1, pars(1,inxfail1),'k.')
% ylim([0 100])
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail2, pars(3,inxfail2),'k.')
ylim([0 100])

%------------
inxfail3 =[(1:741) (2810:2958) (3186:3359)  (5708:6100)];
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail3, pars(3,inxfail3),'k.')


thrs =[20 30];
[out, inxfails,outtmps2] = reest_pupcms(out,self.data,out.CM,inxfail3,thrs,1,PL);
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM0,sPIX,100,3,selfr,3);
fitInfo3 = apply_fun_chunk (hf2, chunksize, fitInfo2, self.data, inxfail3,out.CM, out.sPIX);
[inxfail3, pars] =get_badfit(fitInfo3,[20 140 15 5]);



% [~,inxfail4] = detect_spk(pars(1:2,:),5,10);
% [~,inxfailx] = detect_spk(pars(3,:),5,10);
% inxfail4 = union(inxfail4, inxfailx);
% inxfail4 = inxfail4(find(inxfail4<6500));
% figure; plot([pars(3,:)],'.'); hold on; plot(inxfail4, pars(3,inxfail4),'k.')
% thrs =[30 60];
% [out, inxfails,outtmps2] = reest_pupcms(out,self.data,CM0,inxfail4,thrs,1,PL);
% hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,100,3,selfr,3);
% fitInfo4 = apply_fun_chunk (hf2, chunksize, fitInfo3, self.data, inxfail4,out.CM, out.sPIX);
% [~, pars] =get_badfit(fitInfo4,[20 140 15 5]);
% figure; plot([pars(3,:)],'.'); hold on; plot(inxfail4, pars(3,inxfail4),'k.')

%%%
opt.winsize=10;
opt.thr = 30;
opt.bdisp = false;

hf = @(data,CM,selfr)track_rawdata2(data,opt,CM,6,selfr);
out = apply_fun_chunk (hf, chunksize,out,self.data,inxfail3, CM0);
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,300,3,selfr,3);
fitInfo3 = apply_fun_chunk (hf2, chunksize, fitInfo2, self.data, inxfail3,out.CM, out.sPIX);
[inxfail3, pars] =get_badfit(fitInfo3,[20 140 15 5]);
% inxfail3 = find(pars(3,:)<25 | pars(3,:)>120);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail3, pars(3,inxfail3),'k.')

 inxfail3 = find(pars(3,:)>100 );
 opt.thr=30
hf = @(data,CM,selfr)track_rawdata2(data,opt,CM,6,selfr);
out = apply_fun_chunk (hf, chunksize,out,self.data,inxfail3, out.CM);
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,300,3,selfr,3);
fitInfo3 = apply_fun_chunk (hf2, chunksize, fitInfo3, self.data, inxfail3,out.CM, out.sPIX);
[inxfail4, pars] =get_badfit(fitInfo3,[20 140 15 5]);
% inxfail3 = find(pars(3,:)<25 | pars(3,:)>120);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail4, pars(3,inxfail4),'k.')



%-----------------------------
opt.thr=30;
hf = @(data,CM,selfr)track_rawdata2(data,opt,CM,6,selfr);
out = apply_fun_chunk (hf, chunksize,out,self.data,inxfail3, out.CM);
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,300,3,selfr,3);
fitInfo4 = apply_fun_chunk (hf2, chunksize, fitInfo3, self.data, inxfail4,out.CM, out.sPIX);
[inxfail4, pars] =get_badfit(fitInfo4,[20 140 15 5]);

%----- manual correction


[~,inxfail6] = detect_spk(pars(3,:),5,10);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail6, pars(3,inxfail6),'k.')


opt.winsize=10;
opt.thr=60;
opt.bdisp=false;
hf = @(data,CM,selfr)track_rawdata2(data,opt,CM,6,selfr);
out2 = apply_fun_chunk (hf, chunksize,out,self.data,inxfail6, out.CM);

hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,300,3,selfr,3);
fitInfo6 = apply_fun_chunk (hf2, chunksize, fitInfo5, self.data, inxfail6,out.CM, out.sPIX);

[inxfail6, pars] =get_badfit(fitInfo6,[35 120 15 5]);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail6, pars(3,inxfail6),'k.')


% --- plot
iimg=inxfail3(10)
msk=self.data(:,:,iimg);
msk(out.sPIX{iimg})=255;
msk(fitInfo3(iimg).data)=255;
figure; imagesc(msk)


iimg=15400
msk=self.data(:,:,iimg);
% msk(out.sPIX{iimg})=255;
% msk(msk(:)<20)=255;
msk(fitInfo1(iimg).data)=255;
figure; imagesc(msk)

%---------------------------------

fitInfo = fitInfo1;
fitout.fitInfo = fitInfo;
fitout.out3 = out;
opt = fitout.opt2;
opt.thr =[3 5 10 15 20 30 60 80];
fitout.opt3=opt
[~, pars] =get_badfit(fitInfo,[25 120 15 5]);
inxfail = find(abs(fitout.out3.CM(1,:)-pars(1,:))>20);
fitout.error_frames = inxfail;

% save fitInfo fitInfo -v7.3

[locpath,sfn0]= fileparts(videofn);
sfn = [sfn0 '_fitinfo-2.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');
mkdir(locpath)
save(fullfile(locpath,sfn),'fitout','-v7.3');
imgseg =fitout.imgseg;
save(fullfile(locpath,[sfn0 '_smallfitinfo-2.mat']),'fitInfo','imgseg','-v7.3');


%% rerun est_puppar
% 
% %####################### rerun est_pupar
% 
% clear all
% 
% 
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF_2/eye11082015_2';
% recordnum = 5
% videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));
% 
% [locpath,sfn]= fileparts(videofn);
% sfn = [sfn '_fitinfo.mat'];
% locpath = fullfile(fileparts(locpath),'matlab','eyepar');
% load(fullfile(locpath,sfn))
% Xr = fitout.imgseg.Xr;
% Yr = fitout.imgseg.Yr;
% fitInfo =fitout.fitInfo;
% 
% self= filt_video(videofn);
% Nframe = self.NumberOfFrames;
% 
% %----- loading data
% 
% 
% 
% Xr = fitout.imgseg.Xr;
% Yr = fitout.imgseg.Yr;
% self =read_subsample(self, [],Xr,Yr);
% heq = @(self) histeq_(self);
% H=fspecial('gaussian',5,5);
% sF = @(self) filt(self,H);
% self = addfiltstep(self,heq);
% self = addfiltstep(self,sF);
% applyfilt(self);  
% 
% 
% %------cal
% out2 = fitout.out2;
% Nframe = self.NumberOfFrames;
% chunksize=1000;
% s = struct('par',[],'data',[]);
% fitInfo = repmat(s,Nframe,1);
% fail = zeros(Nframe,1);
% nchunk = ceil(Nframe/chunksize);
% for ichunk = 1 : nchunk
%     framelist = (ichunk-1)*chunksize + 1 : ichunk*chunksize;
%     if framelist(end)>Nframe
%         framelist = (ichunk-1)*chunksize + 1: Nframe;
%     end    
%     subdata = self.data(:,:,framelist);
%     CM=out2.CM(:,framelist);
%     sPIX = out2.sPIX(framelist);
%     [fitInfoi, faili] = est_puppar(subdata, CM, sPIX,100);
%     fitInfo(framelist) = fitInfoi(:);
%     fail(framelist) = faili;    
% end
% 
% 
% % save fitInfo fitInfo -v7.3
% 
% fitout.fitInfo = fitInfo;
% fitout.error_frames = find(fail==1);
% 
% 
% [locpath,sfn]= fileparts(videofn);
% sfn = [sfn '_fitinfo-2.mat'];
% locpath = fullfile(fileparts(locpath),'matlab','eyepar');
% mkdir(locpath)
% save(fullfile(locpath,sfn),'fitout','-v7.3');

