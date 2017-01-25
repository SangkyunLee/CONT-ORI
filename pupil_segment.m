clear all

eye_dir = '/media/sdb_WD4T/data_2photon/AWAKE/thy1_0129B3/05102016_thy1_0129B3_1/eye0129B3_1'

recordnum = 2
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));

self= filt_video(videofn);
self =read_subsample(self, 1:10: self.NumberOfFrames);
heq = @(self) histeq_(self);
H=fspecial('gaussian',5,5);
sF = @(self) filt(self,H);
self = addfiltstep(self,heq);
self = addfiltstep(self,sF);
applyfilt(self);  



%% pick up the threshold to find the pupil in the selected frame
baseimg = single(self.data(:,:,i));
par.parname={'Mean'};
par.op1={'>p'};
par.pmaps = single(self.data(:,:,i));
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


%% load all data

Vw = self.Width;
Vh = self.Height;
CM=out1.CM;
Nframe = self.NumberOfFrames;
hSize=100;
Xr=floor([min(CM(1,:))-hSize max(CM(1,:))+hSize]);
Yr=floor([min(CM(2,:))-hSize max(CM(2,:))+hSize]);

% Xr=self.Xr;
% Yr=self.Yr
% - loading data
self =read_subsample(self, [],Xr,Yr);
applyfilt(self);  



CM0 = mean(out1.CM,2);
D = sqrt(sum(bsxfun(@minus,out1.CM,CM0).^2,1));
thr = 1*std(D);
inx = find((D-mean(D))<thr);
CM0 = mean(out1.CM(:,inx),2);
inxbadfit = 1:Nframe;
hf = @(data,CM0,selfr)track_rawdata2_PL(data,opt1,CM0,2,selfr);
out2.sPIX=cell(1,Nframe);
out2.CM=zeros(2,Nframe);

chunksize =1000;

out2 = apply_fun_chunk (hf, chunksize,out2,self.data,inxbadfit, CM0);


%%
ref_fr = 1:100;
CM0=mean(out1.CM(:,ref_fr),2);
D = sqrt(sum(bsxfun(@minus,out2.CM,CM0).^2,1));

thr = 1*std(D);
inxfail = find((D-mean(D))>thr);

thrs =[3 5 10 15 30 60];
[out3, inxfails,outtmps] = reest_pupcms(out2,self.data,CM0,inxfail,thrs,1,1);
[out3, inxfails,outtmps] = reest_pupcms(out3,self.data,CM0,inxfails{1},thrs,1,0);
[out3, inxfails,outtmps] = reest_pupcms(out3,self.data,CM0,inxfails{1},thrs,1,0);


a = track_rawdata2_PL(self.data,opt1,CM0*ones(1,Nframe),2,5512);

chunksize=1000;
hf2 = @(data,CM,sPIX,selfr) est_puppar(data,CM,sPIX,150,6,selfr,1);
fitInfo1 = apply_fun_chunk (hf2, chunksize, [], self.data, 1:Nframe,out3.CM,out3.sPIX);

[inxfail1, pars] =get_badfit(fitInfo1,[10 100 10 5]);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail1, pars(3,inxfail1),'k.')
ylim([0 100])

inxbadfit=5610:6041;
thrs =[10 30 60];
[out3, inxfails,outtmps] = reest_pupcms(out3,self.data,out3.CM,inxbadfit,thrs,1,0);
fitInfo1 = apply_fun_chunk (hf2, chunksize, fitInfo1, self.data, inxbadfit,out3.CM,out3.sPIX);
[inxfail1, pars] =get_badfit(fitInfo1,[10 100 10 5]);
figure; plot([pars(3,:)],'.'); hold on; plot(inxfail1, pars(3,inxfail1),'k.')
ylim([0 100])


fitInfo = fitInfo1;
fitout.fitInfo = fitInfo1;
fitout.out = out3;

opt1.thr =[3 5 10 15 30 60];
fitout.opt=opt1;
[~, pars] =get_badfit(fitInfo1,[25 120 15 5]);
inxfail = find(abs(fitout.out.CM(1,:)-pars(1,:))>30);



fitout.error_frames = [];
figure; plot([pars(1,:)],'.'); hold on; plot(inxfail, pars(1,inxfail),'k.')

% save fitInfo fitInfo -v7.3

[locpath,sfn0]= fileparts(videofn);
sfn = [sfn0 '_fitinfo.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');
mkdir(locpath)
save(fullfile(locpath,sfn),'fitout','-v7.3');

save(fullfile(locpath,[sfn0 '_smallfitinfo-2.mat']),'fitInfo','-v7.3');



