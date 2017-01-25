

clear all
addpath(genpath('/home/slee/data/codes2P'));
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015';
eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF_2/eye11082015_2';
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/thy1_0822BN/11052015_thy1_0822BN_2/eye11052015_2'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0828BLF/11082015_0828BLF/eye11082015_PM_1';
recordnum = 2;
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));

[locpath,sfn]= fileparts(videofn);
sfn = [sfn '_fitinfo.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');
load(fullfile(locpath,sfn))
fitInfo =fitout.fitInfo;

self= filt_video(videofn);


%----- loading data
Nframe = self.NumberOfFrames;
framelist = 1:Nframe;


PARs= zeros(3,Nframe);
list=zeros(Nframe,1);
for ii=1:Nframe
    if isempty(fitInfo(ii).par)
        list(ii)=1;
    else
        PARs(:,ii)=fitInfo(ii).par;
    end
end
df_PARs=[0 diff(PARs(3,:));diff(PARs(3,:))  0];
dist = sqrt(sum(PARs(1:2,:).^2,1));
rdist =dist -mean(dist);
inxoutCM = find(rdist>2*std(rdist));

inxP=abs(df_PARs(1,:))>20 & (df_PARs(1,:).*df_PARs(2,:))<0;
inxbadfit = find((PARs(3,:)<10 | PARs(3,:)>60) | inxP | rdist>2*std(rdist) );

framelist =union(framelist,inxbadfit);

Xr = fitout.imgseg.Xr;
Yr = fitout.imgseg.Yr;
self =read_subsample(self, framelist ,Xr,Yr);


dSize =size(self.data);

subdata=self.data(51:250,131:400,1:10:end);
subdata=imresize(subdata,0.5);
dSize2 = size(subdata);
X= single(reshape(subdata,[dSize2(1)*dSize2(2) dSize2(3)]));
X=bsxfun(@minus,X, mean(X,2));
[U, S, V]=svd(X);
figure; imagesc(mean(subdata,3))
