clear all
close all

addpath('/home/slee/data/codes2P');
eye_dir = 'Y:\data_2photon\AWAKE\0822BRLF\11082015_0822BRLF_2\eye11082015_2'
eye_dir = 'Y:\data_2photon\AWAKE\0828BLF\11082015_0828BLF\eye11082015_PM_1'
% eye_dir ='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015';
% since recordnum=4, most reliable 
recordnum=5;
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));

[locpath,sfn]= fileparts(videofn);
sfn = [sfn '_smallfitinfo-2.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');
M=load(fullfile(locpath,sfn));
fitInfo = M.fitInfo;
Xr = M.imgseg.Xr;
Yr = M.imgseg.Yr;
% fitInfo =fitout.fitInfo;

self= filt_video(videofn);
self =read_subsample(self, [],Xr,Yr);
Nframe = self.NumberOfFrames;


[~, pars] =get_badfit(fitInfo,[25 80 15 5]);
figure; plot([ pars(3,:)'],'.');


framelist = 1:self.NumberOfFrames;


piterr=plotfit(fullfile(locpath,[sfn(1:end-4) '.mp4']), self.data,framelist, fitInfo );

%%
% 
% clear all
% 
% 
% eye_dir = 'Y:\data_2photon\AWAKE\0822BRLF\11082015_0822BRLF_2\eye11082015_2'
% recordnum=2
% videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));
% 
% [locpath,sfn]= fileparts(videofn);
% sfn = [sfn '_fitinfo-2.mat'];
% locpath = fullfile(fileparts(locpath),'matlab','eyepar');
% load(fullfile(locpath,sfn))
% fitInfo =fitout.fitInfo;
% 
% self= filt_video(videofn);
% 
% 
% %----- checking the fitting
% Nframe = self.NumberOfFrames;
% framelist = fitout.error_frames;
% PARs= zeros(3,Nframe);
% 
% for ii=1:Nframe
%     if~isempty(fitInfo(ii).par)        
%         PARs(:,ii)=fitInfo(ii).par;
%     end
% end
% inx_badfit = find(PARs(3,:)>90);
% framelist = union(framelist,inx_badfit);
% 
% %----- loading data
% Xr = fitout.imgseg.Xr;
% Yr = fitout.imgseg.Yr;
% self =read_subsample(self, framelist ,Xr,Yr);
% 
% 
% 
% heq = @(self) histeq_(self);
% H=fspecial('gaussian',5,5);
% sF = @(self) filt(self,H);
% self = addfiltstep(self,heq);
% self = addfiltstep(self,sF);
% 
% applyfilt(self);  
% 
% % --- estimate
% iimg = 10;
% tmpImg = self.data(:,:,iimg);
% baseimg = single(tmpImg);
% par.parname={'Mean'};
% par.op1={'>p'};
% par.pmaps = single(tmpImg);
% par.fnsave='THRtemp.mat';
% validate_thresholds(baseimg,par);
% 
% fldnames = fields(fitout);
% inxout=0;
% for ifn = 1: length(fldnames)
%     if strfind(fldnames{ifn},'out')
%         inxout = inxout+1;
%     end
% end
% 
% SEL.thr = 4;
% SEL.op ='>';
% InitPar = gen_manual_pupilmask(self.data,10, SEL);
% eval(sprintf('opt%d = fitout.opt%d;',inxout+1,inxout));
% eval(sprintf('opt%d.thr =SEL.thr;', inxout+1));
% eval(sprintf('out%d = track_rawdata(self.data,opt%d,InitPar);',inxout+1,inxout+1));  
% eval(sprintf('CM=out%d.CM;',inxout+1));
% eval(sprintf('sPIX = out%d.sPIX;',inxout+1))
% [fitInfoi, faili] = est_puppar(self.data, CM, sPIX,100);
% fitInfo(framelist) = fitInfoi(:);
% fail(framelist) = faili;    
% 
% eval(sprintf('fitout.opt%d =opt%d;',inxout+1,inxout+1))
% eval(sprintf('fitout.out%d =out%d;',inxout+1,inxout+1))
% fitout.fitInfo = fitInfo;
% fitout.error_frames = find(fail==1);
% 
% 
% % save fitInfo fitInfo -v7.3
% 
% [locpath,sfn]= fileparts(videofn);
% sfn = [sfn '_fitinfo-2.mat'];
% locpath = fullfile(fileparts(locpath),'matlab','eyepar');
% mkdir(locpath)
% save(fullfile(locpath,sfn),'fitout','-v7.3');
% 
% 
% 
% 
% 
% 
% %%
% figure;
%  for i = 4140
%     iimg = i;
%     self =read_subsample(self, i,Xr,Yr);
%     A = self.data;
%  
%     A1 =histeq(uint8(A));
% 
%     boundary = fitInfo(iimg).data;
%     P = fitInfo(iimg).par;
%     if ~isempty(P)
%         t = linspace(0, 2*pi, 100);
%        
%             msk = A1;
%             msk(boundary)=0;
%             imagesc(msk)
%             hold on;
%             plot(P(1)  + P(3)  * cos(t), P(2)  + P(3) * sin(t), 'r')
%             axis equal
%             title(sprintf('%d',iimg));
%             pause(0.002);
%        
%     end
% end
% 
% %%
% 
% framelist=7536 %4960
% 
% out2 = fitout.out2
% 
% self =read_subsample(self, framelist,Xr,Yr);
% heq = @(self) histeq_(self);
% H=fspecial('gaussian',5,5);
% sF = @(self) filt(self,H);
% self = addfiltstep(self,heq);
% self = addfiltstep(self,sF);
% applyfilt(self);  
% 
% A = self.data(:,:,1);
% CM=out2.CM(:,framelist);
% % sPIX = out2.sPIX(framelist);
% % msk = A;
% % msk(sPIX{1})=100;
% % figure; imagesc(msk)
% 
% 
% center=[round(CM(2)) round(CM(1))];
% %[pixel_list,~,bound] = disk_roi(255-A,center,pupil_radius,[1 1],0);
% [pixel_list,boundary_list] = search_pupbound(255-A,center,100,'threshold');
% msk = A;
% msk(boundary_list)=100;
% figure; imagesc(msk)
% 
%     
% %%
% full_radius=100;
% 
% ntheta=90;
% nsample=pupil_radius;
% t=(1:ntheta)/ntheta*2*pi;
% x=round(center(2));
% y=round(center(1));
% lprofile=zeros(nsample,ntheta);
% lprofile_mod=zeros(nsample,ntheta);
% 
% r_first_min=zeros(1,ntheta);
% r_threshold_cross=zeros(1,ntheta);
% for i=1:ntheta
%     f=improfile(255-A,[x,x+full_radius*cos(t(i))],[y,y+full_radius*sin(t(i))],nsample,'bilinear')';    
%     lprofile(:,i)=f;%-[zeros(1,nsample/3),(1:(nsample*2/3))/nsample*max(f)/2];
%     
%     r_threshold_cross(i)=find_first_cross(f,0.8);       
%     lprofile_mod(:,i)=-1*abs((1:nsample)-r_threshold_cross(i))+nsample;
% end
% 
% msk = A;
% % msk(sPIX{1})=100;
% 
% 
% [pixel_list,~,bound] = disk_roi(255-A,center,39,[1 1],0);
% 
% msk(bound)=200;
% figure; imagesc(msk)