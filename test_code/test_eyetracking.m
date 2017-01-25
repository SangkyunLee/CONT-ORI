clear all
eye_dir = 'X:\Awake\EYE\eye0129B2_6'
% eye_dir='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015'
recordnum=1
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));
Vreader = VideoReader(videofn);

Nframe = Vreader.NumberOfFrames;
Frlist = 1 : 10 :Nframe;
j = 1;
for k = 1 : Nframe
    sMOV1 = read(Vreader,k);
    if mod(k,10)==1,
        if k==1,
            sMOV = uint8(zeros(size(sMOV1,1),size(sMOV1,2),length(Frlist)));
        end
        sMOV(:,:,j) = histeq(sMOV1);
        j = j+1;
    end
end

%% pick up the threshold to find the pupil in the selected frame
baseimg = single(sMOV(:,:,1));
par.parname={'Mean'};
par.op1={'>p'};
par.pmaps = single(sMOV(:,:,1));
par.fnsave='THRtemp.mat';
validate_thresholds(baseimg,par);


%---------------
opt1.winsize=10;
opt1.thr =18;
opt1.bdisp =false;

out1 = track_rawdata(sMOV,opt1);



%% second pass filter
Vw = Vreader.Width;
Vh = Vreader.Height;
CM=out1.CM;
Nframe = Vreader.NumberOfFrames;
hSize=100;
Xr=floor([min(CM(1,:))-hSize max(CM(1,:))+hSize]);
Yr=floor([min(CM(2,:))-hSize max(CM(2,:))+hSize]);
if Xr(1)<1, Xr(1)=1; end
if Xr(2)>Vw, Xr(2)=Vw; end
if Yr(1)<1, Yr(1)=1; end
if Yr(2)>Vh, Yr(2)=Vh; end

H=fspecial('gaussian',5,5);


iimg = 1;
IMG = histeq(read(Vreader, iimg));
tmpImg= IMG(Yr(1):Yr(2),Xr(1):Xr(2));
tmpImg=imfilter(tmpImg,H,'same');

baseimg = single(tmpImg(:,:,1));
par.parname={'Mean'};
par.op1={'>p'};
par.pmaps = single(tmpImg(:,:,1));
par.fnsave='THRtemp.mat';
validate_thresholds(baseimg,par);


chunksize=1000;
nchunk = ceil(Nframe/chunksize);
for ichunk = 1 : nchunk
    framelist = [(ichunk-1)*chunksize + 1 ichunk*chunksize];
    if framelist(end)>Nframe
        framelist = [(ichunk-1)*chunksize + 1 Nframe];
    end
    IMG = squeeze(read(Vreader, framelist));
    
    for ii = 1 : size(IMG,3)         
        tmpImg= histeq(IMG(Yr(1):Yr(2),Xr(1):Xr(2),ii));
        if ii ==1 && ichunk ==1,    
            fImg=uint8(zeros(size(tmpImg,1),size(tmpImg,2),Nframe));
        end
        fImg(:,:,iimg)=imfilter(tmpImg,H,'same');    
        iimg = iimg+1;
    end
end

imgseg.hSize = 100;
imgseg.Xr = Xr;
imgseg.Yr = Yr;


opt2.winsize=10;
opt2.thr =25;
opt2.bdisp =false;

out2 = track_rawdata(fImg,opt2);


% save out out2 fImg -v7.3



s = struct('par',[],'data',[]);
fitInfo = repmat(s,Nframe,1);

% figure;
H=fspecial('gaussian',5,5);
error_list = zeros(Nframe,1);
parfor iimg = 1:Nframe   
    CM = out2.CM(:,iimg);
    A = fImg(:,:,iimg);
    A1 = 255- A;
%     
%     A = read(Vreader,iimg);
%     A =A((Yr(1):Yr(2))-20,Xr(1):Xr(2));
%     A = histeq(A);
%     A1 =255- imfilter(A,H,'same'); 

    map =edge(A1);
    center=[round(CM(2)) round(CM(1))];
    pupil_radius=100;
    try
        [pixel_list,~,boundary]=disk_roi(A1,center,pupil_radius,[1 1],0);
   
%     msk = A1;
%     msk(boundary)=0;
%     imagesc(msk)
    
        msk=zeros(size(A));
        msk(pixel_list)=1;    
        [X0, Y0]=cmass(msk);
        pupil_radius=100;
    
        [pixel_list,~,boundary]=disk_roi(A1,[Y0, X0],pupil_radius,[1 1],0);
    
%     msk = A1;
%     msk(boundary)=0;
%     figure;imagesc(msk)
    

        pixel_list =union(pixel_list, out2.sPIX{iimg});    

        msk = zeros(size(A1));
        msk(pixel_list)=1;

        bw=bwperim(msk,4);
        boundary=find(bw(:));

        [y, x]= ind2sub(size(A1), boundary);

        [z, r, residual] = fitcircle([x y]');
        fitInfo(iimg).par=[z;r];
        fitInfo(iimg).data = boundary;
    catch
        error_list(iimg)=1;
    end
end
% save fitInfo fitInfo -v7.3
fitout.opt1=opt1;
fitout.out1=out1;
fitout.imgseg = imgseg;
fitout.opt2=opt2;
fitout.out2=out2;
fitout.fitInfo = fitInfo;
fitout.error_frames = find(error_list==1);


[locpath,sfn]= fileparts(videofn);
sfn = [sfn '_fitinfo.mat'];
locpath = fullfile(fileparts(locpath),'matlab','eyepar');
mkdir(locpath)
save(fullfile(locpath,sfn),'fitout','-v7.3');


%%
% figure;
% H=fspecial('gaussian',5,5);
% for iimg = 1:Nframe
%     CM = out2.CM(:,iimg);
%     A = fImg(:,:,iimg);
%     A1 = 255- A;
%     boundary = fitInfo(iimg).data;
%     P = fitInfo(iimg).par;
%     
%     t = linspace(0, 2*pi, 100);
%     if mod(iimg,3)==1.
%         msk = A1;
%         msk(boundary)=0;
%         imagesc(msk)
%         hold on;
%         plot(P(1)  + P(3)  * cos(t), P(2)  + P(3) * sin(t), 'r')
%         axis equal
%         title(sprintf('%d',iimg));
%         pause(0.002);
%     end
% end




%%


% iimg = 20;
% I =fImg(:,:,iimg);
% [BW th]=edge(I);
% figure; imagesc(BW)
% 
% figure; imagesc(255-I)
% 
% A = zeros(size(fImg,1),size(fImg,2));
% A(out2.sPIX{iimg})=1;
% A(round(out2.CM(2,iimg)),round(out2.CM(1,iimg)))=2;
% figure; imagesc(A)
% 
% 
% 
% 
% baseimg = single(fImg(:,:,iimg));
% par.parname={'Mean'};
% par.op1={'>p'};
% par.pmaps = single(baseimg(:,:,1));
% par.fnsave='THRtemp.mat';
% validate_thresholds(baseimg,par);

%%

% hfig = figure;
% [dy, dx,~]=size(fImg);
% left=0.1;
% if dy>=dx,
%     fheight = 0.7;
%     fwidth = fheight*dx/dy;
% else
%     fwidth = 0.7;
%     fheight = fwidth*dy/dx;
% end
% fbottom = 1 - (fheight+0.1);
% hax=axes('Position',[left fbottom fwidth fheight],'Parent',hfig);
% for iimg = 1 : Nframe
%     baseimg = fImg(:,:,iimg);
%     baseimg = repmat(single(baseimg), [1 1 3]);
%     fmap = zeros(size(baseimg));
%     fmap(out2.sPIX{iimg})=1;
%     overlayImage(baseimg,fmap, hax,'Continuous');
%     title(sprintf('%d',iimg));
%     pause(0.002);
% end

%%

% figure; 
% for iimg = 1 : Nframe
%     
%     IMG = histeq(read(Vreader, iimg));
%     tmpImg= IMG(Yr(1):Yr(2),Xr(1):Xr(2));
%     if iimg ==1,
%        % PImg=uint8(zeros(size(tmpImg,1),size(tmpImg,2),Nframe));
%         fImg=uint8(zeros(size(tmpImg,1),size(tmpImg,2),Nframe));
%     end
%     %PImg(:,:,iimg)=tmpImg;
%     fImg(:,:,iimg)=imfilter(tmpImg,H,'same');
% %     subplot(121);
% %     imagesc((PImg(:,:,iimg)));axis image;
% %     subplot(122);
% %     imagesc((fImg(:,:,iimg)));axis image;
% %     pause(0.002);
%     
% end

