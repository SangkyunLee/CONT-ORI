function get_datawithPA(iexp_type, DATA_thr_str, comcont, ctm)




[~, ORI_list, ~, nses, seslist] =get_expinfo(iexp_type);


%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
fntmp = {'AN1-16','AN17-22','','',''};

dtype='Xsel';

%----------------
if ~exist('ctm','var')
ctm=0.6; 
end

cell_sel_method = 'UNION_CONTRSP'; 
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);


currentFolder = pwd;

fnsave = sprintf('PADEPDATA_%s-%s_ctm%0.2f.mat',fntmp{iexp_type},dtype,ctm); %L2norm -regularization
fullfnsav = fullfile(fileparts(currentFolder),'GRP_data',exp_type{iexp_type},DATA_thr_str, fnsave);


% population activity thresholds for high vs low population activity
PATHR.L=[0 0.5]';
PATHR.H=[0.5 1]';


scriptname = mfilename('fullpath');

CELLINX_SEL = cell(length(comcont),length(ORI_list),nses);
DPA = cell(length(comcont),length(ORI_list),nses);
%% ---------------------------------------------------
for ises = seslist 
  
    
    fndata1 = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    [D1]= loadData(data_path,fndata1);
    
    if iexp_type==5
        D1.events_ORI(D1.events_ORI(:)==-15)=-10;
    end

    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    for icomp = 1 : length(ORI_list)
 
        for icont = 1 : length(comcont)
            
            fprintf('ises: %d, icont:%d, icomp:%d\n',ises,icont, icomp);


            continfo = strsplit(comcont{icont},',');
            sevts{1} = str2double(continfo{1});
            sevts{2} = ORI_list(icomp);
            eORI = D1.events_ORI(:);
            if strcmp(continfo{2},'L')
                dtr = collect_subdata_PA(D1.(dtype), evt1, sevts,eORI,PATHR.L);
            elseif strcmp(continfo{2},'H')
                dtr = collect_subdata_PA(D1.(dtype), evt1, sevts,eORI,PATHR.H);
            end           
           
        
            DPA{icont,icomp, ises} = dtr;
            CELLINX_SEL{ises}=D1.cellinx_sel;
            save(fullfnsav,'DPA',...
                'CELLINX_SEL', 'PATHR','comcont',...
                'scriptname');
            pause(1);


        end % for icont  
    end %for icomp0 

end % ises


end
%-------------------------------

function varargout = loadData(data_path,varargin)
    ninput = length(varargin);
    varargout = cell(1,ninput);
    for i = 1: ninput
        fndata = varargin{i};
        fullfndata = fullfile(data_path,fndata);
        disp(['fndata= ' fullfndata]);
        varargout{i} = load(fullfndata);
    end
end


