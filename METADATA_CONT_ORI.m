%% list all experiments successfully performed
function [path_list, eyeinfo] =METADATA_CONT_ORI (exp_type,mainpath1, mainpath2,mainpath3,mainpath4)

switch exp_type
    case {'AWAKE','AWAKE_EYE'}
    
        % only data4 is ok
        path_list{1} = fullfile(mainpath3,'0705BL','09142015_0705BL');%spiral
        path_list{2} = fullfile(mainpath3,'0705BL','09142015_0705BL_2');%sprial
        path_list{3} = fullfile(mainpath3,'0705BR','09182015_0705BR_2');%reso
        path_list{4} = fullfile(mainpath3,'0705BN','09152015_0705BN');%spiral

        path_list{5} = fullfile(mainpath3,'0705BR','09212015_0705BR_Nat_Grating'); %reso

        path_list{6} = fullfile(mainpath3,'0705BL','10072015_0705BL_Nat_Grating');%spiral
        path_list{7} = fullfile(mainpath3,'0705BR','10042015_0705BR_Nat_Grating');%spiral
        path_list{8} = fullfile(mainpath3,'0705BR','10122015_0705BR_Nat_Grating');%spiral
        path_list{9} = fullfile(mainpath3,'0705BRL','10032015_0705BRL_Nat_Grating');%spiral
        path_list{10} = fullfile(mainpath3,'0705BRL','10112015_0705BRL_Nat_Grating');%spiral

        %----------------------
        path_list{11} = fullfile(mainpath3,'0822BRLF','11082015_0822BRLF');
        eyeinfo.l2r(11) = 588;
        eyeinfo.MPD(11) = 80*2;
        path_list{12} = fullfile(mainpath3,'0822BRLF','11082015_0822BRLF_2');

        path_list{13} = fullfile(mainpath3,'0828BLF','11082015_0828BLF');
        path_list{14} = fullfile(mainpath3,'0828BLF','11082015_0828BLF_2');

        path_list{15} = fullfile(mainpath3,'thy1_0822BL','11042015_thy1_0822BL_2');   
        path_list{16} = fullfile(mainpath3,'thy1_0822BN','11052015_thy1_0822BN');
        path_list{17} = fullfile(mainpath3,'thy1_0822BN','11052015_thy1_0822BN_2');

        %------------------------- I have to process
        % path_list{18}= fullfile(mainpath3,'0705BRL','100112015_0705BRL');%spiral
        
        
        path_list{21} = fullfile(mainpath3,'thy1_0129B2','04192016_thy1_0129B2_3');
        path_list{22} = fullfile(mainpath3,'thy1_0129B1','04192016_thy1_0129B1_3');

        
        
        % -------------------------
        path_list{23} = fullfile(mainpath3,'thy1_0129B1','04272016_thy1_0129B1_4');
        path_list{24} = fullfile(mainpath3,'thy1_0129B1','04272016_thy1_0129B1_5'); % two sessions, cell not identified yet
        path_list{25} = fullfile(mainpath3,'thy1_0129B1','05122016_thy1_0129B1_9'); % two sessions
        path_list{26} = fullfile(mainpath3,'thy1_0129B2','04282016_thy1_0129B2_4');
        path_list{27} = fullfile(mainpath3,'thy1_0129B2','04282016_thy1_0129B2_5'); % two sessions,
        path_list{28} = fullfile(mainpath3,'thy1_0129B2','05062016_thy1_0129B2_6'); % two sessions

        path_list{29} = fullfile(mainpath3,'thy1_0129B3','05102016_thy1_0129B3_1');% two sessions
        path_list{30} = fullfile(mainpath3,'thy1_0129B3','05112016_thy1_0129B3_2'); % two sessions, 
        path_list{31} = fullfile(mainpath3,'thy1_0129B3','05182016_thy1_0129B3_3'); % one session(0:30:330)

        path_list{32} = fullfile(mainpath3,'022216BL','05132016_022216BL_2');% two sessions
        path_list{33} = fullfile(mainpath3,'022216BL','05152016_022216BL_3');% two sessions
        path_list{34} = fullfile(mainpath3,'022216BL','05182016_022216BL_4');% two sessions (0:30:330)

        path_list{35} = fullfile(mainpath3,'022216BN','05132016_022216BN_2');% two sessions
        path_list{36} = fullfile(mainpath3,'022216BN','05152016_022216BN_3');% two sessions 
        %---------
        path_list{37} = fullfile(mainpath3,'022016BN','05162016_022016BN_1');% two sessions
        path_list{38} = fullfile(mainpath3,'022016BN','05182016_022016BN_2');% two sessions 
        path_list{39} = fullfile(mainpath3,'022016BRF','05172016_022016BRF');% two sessions
        path_list{40} = fullfile(mainpath3,'thy1_032016BR','05172016_Thy1_032016BR_2');% two sessions

        
        eyeinfo.l2r(23) = 390;
        eyeinfo.MPD(23) = 60*2;
        
        eyeinfo.l2r(25) = 390;
        eyeinfo.MPD(25) = 60*2;
        
        eyeinfo.l2r(26) = 420;
        eyeinfo.MPD(26) = 62*2;
        
        eyeinfo.l2r(27) = 410;
        eyeinfo.MPD(27) = 90*2;
        
        eyeinfo.l2r(28) = NaN;
        eyeinfo.MPD(28) = NaN;
        
        eyeinfo.l2r(29) = 390;
        eyeinfo.MPD(29) = 60*2;
        
        
        eyeinfo.l2r(30) = 325;
        eyeinfo.MPD(30) = 41*2;
        
        eyeinfo.l2r(32) = 410;
        eyeinfo.MPD(32) = 81*2;
        
        eyeinfo.l2r(33) = 440;
        eyeinfo.MPD(33) = 87*2;
        
        eyeinfo.l2r(35) = NaN;
        eyeinfo.MPD(35) = NaN;
        
        eyeinfo.l2r(36) = 420;
        eyeinfo.MPD(36) = 87*2;
        
        eyeinfo.l2r(40) = 420;
        eyeinfo.MPD(40) = 84*2;
        %varout{1}=eyeinfo;

    
	case {'AN_0TO150'}
        % 1, 7 are best, 2, 5 ok
    
        path_list{1} = fullfile(mainpath2,'Thy1_4P3_05042015BN','07102015_thy1_05042015BN');
        path_list{2} = fullfile(mainpath2,'Thy1_4P3_05042015BN','07102015_thy1_05042015BN_2');
        path_list{3} = fullfile(mainpath2,'Thy1_4P3_05042015BR','07092015_thy1_05042015BR_1');
        path_list{4} = fullfile(mainpath2,'Thy1_4P3_05042015BR','07122015_thy1_05042015BR_2');
        path_list{5} = fullfile(mainpath2,'Thy1_4P3_05102015BL','07142015_thy1_0510BL');
        path_list{6} = fullfile(mainpath2,'Thy1_4P3_05102015BL','07142015_thy1_0510BL_2');
        path_list{7}  = fullfile(mainpath1,'150507BRL','07222015_0507BRL');
         %  I have to process data
        %fn_list{8}  = fullfile(mainpath1,'150507BRL','08142015_0507BRL_Nat3_Grating')
        
    case {'AN'}
        %------ main BCM
        path_list{1} = fullfile(mainpath1,'140226 BL','20140327');
        path_list{2} = fullfile(mainpath1,'R140123 BR','20140403');
        path_list{3} = fullfile(mainpath1,'140226 BL','20140325');
        path_list{4} = fullfile(mainpath1,'140226 BL','20140329');
        path_list{5} = fullfile(mainpath1,'140226 BL','20140401');
        path_list{6} = fullfile(mainpath1,'140328BL','140425');
        path_list{7} = fullfile(mainpath1,'140505BN','140611');
        path_list{8} = fullfile(mainpath1,'140506BR','140527');

        % I have to process data
        %path_list{9} = fullfile(mainpath1,'140603BRM','140622')
        %path_list{10} = fullfile(mainpath1,'140618BRM','140717')
        %------ left eye not blocked
        path_list{11} = fullfile(mainpath2,'110315BN','02252016_110315BN');%spiral512 
        path_list{12} = fullfile(mainpath2,'110315BN','02252016_110315BN_2');%sprial512
        
        path_list{13} = fullfile(mainpath4,'thy1_012916B1','04172016_thy1_0129');%spiral512, High lum - bad light shielding
        path_list{14} = fullfile(mainpath4,'thy1_012916B1','04172016_thy1_0129_2');%sprial256x256, high lum - bad light shielding
        % hereafter left eye blocked
        path_list{15} = fullfile(mainpath4,'thy1_012916B2','04192016_thy1_0129B2');%sprial-lefteyeblock
        path_list{16} = fullfile(mainpath4,'thy1_012916B2','04192016_thy1_0129B2_2');%sprial -lefteyeblock
        
        path_list{17} = fullfile(mainpath4,'thy1_012916B2','05072016_thy1_0129B2_7');%sprial -lefteyeblock, two-sessions
        path_list{18} = fullfile(mainpath4,'thy1_012916B1','05062016_thy1_0129B1_8');%sprial -lefteyeblock, two-sessions
        path_list{19} = fullfile(mainpath4,'thy1_030416BN','05052016_thy1_030416BN');%sprial -lefteyeblock, two-sessions
        path_list{20} = fullfile(mainpath4,'thy1_032016BR','05102016_thy1_032016BR');%sprial -lefteyeblock, two-sessions
        path_list{21} = fullfile(mainpath4,'022216BL','05112016_022216BL_1');%sprial -lefteyeblock, two-sessions
        path_list{22} = fullfile(mainpath4,'022216BN','05122016_022216BN_1');%sprial -lefteyeblock, two-sessions


        

        
        %----------- ORI=0:30:330 and contrast=[40 100]
    case {'AN_FULLORI'}
        path_list{1} = fullfile(mainpath4,'thy1_012916B1','04292016_thy1_0129B1_6');%sprial -lefteyeblock, medium lum
        path_list{2} = fullfile(mainpath4,'thy1_012916B1','05032016_thy1_0129B1_7');%sprial -lefteyeblock, mainBCM
        
  
        
    otherwise
        error('not specified exp_type');
end