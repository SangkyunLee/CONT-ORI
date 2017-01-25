function session = get_scanspersession(exp_type)

switch exp_type
    
    case 'AWAKE'
        session{10,1}=18:20; session{10,2}=[21 22 24];
        session{11,1}=2:4; session{11,2}=5:7;
        session{13,1}=4:6; session{13,2}=7:9;
        session{15,1}=2:4; session{15,2}=5:7;
        session{21,1}=1:3; session{21,2}=4:7;
        session{22,1}=1:3; session{22,2}=4:7;
        %-----------------------------------------------------------
        % session 26, 27,28  look not very good response time course
        
        
        session{23,1}=[1 2 ]; session{23,2}=[4 5 6];
        %session{23,1}=[1 2 4 5 6]; session 23 is originally 1 session but
        % 
        % scan 5: eyemovie frozen, lots of movement

        session{25,1} = 1:3; session{25,2} = 5:7;
        % scan4: lots of movement and eyemovie frozen, scan 6: eye closed

        session{26,1} = 1:3;session{26,2} = 4:5;
        %session{26,1} = 1:5;

        session{27,1} = 2:3;session{27,2} = 4:5;
        %session{27,1} = 2:5;        
        %originally two sessions acquired but the second session was discard due to bad eye

        session{28,1} = 1:4; session{28,2} = [5 7 8];

        session{29,1} = 1:3; session{29,2} = 4:6;

        session{30,1} = 1:4; session{30,2} = 5:8;
        %scan3: eyefrozen, 2nd session not good.

        session{32,1} = 1:4; session{32,2} = 5:7;

        session{33,1} = 1:4; session{33,2} = 5:8;
        %scan3: eyemovie frozen in the middle of scan

        session{35,1} = 1:4; session{35,2} = 6:9; 
        %6:lots of movement

        session{36,1} = 1:3; session{36,2} = 5:9;

        session{37,1} = 1:4; session{37,2} = 5:8; 
        % scan 7:  not good eye,

        session{40,1} = 1:3; session{40,2} = 4:6; 
        % scan 6: pupil contracted a lot



    case 'AN'
        session{1,1} = 2:6; session{1,2} = 7:11; 
        session{2,1} = 4:8; session{2,2} = 9:13; 
        session{3,1} = 3:5; session{3,2} = 6:8; 
        session{4,1} = 2:7; session{4,2} = 8:13; 
        session{5,1} = 3:7; session{5,2} = 8:12 ; 
        session{6,1} = 2:7; session{6,2} = [8 (10:14)]; 
        session{7,1} = 6:10; session{7,2} = 11:15; 
        session{8,1} = 4:9; session{8,2} = 10:14; 
        
        session{11,1} =2:3; session{11,2} = 4:5; 
        session{12,1} =1:3; session{12,2} = 4:5; 
        
        session{15,1} =1:3; session{15,2} = 4:5; 
        session{16,1} =1:3; session{16,2} = 4:5;
        %---------------------------------        
        session{17,1} = 1:3; session{17,2} = 4:6; 
        session{18,1} = 1:3; session{18,2} = 4:6; 
        session{19,1} = 1:3; session{19,2} = 4:7; 
        session{20,1} = 1:3; session{20,2} = 4:6; 
        session{21,1} = 1:3; session{21,2} = 4:6; 
        session{22,1} = 1:3; session{22,2} = 4:6; 
end