function session = get_scanspersession_eye(exp_type)

switch exp_type
    
    case 'AWAKE_EYE'

        %-----------------------------------------------------------
              
        
        session{23,1}=[1 2 4 5 6]; %session 23 is originally 1 session but        
        % scan 5: eyemovie frozen, lots of movement

        session{25,1} = [(1:3) (5:7)];
        % scan4: lots of movement and eyemovie frozen, scan 6: eye closed

        
        session{26,1} = 1:5;

        
        session{27,1} = 2:5;        
        %originally two sessions acquired but the second session was discard due to bad eye

        session{28,1} = [(1:4) 5 7 8];

        session{29,1} = 1:6;

        session{30,1} = [1 2 (4:6)];
        %scan3: eyefrozen, 2nd session not good.

        session{32,1} = 1:6;

        session{33,1} = 1:7;
        %scan3: eyemovie frozen in the middle of scan
        % scans acquired: 1:8, but I execlude 8th scan due to small pupil
        % size and difficulty of eyetracking.

        session{35,1} = [(1:4) (6:9)]; 
        %6:lots of movement

        session{36,1} = [(1:3) (5:9)];

        session{37,1} = 1:8; 
        % scan 7:  not good eye,

        session{40,1} = 1:6; 
        % scan 6: pupil contracted a lot

end