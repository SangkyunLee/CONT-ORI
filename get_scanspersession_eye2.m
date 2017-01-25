function session = get_scanspersession_eye2(exp_type)
% data sets with two sessions timely separated
switch exp_type
    
    case 'AN'
                %---------------------------------        
        session{17,1} = 1:3; session{17,2} = 4:6; 
        session{18,1} = 1:3; session{18,2} = 4:6; 
        session{19,1} = 1:3; session{19,2} = 4:7; 
        session{20,1} = 1:3; session{20,2} = 4:6; 
        session{21,1} = 1:3; session{21,2} = 4:6; 
        session{22,1} = 1:3; session{22,2} = 4:6; 
    
    case 'AWAKE_EYE'

        %-----------------------------------------------------------
        % session 26, 27,28  look not very good response time course
        
        
        
        %session{23,1}=[2 4];  session{23,2}=[1 5 6]; %session 23 is originally 1 session but        
        % scan 5: eyemovie frozen, lots of movement

        session{25,1} = 1:3; session{25,2} = 5:7;
        % scan4: lots of movement and eyemovie frozen, scan 6: eye closed

        
        session{26,1} = 1:3;session{26,2} = 4:5;

        
        %session{27,1} = 2:3;session{27,2} = 4:5;     
        
        
        session{29,1} = [1 2 3]; session{29,2} = [4 5 6];
        %session{29,1} = [1 3]; session{29,1} = [2 4 5 6];

        session{30,1} = [1 2 4]; session{30,2} = 5:6;

        session{32,1} = 1:4; session{32,2} = 5:6;

        session{33,1} = 1:4; session{33,2} = 5:7;
        %scan3: eyemovie frozen in the middle of scan
        % scans acquired: 1:8, but I execlude 8th scan due to small pupil
        % size and difficulty of eyetracking.

        %session{35,1} = 1:4; session{35,2} = 6:9; 
        %6:lots of movement

        session{36,1} = 1:3; session{36,2} = 5:9;

       session{40,1} = 1:3; session{40,2} = 4:6; 
        % scan 6: pupil contracted a lot

end