function r_threshold_cross=find_first_cross(f,th)
    globalmin=min(f);
    
    threshold=(max(f(1:5))-globalmin)*th+globalmin;
    r_threshold_cross=find(f<threshold,1,'first')-1;

end
    