function [e] = eccentricity(apogee, perigee)
    % calculating eccentricity from apogee and perigee
    
    Kr_e = 6378.1363; % earth mean radius
    
    ra = Kr_e + apogee;
    rp = Kr_e + perigee;
    
    e = (ra - rp)/(ra+rp);
    
end