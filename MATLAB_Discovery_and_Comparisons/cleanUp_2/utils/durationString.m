% takes the output of a tic-toc
function [dispStr] = durationString(seconds)
%DURATIONSTRING Summary of this function goes here
%   Detailed explanation goes here
    days = 0; hours = 0; minutes = 0; % defaults
    if seconds >= 60
        minutes = floor(seconds / 60);    
        seconds = rem(seconds, 60);
        
        if minutes >= 60
            hours = floor(minutes / 60);    
            minutes = rem(minutes, 60);
            
            if hours >= 24
                days = floor(hours / 24);
                hours = rem(hours, 24);
            end
        end
    end
    dispStr = sprintf("%d d, %d hr, %d min and %0.3f sec", ...
                days, hours, minutes, seconds);
end

