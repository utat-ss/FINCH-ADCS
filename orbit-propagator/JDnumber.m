function[jd] = JDnumber(year,month,day,hour,minute,second)
% Calculates the JD number for a given year,month and day

j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) + fix(275*month/9) + day + 1721013.5; 

ut = hour + minute/60 + second/3600; % UTC time in hour format

jd = j0 + ut/24; % JD number calculation
end