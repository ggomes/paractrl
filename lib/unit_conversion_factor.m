function [cnv] = unit_conversion_factor(fromunits,tounits)
% units are {time,link length,flw,dty,spd} 
%   si = {sec,meter,veh/sec,veh/meter,meter/sec}
%   us = {hr,mile,veh/hr,veh/mile,mile/hr}
%   metric = {hr,km,veh/hr,veh/km,km/hr}

fromunits = lower(fromunits);
tounits = lower(tounits);

if(~ismember(tounits,{'si','us','metric'}))
    error('bad units given to unit_conversion_factor')
end

cnv = struct('time',1,'length',1,'flw',1,'dty',1,'spd',1);

if(strcmp(fromunits,tounits))
    return
end

meters_per_mile = 1609.34;
% feet_per_mile = 5280;
% meters_per_foot = meters_per_mile/feet_per_mile;
meters_per_km = 1000;
seconds_per_hour = 3600;
km_per_mile = meters_per_mile/meters_per_km;

us2si.time = seconds_per_hour;
us2si.length = meters_per_mile;
us2si.flw = 1/seconds_per_hour;
us2si.dty = 1/meters_per_mile;
us2si.spd = meters_per_mile/seconds_per_hour;
                
metric2si.time = seconds_per_hour;
metric2si.length = meters_per_km;
metric2si.flw = 1/seconds_per_hour;
metric2si.dty = 1/meters_per_km;
metric2si.spd = meters_per_km/seconds_per_hour;

us2metric.time = 1;
us2metric.length = meters_per_mile / meters_per_km;
us2metric.flw = 1;
us2metric.dty =  1/km_per_mile;
us2metric.spd = km_per_mile;
                
switch fromunits
    case 'si'
        switch tounits
            case 'us'
                cnv= inverse(us2si);
            case 'metric'
                cnv= inverse(metric2si);
        end
    case 'us'
        switch tounits
            case 'si'
                cnv = us2si;
            case 'metric'
                cnv = us2metric;
        end
    case 'metric'
        switch tounits
            case 'si'
                cnv = metric2si;
            case 'us'
                cnv= inverse(us2metric);
        end
end

function [B] = inverse(A)
B.time = 1/A.time;
B.length = 1/A.length;
B.flw = 1/A.flw;
B.dty = 1/A.dty;
B.spd = 1/A.spd;

