%% ALINEA ramp metering controller

function u_alinea = alinea (u_alinea,ka,x_des,x)

    u_alinea = u_alinea + ka * (x_des - x);

end