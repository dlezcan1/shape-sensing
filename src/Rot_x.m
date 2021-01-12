function R = Rot_x(phi)
% rotation matrix about x axis
% - written by Jin Seob (Jesse) Kim

R = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];