function R = Rot_y(phi)
% rotation matrix about y axis
% - written by Jin Seob (Jesse) Kim

R = [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)];