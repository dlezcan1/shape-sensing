function R = Rot_z(theta)
% rotation matrix about z axis
% - written by Jin Seob (Jesse) Kim

R = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];