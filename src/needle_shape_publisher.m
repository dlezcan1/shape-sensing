%% needle_shape_publisher.m
%
% script to publish the needle shape (testing script)
%
% - written by: Dimitri Lezcano

clear; 
%% Setup
% setup namespace
namespace = '/needle';
needlepose_pub_name = strcat(namespace, '/pose');

% setup needle shape publisher
needle_mechfile  = fullfile("../../shape-sensing/shapesensing_needle_properties.mat");
needle_paramfile = fullfile("../../data","needle_3CH_3AA",...
                            "needle_params-Jig_Calibration_11-15-20_weighted_weights.json");

%% Configure ROS Nodes
% setup needle pose publisher
needlepose_talker.node = ros2node('TestNeedlePoseTalker');
needlepose_talker.pub  = ros2publisher(needlepose_talker.node, needlepose_pub_name, 'geometry_msgs/Pose');
needlepose_talker.sub  = ros2subscriber(needlepose_talker.node, needlepose_pub_name, @needlepose_cb);

needleshape_talker = NeedleShapePublisher.fromFBGneedlefiles(needle_paramfile, needle_mechfile);

needlepose_talker.shape_sub = ros2subscriber(needlepose_talker.node, '/needle/shape/current', @needleshape_cb);



%% Run through the loop
disp("Press [ENTER] to start publishing.")
pause;
L = 65;
fig = figure(1);
needleshape_talker.current_L = L;
while true
    msg = generate_pose_msg(L, 0);
    needlepose_talker.pub.send(msg);
    L = L + 5;
    needleshape_talker.publish_shape_cb();
    pause(0.1);
end

%% Helper functions
function msg = generate_pose_msg(L, theta_rot)
    msg = ros2message('geometry_msgs/Pose');
    
    msg.position.z = L;
    
    q = rotm2quat(rotz(theta_rot));
    msg.orientation.w = q(1);
    msg.orientation.z = q(4);
    
end

function needlepose_cb(msg)
    fprintf("PoseSub: I heard z: %.2f, ang_z: %.2f\n", msg.position.z, msg.orientation.z);
end

function needleshape_cb(msg)
   % grab the positions
   [pmat, Rmat] = NeedleShapePublisher.posearray2shape(msg);
   
   figure(1);
   plot3(pmat(3,:), pmat(1,:), pmat(2,:));

end