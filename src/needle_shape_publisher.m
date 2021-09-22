
%% needle_shape_publisher.m
%
% script to publish the needle shape (testing script)
%
% - written by: Dimitri Lezcano

clear; 
%% Setup
% setup namespace
namespace = '/needle';
needlepose_pub_name = strcat(namespace, '/state/pose');
entrypoint_pub_name = strcat(namespace, '/state/skin_entry');

% setup needle shape publisher
needle_mechfile  = fullfile("../../shape-sensing/shapesensing_needle_properties.mat");
needle_paramfile = fullfile("../../data","needle_3CH_3AA",...
                            "needle_params-Jig_Calibration_11-15-20_weighted_weights.json");

%% Configure ROS Nodes
% setup needle pose publisher
helper.node = ros2node('TestNeedleHelper');
helper.pose_pub    = ros2publisher(helper.node, needlepose_pub_name, 'geometry_msgs/PoseStamped');
helper.entrypt_pub = ros2publisher(helper.node, entrypoint_pub_name, 'geometry_msgs/PointStamped');
helper.pose_sub    = ros2subscriber(helper.node, needlepose_pub_name, @needlepose_cb);

needleshape_talker = NeedleShapePublisher.fromFBGneedlefiles(needle_paramfile, needle_mechfile);

helper.shape_sub   = ros2subscriber(helper.node, '/needle/state/shape', @needleshape_cb);


%% Run through the loop
disp("Press [ENTER] to start publishing.")
pause;
L = 65; entrypt = [1;0;0];
fig = figure(1);
needleshape_talker.current_L = L;
while true
    pose_msg = generate_pose_msg(L, 0);
    skinentry_msg = generate_skinentry_msg(entrypt);
    helper.pose_pub.send(pose_msg);
    helper.entrypt_pub.send(skinentry_msg);
%     L = L + 5;
    needleshape_talker.publish_shape_cb();
    pause(0.1);
end

%% Helper functions
function msg = generate_pose_msg(L, theta_rot)
    msg = ros2message('geometry_msgs/PoseStamped');
    
    msg.pose.position.z = L;
    
    q = rotm2quat(rotz(theta_rot));
    msg.pose.orientation.w = q(1);
    msg.pose.orientation.z = q(4);
    
end

function needlepose_cb(msg)
%     fprintf("PoseSub: I heard z: %.2f, ang_z: %.2f\n", msg.position.z, msg.orientation.z);
end

function needleshape_cb(msg)
   % grab the positions
   [pmat, ~] = NeedleShapePublisher.posearray2shape(msg);
   
   figure(1);
   plot3(pmat(3,:), pmat(1,:), pmat(2,:));
   axis equal;

end

function msg = generate_skinentry_msg(pt)
    msg = ros2message('geometry_msgs/PointStamped');
    msg.point.x = pt(1);
    msg.point.y = pt(2);
    msg.point.z = pt(3);

end