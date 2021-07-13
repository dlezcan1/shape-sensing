%% needle_shape_publisher.m
%
% script to publish the needle shape (testing script)
%
% - written by: Dimitri Lezcano

clear; 
%% Setup
% setup topic names
namespace = '/needle';
node_name = strcat(namespace,'/needle_shape');
shapecurr_pub_name = strcat(namespace, '/shape/current');
shapepred_pub_name = strcat(namespace, '/shape/predicted');
sensor_sub_name = strcat(namespace, '/sensor/processed');

% configure node and topics
global needleshape_pub;
needleshape_pub.node = ros2node(node_name);
needleshape_pub.curr_pub = ros2publisher(needleshape_pub.node, shapecurr_pub_name,...
                                         'geometry_msgs/PoseArray');
needleshape_pub.pred_pub = ros2publisher(needleshape_pub.node, shapepred_pub_name,...
                                         'geometry_msgs/PoseArray');
needleshape_pub.sensor_sub = ros2subscriber(needleshape_pub.node, sensor_sub_name);


publish_needle_shape()

%% Helper functions
function publish_needle_shape()
    global needleshape_pub
    disp(needleshape_pub)

end
