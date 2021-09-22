%% NeedleShapePublisher.m
%
% class to publish the needle shape in ROS 2
%
% TODO
%   - join the needle parameters and the sensing stuff together into 1
%       struct
%   - subscriber to needle insertion length
%   - subscriber to theta0
%   - needle shape with rotation
%
% - written by: Dimitri Lezcano

classdef NeedleShapePublisher < handle
    properties
        % ROS 2 stuff
        node                    ros2node;
        current_pub             ros2publisher;
        predicted_pub           ros2publisher;
        sensor_sub              ros2subscriber;
        needlepose_sub          ros2subscriber;
        entrypoint_sub          ros2subscriber;
        sub_timeout             double {mustBePositive} = 20;
        % FBG Parameters
        num_samples             {mustBeInteger, mustBePositive} = 200; % number of samples to be gathered
        needleLength            double {mustBePositive}; % length of the entire needle
        num_channels            {mustBeInteger, mustBePositive};
        num_activeAreas         {mustBeInteger, mustBePositive};
        sensorLocations         (1,:) {mustBePositive}; % from the tip
        sensorCalMatrices       (2,:,:);
        aaReliabilityWeights    (1,:);
        % shape sensing parameters
        needleMechParams        struct;
        current_L               double = 0;
        kc_i                    double = 0.0025;
        w_init_i                (3,1)  = [0.0025; 0; 0];
        current_entry_point            = [];
    end
    methods
        % constructor
        function obj = NeedleShapePublisher(num_chs, num_aas, slocs, calMats, ...
                            needle_length, needle_mech_params, options)
            arguments
                num_chs {mustBeInteger, mustBePositive};
                num_aas {mustBeInteger, mustBePositive};
                slocs (1,:); % from the tip
                calMats (2,:,:);
                needle_length double {mustBePositive};
                needle_mech_params struct; % mechanical properties of th needle
                % ros options
                options.ns string = '/needle';
                options.node_name string = '/NeedleShape';
                options.num_samples {mustBeInteger, mustBePositive} = 200;
                options.timeout {mustBePositive} = 20;
                % shape sensing options
                options.kc_i double = 0.0025;
                options.w_init_i (3,1) = [0.0025; 0; 0];
            end
            % argument checking
            assert(numel(slocs) == num_aas);
            assert(size(calMats,1) == 2); % kx, ky
            assert(size(calMats,2) == num_chs); % CH1, CH2, ...
            assert(size(calMats,3) == num_aas); % AA1, AA2, ...
            
            % configure needle parameters
            obj.needleMechParams = needle_mech_params;
            obj.needleLength = needle_length;
            obj.num_channels = num_chs;
            obj.num_activeAreas = num_aas;
            obj.sensorLocations = slocs;
            obj.num_samples = options.num_samples;
            obj.sensorCalMatrices = calMats;
            
            % initial kc and w_init
            obj.kc_i = options.kc_i;
            obj.w_init_i = options.w_init_i;            
            
            % configure node, publishers and subscribers
            node_name = strcat(options.ns, options.node_name);
            curr_pub_name = strcat(options.ns, '/state/shape');
            pred_pub_name = strcat(options.ns, '/state/predicted_shape');
            sensor_sub_name = strcat(options.ns, '/sensor/processed');
            needlepose_sub_name = strcat(options.ns, '/state/pose');
            entrypoint_sub_name = strcat(options.ns, '/state/skin_entry');
            
            % configure node, publishers and subscibers
            obj.node = ros2node(node_name);
            obj.current_pub    = ros2publisher(obj.node, curr_pub_name,...
                                                'geometry_msgs/PoseArray');
            obj.predicted_pub  = ros2publisher(obj.node, pred_pub_name,...
                                                'geometry_msgs/PoseArray');
                                            
            obj.sensor_sub     = ros2subscriber(obj.node, sensor_sub_name,...
                                                'std_msgs/Float64MultiArray');
            obj.needlepose_sub = ros2subscriber(obj.node, needlepose_sub_name, @obj.needlepose_cb,...
                                                'geometry_msgs/PoseStamped'); 
            obj.entrypoint_sub = ros2subscriber(obj.node, entrypoint_sub_name, @obj.skin_entry_cb,...
                                                'geometry_msgs/PointStamped');
            obj.sub_timeout = options.timeout;
        end
        % callback to update needle shape parateters
        function needlepose_cb(obj, pose_msg)
           obj.current_L = pose_msg.pose.position.z; 
           
           current_orientation = pose_msg.pose.orientation.z;
           
        end
        % callback to publish the needle shape
        function obj = publish_shape_cb(obj)
            % check for current length viable
            if all(obj.current_L <= obj.sensorLocations)
                return;
            end
            
            % see if there is a current entry point
            if isempty(obj.current_entry_point)
                return;
            end
            
            % gather FBG samples and average
            fbg_peaks = zeros(obj.num_channels, obj.num_activeAreas);
            for i = 1:obj.num_samples
                fbg_msg = receive(obj.sensor_sub);
                fbg_peaks = fbg_peaks + obj.parseFBGPeaks(fbg_msg);

            end
            fbg_peaks = fbg_peaks/obj.num_samples;
            
%             % tensor way includes timeout reading
%             fbg_peaks = zeros(obj.num_channels, obj.num_activeAreas, 0);  
%             for i = 1:obj.num_samples
%                 fbg_msg = receive(obj.sensor_sub, obj.sub_timeout);
%                 fbg_peaks(:,:,i) = obj.parseFBGPeaks(fbg_msg); % tensor
%             end
%             fbg_peaks = mean(fbg_peaks,3); % tensor            

            % calibrate sensors and determine needle shape
            theta0 = 0; % TODO: change to current pose
            curvatures = calibrate_fbgsensors(fbg_peaks, obj.sensorCalMatrices);
            [pmat, wv, Rmat, kc, w_init] = singlebend_needleshape(curvatures, ...
                                            obj.sensorLocations, obj.current_L,...
                                            obj.kc_i, obj.w_init_i, theta0);
                                        
            % parse the position and Rmat into a geometry_msgs/PoseArray msg
            poseArr_msg = NeedleShapePublisher.shape2posearray(pmat, Rmat);
            
            % publish the message
            obj.current_pub.send(poseArr_msg);
            
            % updates
            obj.kc_i = kc;
            obj.w_init_i = w_init;
        end
        % parse FBG peak message
        function [peaks, peaks_struct] = parseFBGPeaks(obj, fbg_msg)
           % fbg_msg is of type std_msgs/Float64MultiArray 
            peaks = zeros(obj.num_channels, obj.num_activeAreas);
            peaks_struct = struct();
            idx = 1;
            for i = 1:numel(fbg_msg.layout.dim)
               ch_i = fbg_msg.layout.dim(i).label;
               size_i = fbg_msg.layout.dim(i).size/fbg_msg.layout.dim(i).stride;
               if size_i > 0
                    peaks(i,:) = fbg_msg.data(idx:idx + size_i - 1);
               
                    peaks_struct.(ch_i) = fbg_msg.data(idx:idx + size_i - 1);
               end
               idx = idx + size_i;
            end
        end
        % callback to update skin entry point
        function obj = skin_entry_cb(obj, msg)
            obj.current_entry_point = [msg.point.x; msg.point.y; msg.point.z];
            
        end
    end
    methods(Static)
       % turn needle shape pose into PoseArray
       function poseArray = shape2posearray(pmat, Rmat)
          arguments
              pmat (3,:);
              Rmat (3,3,:);
          end
          
          assert(size(pmat,2) == size(Rmat,3)); % they must have the same arclength points
          poseArray = ros2message('geometry_msgs/PoseArray');
          pose_i = ros2message('geometry_msgs/Pose');
          for i = 1:size(pmat,2)
              p_i = pmat(:,i);
              q_i = rotm2quat(Rmat(:,:,i));
              
              pose_i.position.x = p_i(1);
              pose_i.position.y = p_i(2);
              pose_i.position.z = p_i(3);
              
              pose_i.orientation.w = q_i(1);
              pose_i.orientation.x = q_i(2);
              pose_i.orientation.y = q_i(3);
              pose_i.orientation.z = q_i(4);
              
              poseArray.poses(i) = pose_i;
          end
       end
       function [pmat, Rmat] = posearray2shape(msg)
           N = numel(msg.poses);
           pmat = zeros(3,N);
           Rmat = zeros(3,3,N);
           
           for i = 1:N
               pmat(:,i) = [msg.poses(i).position.x; 
                            msg.poses(i).position.y; 
                            msg.poses(i).position.z];
           
               quat = [msg.poses(i).orientation.w;
                       msg.poses(i).orientation.x;
                       msg.poses(i).orientation.y;
                       msg.poses(i).orientation.z]';
                   
               Rmat(:,:,i) = quat2rotm(quat);
           end
       end
       % function to parse FBG json into a new object
       function obj = fromFBGneedlefiles(json_filename, needle_mechanics_file, options)
           arguments
               json_filename         string;
               needle_mechanics_file string;
               options.ns string = '/needle';
               options.node_name string = '/NeedleShape';
               options.num_samples {mustBeInteger, mustBePositive} = 200;
               options.timeout {mustBePositive} = 20;
               % shape sensing options
               options.kc_i double = 0.0025;
               options.w_init_i (3,1) = [0.0025; 0; 0];
           end
           % read in the data
           needle_params    = parse_fbgneedle_json(json_filename);
           needle_mechanics = load(needle_mechanics_file);
           
           % check for fields
           assert(isfield(needle_mechanics, 'B'));
           if ~isfield(needle_mechanics, 'Binv')
               needle_mechanics.Binv = inv(needle_mechanics.B);
           end
           obj = NeedleShapePublisher(needle_params.num_channels,...
                                      needle_params.num_activeAreas,...
                                      needle_params.slocs_tip,...
                                      needle_params.aa_calMats,...
                                      needle_params.length,...
                                      needle_mechanics, ...
                                      'ns', options.ns,...
                                      'node_name', options.node_name,...
                                      'num_samples', options.num_samples,...
                                      'timeout', options.timeout,...
                                      'kc_i', options.kc_i,...
                                      'w_init_i', options.w_init_i);
           obj.aaReliabilityWeights = needle_params.aa_weights;
       end
    end
    
end 