function lengthErr = errorLength3D(ref1, ref2)
% take in two 3 by N matrices called ref1 and ref2
% calculate the distance between the matrices


if size(ref1,2) == 3
    ref1 = ref1'; % change matrix to 3 x N
elseif size(ref1,2) ~=3 && size(ref1,1) ~=3
    error('1st matrix must be 3 x N or N x 3');
end

if size(ref2,2) == 3
    ref2 = ref2';
elseif size(ref2,2) ~=3 && size(ref2,2) ~=3
    error('2nd matrix must be 3 x N or N x 3');
end
%%
% sort by z axis
ref1 = sortrows(ref1',3)';
ref2 = sortrows(ref2',3)';

if size(ref1,2) < size(ref2,2)
    shorter_ref = ref1;
    longer_ref = ref2;
else
    shorter_ref = ref2;
    longer_ref = ref1;
end
% distance calculation based on closest point
lengthErr = zeros(1,length(shorter_ref));
for i = 1:length(shorter_ref)
    templength = 900;
    for j = 1:length(longer_ref)
        temptemp = norm(longer_ref(:,j)-shorter_ref(:,i));
        if temptemp < templength
            % make new length, the shortest length bt two points
            templength = temptemp;
        else 
            break % exit the for loop
        end
    end
    lengthErr(i) = templength;
end

%% distance calculation based on z-coordinates
lengthErrz = zeros(1,length(shorter_ref));
for i = 1:length(shorter_ref)
    temp_z = shorter_ref(3,i);
    % find closest z-coordinate to the shorter array, on the longer array
    [~, closest_ind] = min(abs(longer_ref(3,:)-temp_z));
    lengthErrz(i) = norm(longer_ref(:,closest_ind)-shorter_ref(:,i));
end


end