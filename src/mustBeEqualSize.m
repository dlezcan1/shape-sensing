function mustBeEqualSize(a, b, axes)
    arguments
        a;
        b;
        axes = [];
    end
    % Test for equal size
    if isempty(axes)
        if ~isequal(size(a),size(b))
            eid = 'Size:notEqual';
            msg = 'Inputs must have equal size.';
            throwAsCaller(MException(eid,msg))
        end
    elseif length(axes) == 1
        if ~isequal(size(a, axes),size(b, axes))
            eid = 'Size:notEqual';
            msg = "Input axes must have equal size along axis: " + axes(1);
            throwAsCaller(MException(eid,msg))
        end
    else
        if ~isequal(size(a, axes(1)),size(b, axes(2)))
            eid = 'Size:notEqual';
            msg = sprintf("Input axes a[%d] and b[%d] must have equal size.", axes(1:2));
            throwAsCaller(MException(eid,msg))
        end
    end
end