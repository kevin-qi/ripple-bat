function intersections = findLineBoxIntersections(m, c, box)
    % m: Slope of the line
    % c: y-intercept of the line
    % box: Rectangle boundaries [x_min, y_min, x_max, y_max]

    % Unpack the box boundaries
    x_min = box(1);
    y_min = box(2);
    x_max = box(3);
    y_max = box(4);
    
    % Initialize the list of intersections
    intersections = [];
    
    % Check intersections with vertical boundaries (left and right edges)
    for x = [x_min, x_max]
        y = m * x + c;
        if y >= y_min && y <= y_max
            intersections = [intersections; x, y];
        end
    end
    
    % Check intersections with horizontal boundaries (bottom and top edges)
    for y = [y_min, y_max]
        x = (y - c) / m;
        if x >= x_min && x <= x_max
            intersections = [intersections; x, y];
        end
    end
end