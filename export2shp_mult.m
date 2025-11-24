% define name and step


% After running "t = [start: interval : final]; findstep;" 
% Preallocate the struct array for efficiency
S = repmat(struct('Geometry', 'Line', 'X', [], 'Y', [], 'year', []), length(step), 1);

for i = 1:length(step)
    % Access the correct TransientSolution step
    current_step = step(i);
    
    % Extract glacier front line positions for the current step
    [glx, gly] = gl_position(md, current_step, 0);
    
    % Populate the struct with line data and the corresponding year
    S(i).Geometry = 'Line';
    S(i).X = glx;
    S(i).Y = gly;
    S(i).year = md.results.TransientSolution(current_step).time;
end

% Write the shapefile
shapewrite(S, name);

disp('Shapefile with multiple lines and year attribute has been written successfully');

