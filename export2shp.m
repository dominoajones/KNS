S = struct('Geometry', 'Line', ...
           'X', {glx}, ...
           'Y', {gly});

% Specify the filename for the shapefile

% Write the shapefile
shapewrite(S, name);

disp('Shapefile has been written successfully');
