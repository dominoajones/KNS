function export2nc(md, dataset, filename, varname)
% export2nc - Saves data as netcdf file with multiple variables support
% Usage:
%   export2nc(md, dataset, filename, varname)
%
% md: model struct containing mesh info
% dataset: data array to write
% filename: name of netcdf file
% varname: string, variable name inside netcdf (e.g., 'Thickness')

if nargin < 4
    varname = 'variable';
end

% Define mesh/grid parameters
index = md.mesh.elements;
x = md.mesh.x;
y = md.mesh.y;
data = dataset;

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

grid_resolution = 100; % meters
xgrid = x_min:grid_resolution:x_max;
ygrid = y_min:grid_resolution:y_max;
default_value = NaN;

% Interpolate data from mesh to grid
grid = InterpFromMeshToGrid(index, x, y, data, xgrid, ygrid, default_value);

% Check if file exists
if exist(filename, 'file') ~= 2
    % Create new file and dimensions
    nccreate(filename, 'x', 'Dimensions', {'x', length(xgrid)});
    nccreate(filename, 'y', 'Dimensions', {'y', length(ygrid)});

    ncwrite(filename, 'x', xgrid);
    ncwriteatt(filename, 'x', 'units', 'meters');
    ncwriteatt(filename, 'x', 'standard_name', 'projection_x_coordinate');
    ncwriteatt(filename, 'x', 'long_name', 'x coordinate of projection');

    ncwrite(filename, 'y', ygrid);
    ncwriteatt(filename, 'y', 'units', 'meters');
    ncwriteatt(filename, 'y', 'standard_name', 'projection_y_coordinate');
    ncwriteatt(filename, 'y', 'long_name', 'y coordinate of projection');

    % Global attributes for CRS
    ncwriteatt(filename, '/', 'Conventions', 'CF-1.6');
    ncwriteatt(filename, '/', 'title', 'Model Data Export');
    ncwriteatt(filename, '/', 'institution', 'Your Institution');
    ncwriteatt(filename, '/', 'source', 'Model Data');
    ncwriteatt(filename, '/', 'references', 'Your References');

    % CRS variable
    nccreate(filename, 'crs', 'Datatype', 'char', 'Dimensions', {'string', 1});
    ncwriteatt(filename, 'crs', 'grid_mapping_name', 'polar_stereographic');
    ncwriteatt(filename, 'crs', 'longitude_of_projection_origin', 0.0);
    ncwriteatt(filename, 'crs', 'latitude_of_projection_origin', 90.0);
    ncwriteatt(filename, 'crs', 'standard_parallel', 70.0);
    ncwriteatt(filename, 'crs', 'straight_vertical_longitude_from_pole', -45.0);
    ncwriteatt(filename, 'crs', 'false_easting', 0.0);
    ncwriteatt(filename, 'crs', 'false_northing', 0.0);
    ncwriteatt(filename, 'crs', 'EPSG_code', 'EPSG:3413');
end

% Create variable with variable name
nccreate(filename, varname, 'Dimensions', {'x', length(xgrid), 'y', length(ygrid)});
ncwrite(filename, varname, grid');

% Add variable attributes (customize per variable)
ncwriteatt(filename, varname, 'long_name', varname);
ncwriteatt(filename, varname, 'units', 'unknown'); % You can customize this per variable
ncwriteatt(filename, varname, 'coordinates', 'x y');
ncwriteatt(filename, varname, 'grid_mapping', 'crs');

disp(['Variable "' varname '" written to ' filename]);
end
