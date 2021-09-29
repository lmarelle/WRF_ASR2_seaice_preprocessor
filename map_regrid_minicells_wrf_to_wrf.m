function [mapping_regrid] = map_regrid_minicells_wrf_to_wrf(input_proj, output_proj, ncells)
%
% 2020/02/06, Louis Marelle
%
% -------- Create a "mapping" linking each grid box of an input --------
% -------- WRF grid to boxes of a different output WRF grid for regridding --------

% Program flow:
% - Check inputs
% - Create the empty mapping array mapping_regrid (dimensions of the output
%   grid)
% - Loop over the (i_input, j_input) elements of the input grid
%     - Break each element of the input grid into ncells * ncells minicells
%     - Loop over the ncells * ncells minicells in (i_input, j_input)
%         - Find the output grid cell containing the
%           minicell center, and add the coordinates (i_input, j_input) to
%           the element corresponding to this grid cell in mapping_regrid

% Variables:
%
% - Inputs
% input_proj: WRF input grid projection properties (matlab struct)
% output_proj: WRF output grid projection properties (matlab struct)
% ncells : number of minicells per cell in any direction (ncells in x, ncells in y,
%     ncells * ncells minicells per cell)

% - Outputs
% mapping_regrid: cell array containing lat/lon to WRF output grid mapping, with the 
%     same dimensions than the destination grid

% - Local
% i_cells, j_cells: WRF (i,j) coordinates 
%     of the center of the minicells of input grid box (i_input, j_input)
% icell: Minicell loop index

%-------- Check inputs --------
% Output grid coordinates should be floats ...
% Ncells should be a > 0 integer
if(floor(ncells) ~= ncells | ncells <= 0)
  error('Error, ncells should be a > 0 integer')
end


%-------- Compute some input grid and minicell properties
if(ncells == 1)
  minicells_0_center = 0;
else
  minicells_0 = linspace(-1/2, +1/2, ncells * 2 + 1);
  minicells_0_center = double(minicells_0(2:2:end));
end


%-------- Create the mapping --------
mapping_regrid = cell(output_proj.imax, output_proj.jmax);

for i_input = 1:input_proj.imax
  for j_input = 1:input_proj.jmax
    % For WRF input grid box (i_input, j_input), compute minicell center coordinates in WRF input i,j
    i_cells = repmat(double(i_input) + minicells_0_center, 1, ncells);
    j_cells = reshape(repmat(double(j_input) + minicells_0_center, ncells, 1), 1, ncells * ncells);

    % Convert WRF input minicell coordinates to lat,lon, and then to WRF output i,j
    [lat_cells, lon_cells] = ijll(i_cells, j_cells, input_proj.truelat1, input_proj.truelat2,...
                                  input_proj.hemi, input_proj.stdlon, input_proj.ref_lat,...
                                  input_proj.ref_lon, input_proj.ref_x, input_proj.ref_y,...
                                  input_proj.dx, input_proj.map_proj);
    clear i_cells j_cells
    [i_cells, j_cells] = llij(lat_cells, lon_cells, output_proj.truelat1,...
                              output_proj.truelat2, output_proj.hemi, output_proj.stdlon,... 
                              output_proj.ref_lat, output_proj.ref_lon, output_proj.ref_x,...
                              output_proj.ref_y, output_proj.dx, output_proj.map_proj);

    % Find the nearest WRF grid cell to each minicell
    i_cells = round(i_cells);
    j_cells = round(j_cells);

    % See if the minicell is in (true) or out (false) of the output grid
    valid_cells = j_cells >= 1 & j_cells <= output_proj.jmax & ...
                  i_cells >= 1 & i_cells <= output_proj.imax;

    % Update the mapping
    for icell = 1:ncells*ncells
      if(valid_cells(icell)) % Test that the minicell is within the bounds of the output grid
        % Append this minicell to the corresponding output mapping element
        mapping_regrid{i_cells(icell), j_cells(icell)} = [mapping_regrid{i_cells(icell), j_cells(icell)}, i_input j_input];
      end
    end
    
  end %j_input
end %i_input

end


