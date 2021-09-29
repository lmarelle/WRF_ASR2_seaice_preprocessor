%-------- Read ASR2 ice depth obtained from NCAR-RDA and interpolate onto WRF sea ice --------
%
% Louis Marelle for LATMOS and IGE, 2020/02/06
%
clear
close all

%-------- Input --------
ASR2_DIR = 'path to ASR2 data';
WPS_DIR = 'path to WPS data'
USE_ASR2_SEAICE = true; % If set to false, default constant values are used instead of ASR2
MAX_DOMAINS = 1;  % number of WRF domains for nested runs


%-------- Parameters --------
SNOWDEPTH_SEAICE = 0.3; % snow depth in meters, from N-ICE campaign (before the melt season begins)
ICEDEPTH_SEAICE = 1.5; % seaice depth in meters, only used if USE_ASR2_SEAICE = false
ALBEDO_SEAICE = 0.82; % sea ice albedo, from N-ICE (before the melt season begins)
SNOWDENSITY_SEAICE = 200.0; % snow density in kg/m3


%-------- Initialize --------
for idomain = 1:MAX_DOMAINS
  domain = ['d0', num2str(idomain)];
  disp(domain)

  % Get met_em dates from the met_em_<domain>* filenames in WPS_DIR/
  path_files = dir([WPS_DIR, '/met_em.', domain, '*.nc']);
  met_em_filenames = { path_files.name };
  met_em_filenames = cell2mat(met_em_filenames');
  met_em_dates = met_em_filenames(:, end-21:end-3);
  met_em_dates = datenum(met_em_dates, 'yyyy-mm-dd_HH:MM:SS');
  years_list = year(met_em_dates)';
  months_list = month(met_em_dates)';
  days_list = day(met_em_dates)';
  hours_list = hour(met_em_dates)';

  asr2_dates = met_em_dates(1):1:met_em_dates(end);

  % Get WRF projection parameters
  wrf_proj = get_WRF_proj([WPS_DIR, '/', met_em_filenames(1, :)]);

  if (USE_ASR2_SEAICE)
    % Get ASR2 xlat/xlong
    asr2_filename = [ASR2_DIR, '/ICEDEPTH.asr15km.anl.2D.20120422.nc'];
    ncid = netcdf.open(asr2_filename, 'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid, 'XLAT');
    asr2_xlat = netcdf.getVar(ncid, varid);
    varid = netcdf.inqVarID(ncid, 'XLONG');
    asr2_xlong = netcdf.getVar(ncid, varid);
    netcdf.close(ncid);
    asr2_imax = size(asr2_xlat, 1);
    asr2_jmax = size(asr2_xlat, 2);

    % Set ASR2 projection
    %TODO Warning !! I could not find the actual projection used in ASR2 anywhere in
    % the ASR2 references. So this is not the actual projection, but a very similar one found
    % with trial and error - max errors in the projection are of the order of 1
    % km, do not use for very HR applications with resolutions of 3km or more -
    % get the actual projection from the ASR2 team instead
    asr2_proj.imax = asr2_imax;
    asr2_proj.jmax = asr2_jmax;
    asr2_proj.truelat1 = 60;
    asr2_proj.truelat2 = NaN;
    asr2_proj.hemi = 1;
    asr2_proj.stdlon = -175;
    asr2_proj.ref_lat = 90;
    asr2_proj.ref_lon = -175;
    asr2_proj.ref_x = asr2_imax/2+0.5;
    asr2_proj.ref_y = asr2_jmax/2+0.5;
    asr2_proj.dx = 15;
    asr2_proj.map_proj = 2;

    % Calculate minicell number and initialize regrid
    % To regrid to the same resolution we want at least 4 cells
    ncells = ceil(floor(4*asr2_proj.dx/wrf_proj.dx) / 2) * 2;   ncells(ncells < 2) = 1;
    tic
    mapping_regrid = map_regrid_minicells_wrf_to_wrf(asr2_proj, wrf_proj, ncells);
    toc
  end % if USE_ASR2_SEAICE

  % Loop on days
  ndates = length(met_em_dates);
  for idate = 1:ndates
    met_em_date = met_em_dates(idate);
    asr2_date = floor(met_em_date);
    disp(['  ', num2str(idate), ' ', datestr(met_em_date)])

    % Open WRF met_em SEAICE for the current date
    met_em_filename = [WPS_DIR, '/met_em.', domain, '.', datestr(met_em_date, 'yyyy-mm-dd_HH:MM:SS'), '.nc'];
    ncid = netcdf.open(met_em_filename, 'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid, 'SEAICE');
    wrf_seaice = netcdf.getVar(ncid, varid);
    varid = netcdf.inqVarID(ncid, 'SNOWH');
    wrf_snowh = netcdf.getVar(ncid, varid);
    varid = netcdf.inqVarID(ncid, 'SNOW');
    wrf_snow = netcdf.getVar(ncid, varid);
    netcdf.close(ncid);

    if (USE_ASR2_SEAICE)
      % Open ASR2 ICEDEPTH
      asr2_filename = [ASR2_DIR, '/ICEDEPTH.asr15km.anl.2D.', datestr(asr2_date, 'yyyymmdd'), '.nc'];
      ncid = netcdf.open(asr2_filename, 'NC_NOWRITE');
      varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
      asr2_icedepth = netcdf.getVar(ncid, varid);
      netcdf.close(ncid);

      % Regrid ASR2 field to WRF grid
      input_mask = asr2_icedepth;
      input_mask(input_mask > 0) = 1;
      input_mask(input_mask <= 0) = 0;
      [asr2_icedepth_wrfgrid] = regrid_area_minicells_mask(asr2_icedepth, mapping_regrid, input_mask);
      asr2_icedepth_wrfgrid(wrf_seaice <= 0.0) = 0.0;
      % Default icedepth is set at 2.0 m (Kwok and Rothrock 2009) - would be better
      % to set this based on a regional correlation between ice cover and ice depth
      asr2_icedepth_wrfgrid(wrf_seaice > 0.0 & isnan(asr2_icedepth_wrfgrid)) = 2.0;
      asr2_icedepth_wrfgrid(isnan(asr2_icedepth_wrfgrid)) = 0.0;
    else
      asr2_icedepth_wrfgrid(:) = ICEDEPTH_SEAICE;

    end % if USE_ASR2_SEAICE

    % Set SNOWH on seaice to the value of SNOWDEPTH_SEAICE
    wrf_snowh(wrf_seaice > 0.0) = SNOWDEPTH_SEAICE;
    wrf_snow(wrf_seaice > 0.0) = wrf_snowh(wrf_seaice > 0.0) * SNOWDENSITY_SEAICE;

    % Set ALBSI to the default value of ALBEDO_SEAICE
    wrf_albsi = zeros(size(wrf_seaice));
    wrf_albsi(wrf_seaice > 0.0) = ALBEDO_SEAICE;

    % Create and write ICEDEPTH & ALBSI to all met_em files and write
    % FLAG_ICEDEPTH,FLAG_ALBSI
    ncid = netcdf.open(met_em_filename, 'NC_WRITE');
    % If the ICEDEPTH variable already exists do not create it
    try
      varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
    catch
      netcdf.reDef(ncid)
      dimid_we = netcdf.inqDimID(ncid, 'west_east');
      dimid_sn = netcdf.inqDimID(ncid, 'south_north');
      dimid_time = netcdf.inqDimID(ncid, 'Time');
      varid = netcdf.defVar(ncid, 'ICEDEPTH', 'float', [dimid_we, dimid_sn, dimid_time]);
      varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
      netcdf.putAtt(ncid, varid, 'FieldType', int32(104))
      netcdf.putAtt(ncid, varid, 'MemoryOrder', 'XY')
      netcdf.putAtt(ncid, varid, 'description', 'ICE DEPTH')
      netcdf.putAtt(ncid, varid, 'units', 'm')
      netcdf.putAtt(ncid, varid, 'stagger', '')
      netcdf.putAtt(ncid, varid, 'coordinates', 'XLONG XLAT')
      netcdf.endDef(ncid)
    end
    % If the ALBSI variable exists do not create it
    try
      varid = netcdf.inqVarID(ncid, 'ALBSI');
    catch
      netcdf.reDef(ncid)
      dimid_we = netcdf.inqDimID(ncid, 'west_east');
      dimid_sn = netcdf.inqDimID(ncid, 'south_north');
      dimid_time = netcdf.inqDimID(ncid, 'Time');
      varid = netcdf.defVar(ncid, 'ALBSI', 'float', [dimid_we, dimid_sn, dimid_time]);
      varid = netcdf.inqVarID(ncid, 'ALBSI');
      netcdf.putAtt(ncid, varid, 'FieldType', int32(104))
      netcdf.putAtt(ncid, varid, 'MemoryOrder', 'XY')
      netcdf.putAtt(ncid, varid, 'description', 'SEA ICE ALBEDO')
      netcdf.putAtt(ncid, varid, 'units', ' ')
      netcdf.putAtt(ncid, varid, 'stagger', '')
      netcdf.putAtt(ncid, varid, 'coordinates', 'XLONG XLAT')
      netcdf.endDef(ncid)
    end
    netcdf.reDef(ncid)
    % Write variables and attributes to met_em files
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid, varid, 'FLAG_ICEDEPTH', int32(1));
    netcdf.putAtt(ncid, varid, 'FLAG_ALBSI', int32(1));
    netcdf.endDef(ncid)
    varid = netcdf.inqVarID(ncid, 'ICEDEPTH');
    netcdf.putVar(ncid, varid, asr2_icedepth_wrfgrid);
    varid = netcdf.inqVarID(ncid, 'ALBSI');
    netcdf.putVar(ncid, varid, wrf_albsi);
    varid = netcdf.inqVarID(ncid, 'SNOWH');
    netcdf.putVar(ncid, varid, wrf_snowh);
    varid = netcdf.inqVarID(ncid, 'SNOW');
    netcdf.putVar(ncid, varid, wrf_snow);
    netcdf.close(ncid);

  end % for idate = 1:ndates
end % for idomain

