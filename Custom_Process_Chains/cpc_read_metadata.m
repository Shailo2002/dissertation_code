function meta = cpc_read_metadata()
% CPC_READ_METADATA  Load Custom_Process_Chains/station_metadata.csv if it
% exists. Returns a table with columns station, lat, lon, tectonic_group.
% Returns [] (empty) if the CSV is missing or has no data rows.

here = fileparts(mfilename('fullpath'));
csv  = fullfile(here, 'station_metadata.csv');

meta = [];
if ~exist(csv, 'file')
    fprintf('  [cpc_read_metadata] %s not found -- skipping.\n', csv);
    fprintf('  Copy station_metadata_TEMPLATE.csv to station_metadata.csv ');
    fprintf('and fill in your stations.\n');
    return;
end

% Strip comment lines (starting with #) before reading as table.
fid = fopen(csv,'r');
raw = textscan(fid, '%s', 'Delimiter','\n', 'Whitespace','');
fclose(fid);
lines = raw{1};
keep = true(numel(lines),1);
for i = 1:numel(lines)
    s = strtrim(lines{i});
    if isempty(s) || s(1) == '#'; keep(i) = false; end
end
lines = lines(keep);
if isempty(lines) || numel(lines) < 2
    fprintf('  [cpc_read_metadata] CSV has no data rows -- skipping.\n');
    return;
end

% write filtered content to temporary file and let readtable parse it.
tmp = [tempname '.csv'];
fid = fopen(tmp,'w');
for i = 1:numel(lines); fprintf(fid, '%s\n', lines{i}); end
fclose(fid);
try
    meta = readtable(tmp, 'Delimiter',',', 'TextType','string');
catch ME
    fprintf('  [cpc_read_metadata] failed to parse CSV: %s\n', ME.message);
    meta = [];
end
delete(tmp);
if ~isempty(meta) && height(meta) == 0
    meta = [];
end
end
