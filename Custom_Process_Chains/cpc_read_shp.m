function S = cpc_read_shp(filepath)
% CPC_READ_SHP  Minimal pure-MATLAB shapefile reader (polygon types only).
%
%   S = CPC_READ_SHP(filepath) returns a struct array with fields
%     .X  row-vector of x coordinates (NaN-separated rings)
%     .Y  row-vector of y coordinates (NaN-separated rings)
%
% Works without the Mapping Toolbox by parsing the binary .shp file
% directly. Supports polygon shape types 5 / 15 (PolygonZ) / 25 (PolygonM).
% Z and M values are ignored. Designed for the McCourt (2013) tectonic
% layers used by D_08; not a general shapefile library.

if ~exist(filepath,'file')
    error('cpc_read_shp: file not found: %s', filepath);
end

fid = fopen(filepath, 'rb');
if fid < 0
    error('cpc_read_shp: cannot open %s', filepath);
end
cleanupObj = onCleanup(@() fclose(fid));

% ---- file header ------------------------------------------------------
magic = fread(fid, 1, 'int32', 0, 'ieee-be');
if magic ~= 9994
    error('cpc_read_shp: bad magic %d in %s', magic, filepath);
end
fseek(fid, 24, 'bof');
file_len_words = fread(fid, 1, 'int32', 0, 'ieee-be');
file_len_bytes = file_len_words * 2;
fseek(fid, 32, 'bof');
shape_type = fread(fid, 1, 'int32', 0, 'ieee-le');
if ~ismember(shape_type, [5 15 25])
    error('cpc_read_shp: only polygon (5/15/25) supported, got %d', shape_type);
end

% ---- iterate records --------------------------------------------------
S = struct('X', {}, 'Y', {});
fseek(fid, 100, 'bof');
while ftell(fid) < file_len_bytes
    fread(fid, 1, 'int32', 0, 'ieee-be');                        % rec number
    rec_len_words = fread(fid, 1, 'int32', 0, 'ieee-be');
    rec_start     = ftell(fid);
    rec_len_bytes = rec_len_words * 2;

    rec_type = fread(fid, 1, 'int32', 0, 'ieee-le');
    if rec_type == 0
        S(end+1).X = []; S(end).Y = [];                          %#ok<AGROW>
        fseek(fid, rec_start + rec_len_bytes, 'bof');
        continue;
    end

    fread(fid, 4, 'double', 0, 'ieee-le');                       % bbox
    numParts  = fread(fid, 1, 'int32', 0, 'ieee-le');
    numPoints = fread(fid, 1, 'int32', 0, 'ieee-le');
    parts     = fread(fid, numParts, 'int32', 0, 'ieee-le');
    pts       = fread(fid, [2 numPoints], 'double', 0, 'ieee-le')';

    % --- build NaN-separated X/Y rings ---
    if numParts == 1
        X = pts(:,1).'; Y = pts(:,2).';
    else
        X = []; Y = [];
        parts_end = [parts; numPoints];   % treat as 1-based ends
        for p = 1:numParts
            i1 = parts_end(p) + 1;        % parts is 0-based offset
            i2 = parts_end(p+1);
            X = [X, pts(i1:i2,1).', NaN];     %#ok<AGROW>
            Y = [Y, pts(i1:i2,2).', NaN];     %#ok<AGROW>
        end
        X(end) = []; Y(end) = [];        % strip trailing NaN
    end
    S(end+1).X = X;                      %#ok<AGROW>
    S(end).Y = Y;

    fseek(fid, rec_start + rec_len_bytes, 'bof');
end
end
