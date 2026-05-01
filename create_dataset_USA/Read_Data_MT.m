function[CData, MT] = Read_Data_MT(fname, CData, type)
% Read just Station or the Full Data, depending on the type of call
% First four column are for Z and next 2 are for Tipper, by Z I mean full
% impedance tensor or Phase Tensor, or off diag impedance, or app res and
% phase. even if some data is not there, fill the full four columns, they
% will remain empty, to be consistent with my approach for inversion. One
% must read the full data but want to invert for the selected component
% only. Also, full vertical components always comes in the end. !!!
% Missing data, either site, period or components are also taken care of.
% By default we are following exp{-iwt} sign convention, for exp{+iwt}
% sign convention impedance and vertical TF will be taken as complex
% conjugate while reading and while writing data also. Updated because of
% the Phaser Tensor.

if nargin < 3 || isempty(type)
    type = 'ReadData';
end
if ~strcmp(type,'ReadData')
    type = 'ReadStation';
end

MT.ncomp = 0;
CData.mu0 = 4*pi*1.0E-07;

% Quickly read the data file for finding all unique stations and periods. 
[codes, periods] = Sort_Station_Frequency(fname);
MT.nfreq = length(periods);    % Maximum number of periods
MT.nobs = size(codes,1);       % Maximum number of Stations
MT.T = periods;
% Initialising Data.MT, complete data set
MT.id = ones(MT.nobs,6,MT.nfreq);
MT.data = zeros(MT.nobs,6,MT.nfreq);
MT.err = ones(MT.nobs,6,MT.nfreq);

if exist(fname)
    fid = fopen(fname,'r');  % read file
else
    fprintf('%s %s \n','No file exist by the name : ',fname)
    error('Stop');
end
    
id = 0;
while(1)
    
    id = id + 1;
    % read header
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    tmp = char(tmp{1});
    if isempty(tmp); break; end
    MT.header = tmp;                                % header
    % read the block header
    blockinfo = textscan(fid,'# %s',1,'delimiter','\n');
    header = char(blockinfo{1});
    blockinfo = textscan(fid,'> %s',6,'delimiter','\n');
    blockinfo = char(blockinfo{1});
    dataType = blockinfo(1,:);
    signstr  = blockinfo(2,:);
    if findstr(signstr,'-')
        MT.isign = -1;
        MT.ii = complex(0,-1);
        MT.conjugate = 0;
    else  % complex data as complex conjugate
        MT.isign = -1;
        MT.ii = complex(0,-1);        
        MT.conjugate = 1;
    end
    MT.signstr = signstr;                            % field type
    % need copy only for impedance units, tipper is dimension less.
    typeUnits = strtrim(blockinfo(3,:));                % units
    MT.orientation = sscanf(blockinfo(4,:),'%f',1);  % orientation
    MT.origin = sscanf(blockinfo(5,:),'%f',3);       % origin
    tmp  = sscanf(blockinfo(6,:),'%d',2);               % nfreq and nobs, discarded now

    % All datatype and their corresponding units stored here, in cells
    MT.Nfreq{id} = tmp(1);
    MT.Nobs{id} = tmp(2);
    MT.datatype{id} = dataType;      
    MT.typeUnits{id} = typeUnits;
    MT.lineinfo{id} = header;

    switch strtrim(dataType)
        case 'Off_Diagonal_Impedance'
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            ncomp = 2;  MT.ncomp = MT.ncomp + ncomp;
            comp = ['ZXY';'ZYX'];   MT.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:size(codes,1)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            MT.cord.x(k,1) = data{5}(ind);
                            MT.cord.y(k,1) = data{6}(ind);
                            MT.cord.z(k,1) = data{7}(ind);
                            MT.Stat.code(k,:) = data{2}(ind);
                            MT.Stat.lat(k,1) = data{3}(ind);
                            MT.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                if MT.conjugate == 1   % Data is complex conjugate
                                    MT.data(k,i+1,j) = complex(data{9}(ind), -data{10}(ind));
                                else
                                   MT.data(k,i+1,j) = complex(data{9}(ind), data{10}(ind)); 
                                end
                                MT.err(k,i+1,j) = data{11}(ind);
                                if abs(MT.data(k,i+1,j)) == 0
                                    MT.id(k,i+1,j) = 0;
                                end
                            end
                        else
                            MT.id(k,i+1,j) = 0;
                        end
                    end
                end
            end
        case 'Full_Impedance'                    % Reading Impedance
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            ncomp = 4;  MT.ncomp = MT.ncomp + ncomp;
            comp = ['ZXX';'ZXY';'ZYX';'ZYY'];   MT.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:size(codes,1)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            % disp(data{2}(ind))
                            % disp(data{3}(ind))
                            % disp([i j k])
                            % disp(data{3}(ind))
                            MT.Stat.lat(k,:) = data{3}(ind);
                            MT.Stat.lon(k,:) = data{4}(ind);
                            MT.Stat.code(k,:) = data{2}(ind);
                            MT.cord.x(k,1) = data{5}(ind);
                            MT.cord.y(k,1) = data{6}(ind);
                            MT.cord.z(k,1) = data{7}(ind);
                            if strcmp(type,'ReadData')
                                if MT.conjugate == 1   % Data is complex conjugate
                                    MT.data(k,i,j) = complex(data{9}(ind), -data{10}(ind));
                                else
                                   MT.data(k,i,j) = complex(data{9}(ind), data{10}(ind)); 
                                end
                                MT.err(k,i,j) = data{11}(ind);
                                if abs(MT.data(k,i,j)) == 0
                                    MT.id(k,i,j) = 0;
                                end
                            end
                        else
                            MT.id(k,i,j) = 0;
                        end
                    end
                end
            end
        case 'Off_Diagonal_Rho_Phase'                    % Apparent Resistivity and phase, all data set is rael
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');
            ncomp = 4;  MT.ncomp = MT.ncomp + ncomp;
            comp = ['RHOXY';'PHSXY';'RHOYX';'PHSYX'];   MT.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:size(codes,1)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            MT.cord.x(k,1) = data{5}(ind);
                            MT.cord.y(k,1) = data{6}(ind);
                            MT.cord.z(k,1) = data{7}(ind);
                            MT.Stat.code(k,:) = data{2}(ind);
                            MT.Stat.lat(k,1) = data{3}(ind);
                            MT.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                MT.data(k,i,j) = data{9}(ind);
                                MT.err(k,i,j) = data{10}(ind);
                                if abs(MT.data(k,i,j)) == 0
                                    MT.id(k,i,j) = 0;
                                end
                            end
                        else
                            MT.id(k,i,j) = 0;
                        end
                    end
                end
            end
        case 'Phase_Tensor'                    % Phase Tensor, all data set is rael
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');
            ncomp = 4;  MT.ncomp = MT.ncomp + ncomp;
            comp = ['PTXX';'PTXY';'PTYX';'PTYY'];   MT.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:size(codes,1)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            MT.cord.x(k,1) = data{5}(ind);
                            MT.cord.y(k,1) = data{6}(ind);
                            MT.cord.z(k,1) = data{7}(ind);
                            MT.Stat.code(k,:) = data{2}(ind);
                            MT.Stat.lat(k,1) = data{3}(ind);
                            MT.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                MT.data(k,i,j) = data{9}(ind);
                                MT.err(k,i,j) = data{10}(ind);
                                if abs(MT.data(k,i,j)) == 0
                                    MT.id(k,i,j) = 0;
                                end
                            end
                        else
                            MT.id(k,i,j) = 0;
                        end
                    end
                end
            end
        case 'Full_Vertical_Components'          % Reading Tipper                           
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            ncomp = 2;  MT.ncomp = MT.ncomp + ncomp;
            comp = ['TX';'TY'];     MT.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:size(codes,1)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp	% ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            MT.cord.x(k,1) = data{5}(ind);
                            MT.cord.y(k,1) = data{6}(ind);
                            MT.cord.z(k,1) = data{7}(ind);
                            MT.Stat.code(k,:) = data{2}(ind);
                            MT.Stat.lat(k,1) = data{3}(ind);
                            MT.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                if MT.conjugate == 1   % Data is complex conjugate
                                    MT.data(k,i+4,j) = complex(data{9}(ind), -data{10}(ind));
                                else
                                   MT.data(k,i+4,j) = complex(data{9}(ind), data{10}(ind)); 
                                end
                                MT.err(k,i+4,j) = data{11}(ind);
                                if abs(MT.data(k,i+4,j)) == 0
                                    MT.id(k,i+4,j) = 0;
                                end
                            end
                        else
                            MT.id(k,i+4,j) = 0;
                        end
                    end
                end
            end
        otherwise
            disp('Unknown data type');
            break;
    end
end
MT.omega = 2*pi./MT.T;    % Omega's
MT.id = MT.id;          % Storing component in the MT

% Unit conversion introduced here. we will convert the data to [v/m]/[T]
% as we work in these units and later if required, re convert to the old
% units, unit conversion for data and error as well and only for impedance
% and not for dimension less quantity like VTF and Phase Tensor
MT.NoofDataTypes = length(MT.datatype);  % Number of data types

MT.f1 = -1/MT.ii;
MT.factor = CData.mu0;
MT.factor1 = 1;  MT.factor2 = 1;

for i = 1:MT.NoofDataTypes
    datatype = strtrim(char(MT.datatype{i}));
    switch(datatype)
        case {'Off_Diagonal_Impedance','Full_Impedance'}
            typeUnits = strtrim(char(MT.typeUnits{1}));
            switch(typeUnits)
                case '[V/m]/[T]'
                    MT.factor1 = 1;
                    MT.factor2 = 1;
                case '[mV/km]/[nT]'
                    MT.factor1 = 1000;
                    MT.factor2 = 1/1000;
                case {'[V/m]/[A/m]','Ohm'}
                    MT.factor1 = 1/CData.mu0;
                    MT.factor2 = CData.mu0;
            end
            % Convert data to desired units is required
            if (MT.factor1 ~= 1) && strcmp(type,'ReadData')
                MT.data(:,1:4,:) = MT.factor1 * MT.data(:,1:4,:);
                MT.err(:,1:4,:) = MT.factor1 * MT.err(:,1:4,:);
            end
        case 'Off_Diagonal_Rho_Phase'
            if isfield(CData,'rho')
                % Convert rho to log10(rho) if so desired
                if (strcmp(type,'ReadData') && strcmp(CData.rho,'log10'))
                    MT.data(:,1,:) = log10(MT.data(:,1,:));
                    MT.data(:,3,:) = log10(MT.data(:,3,:));
                    MT.err(:,1,:) = log10(MT.err(:,1,:));
                    MT.err(:,3,:) = log10(MT.err(:,3,:));
                end
            end
    end
end

% z positive downwards, for elevation z is negative,
% checking for topography or bathymetry, if exists;
% Find new location of stations on the basis of new reference, if its
% exists. My new reference is z=0; and we are considering Z positive
% downwards. The station with the maximum elevation is now at Z=0, but
% while writting nothing changes. However, this assumption is not valid
% when topography exceeds my maximum elevation stations hence we need to
% consider this fact also while dealing with them.
MT.refer = 0.0;  % Initialise reference elevation
t1 = min(MT.cord.z);
t2 = max(MT.cord.z);
if t1 == 0 && t2 == 0
    MT.elevation = 'NoElevExist';
else
    MT.elevation = 'ElevExist';
    % This becomes the reference coordinate, can be either positive or
    % negative.
    MT.refer = min(t1, 0);
    MT.znew = MT.cord.z - MT.refer;
end
fclose(fid);

% display information about data and stations
MT = DataSizeDisplay(MT, type, []);
if isfield(CData,'OutputFile')
    MT = DataSizeDisplay(MT, type, CData.OutputFile);
end
% storing all the data of the MT in the variable CData, hence all the next
% time I need to replace CData with CData.MT
CData.MT = MT;
end

function [codes, periods] = Sort_Station_Frequency(fname)
% Intermediary step for finding all the stations and time-periods. This
% is necessary while working for real data set with more than 1 data type
% with different stations and different set of frequencies for each data
% set
if exist(fname)
    fid = fopen(fname,'r');  % read file
else
    fprintf('%s %s \n','No file exist by the name : ',fname)
    error('Stop');
end
id = 0;
while(1)
    % read header
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    tmp = char(tmp{1});
    if isempty(tmp);
        break; 
    end
    id = id + 1;
    
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    blockinfo = textscan(fid,'> %s',6,'delimiter','\n');
    blockinfo = char(blockinfo{1});
    dataType = blockinfo(1,:);   CC.datatype{id} = dataType; 
    tmp  = sscanf(blockinfo(6,:),'%d',2);
    CC.nfreq{id} = tmp(1);                               % nfreq
    CC.nobs{id} = tmp(2);                                % nobs
    switch strtrim(dataType)
        case {'Off_Diagonal_Impedance','Full_Impedance','Full_Vertical_Components'}
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            CC.T{id} = data{1};
            CC.codes{id} = data{2};
        case {'Off_Diagonal_Rho_Phase','Phase_Tensor'}
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');
            CC.T{id} = data{1};
            CC.codes{id} = data{2};
        otherwise
            disp('Unknown data type');
            break;
            
    end
end
clear codes periods
if id == 2
    codes = [CC.codes{1}; CC.codes{2}];
    periods = [CC.T{1}; CC.T{2}];
    % Picking only the unique ones
    codes = sortrows(strtrim(char(unique(codes))));
    periods = unique(periods);
else
    codes = sortrows(strtrim(char(unique(data{2}))));
    periods = unique(data{1});
end
fclose(fid);
end

function MT = DataSizeDisplay(MT, type, logfile)
% Display Stations and data size on the command window oe write it in a log
% file, which ever is suited.

fprintf('\n');
if nargin < 3 || isempty(logfile)
    fprintf('%s \n','======== MT Data Information ========');
    fprintf('%s %d \n',' Total number of stations are    : ', MT.nobs);
    fprintf('%s %d \n',' Total number of frequencies are : ', MT.nfreq);
    fprintf('%s %7.3E %7.3E \n','  Min and Max periods are  : ',MT.T(1),MT.T(end));
else
    fid = fopen(logfile,'a');
    fprintf(fid,'%s \n','======== MT Data Information ========');
    fprintf(fid,'%s %d \n',' Total number of stations are    : ', MT.nobs);
    fprintf(fid,'%s %d \n',' Total number of frequencies are : ', MT.nfreq);
    fprintf(fid,'%s %7.3E %7.3E \n','Min and Max periods are  : ',MT.T(1),MT.T(end));
    fclose(fid);
end
if ~strcmp(type,'ReadData')
    fprintf('%s %d \n','Total data types are  : ',MT.NoofDataTypes);
    for i = 1:MT.NoofDataTypes
        fprintf('%s %s\n','      ',strtrim(char(MT.datatype{i})));
    end
end

% Print information about the data types
if strcmp(type,'ReadData')
    t = zeros(2,1);
    for i = 1:MT.NoofDataTypes
        datatype = strtrim(char(MT.datatype{i}));
        switch(datatype)
            case 'Off_Diagonal_Impedance'
                temp(:,1:2,:) = MT.id(:,2:3,:);  t(i) = 2*sum(temp(:));  clear temp
            case 'Full_Impedance'
                temp(:,1:4,:) = MT.id(:,1:4,:);  t(i) = 2*sum(temp(:));  clear temp
            case{'Off_Diagonal_Rho_Phase','Phase_Tensor'}
                temp(:,1:4,:) = MT.id(:,1:4,:);  t(i) = sum(temp(:));    clear temp
            case 'Full_Vertical_Components'
                temp(:,1:2,:) = MT.id(:,5:6,:);  t(i) = 2*sum(temp(:));  clear temp
        end
    end
    % Check where to print
    if nargin < 3 || isempty(logfile)
        if MT.conjugate == 1
            fprintf('%s \n',' -- > Data File has time dependency of exp{+iwt}');
            fprintf('%s \n',' -- > Complex data read as complex comjugate');
        end
        fprintf('%s %d \n','  Total data types are  : ',MT.NoofDataTypes);
        for i = 1:MT.NoofDataTypes
            datatype = strtrim(char(MT.datatype{i}));
            switch(datatype)
                case 'Off_Diagonal_Impedance'
                    fprintf('%d %s \n',t(i),' Off_Diagonal_Impedance read from the file ');
                case 'Full_Impedance'
                    fprintf('%d %s \n',t(i),' Full_Impedance read from the file ');
                case 'Off_Diagonal_Rho_Phase'
                    fprintf('%d %s \n',t(i),' Off_Diagonal_Rho_Phase read from the file ');
                case 'Phase_Tensor'
                    fprintf('%d %s \n',t(i),' Phase_Tensor read from the file ');
                case 'Full_Vertical_Components'
                    fprintf('%d %s \n',t(i),' Full_Vertical_Components read from the file ');
            end
        end
        for i = 1:MT.NoofDataTypes
            datatype = strtrim(char(MT.datatype{i}));
            switch(datatype)
                case {'Off_Diagonal_Impedance','Full_Impedance'}
                    typeUnits = strtrim(char(MT.typeUnits{1}));
                    switch(typeUnits)
                        case '[V/m]/[T]'
                            fprintf('%s \n',' -- > No unit converting required');
                        case '[mV/km]/[nT]'
                            fprintf('%s \n',' -- > Converting Impedance from [mV/km]/[nT] to [V/m]/[T]');
                        case {'[V/m]/[A/m]','Ohm'}
                            fprintf('%s \n',' -- > Converting Impedance from [V/m]/[A/m] or Ohm to [V/m]/[T]');
                        otherwise
                            fprintf('%s \n',' -- > Error in specifying units for Impedance');
                            fprintf('%s \n',' -- > Converting to default [V/m]/[T]');
                    end
            end
        end
        MT.Ndata = t(1)+t(2);    % Total number of independent data points
        MT.Ndata_impedance = t(1);
        fprintf('%s %d \n',' Total independent data points are : ', MT.Ndata);
    else
        fid = fopen(logfile,'a');
        if MT.conjugate == 1
            fprintf(fid,'%s \n',' -- > Data File has time dependency of exp{+iwt}');
            fprintf(fid,'%s \n',' -- > Complex data read as complex comjugate');
        end
        fprintf(fid,'%s %d \n','  Total data types are  : ',MT.NoofDataTypes);
        for i = 1:MT.NoofDataTypes
            datatype = strtrim(char(MT.datatype{i}));
            switch(datatype)
                case 'Off_Diagonal_Impedance'
                    fprintf(fid,'%d %s \n',t(i),' Off_Diagonal_Impedance read from the file ');
                case 'Full_Impedance'
                    fprintf(fid,'%d %s \n',t(i),' Full_Impedance read from the file ');
                case 'Off_Diagonal_Rho_Phase'
                    fprintf(fid,'%d %s \n',t(i),' Off_Diagonal_Rho_Phase read from the file ');
                case 'Phase_Tensor'
                    fprintf(fid,'%d %s \n',t(i),' Phase_Tensor read from the file ');
                case 'Full_Vertical_Components'
                    fprintf(fid,'%d %s \n',t(i),' Full_Vertical_Components read from the file ');
            end
        end
        for i = 1:MT.NoofDataTypes
            datatype = strtrim(char(MT.datatype{i}));
            switch(datatype)
                case {'Off_Diagonal_Impedance','Full_Impedance'}
                    typeUnits = strtrim(char(MT.typeUnits{1}));
                    switch(typeUnits)
                        case '[V/m]/[T]'
                            fprintf(fid,'%s \n',' -- > No Impedance convertion required');
                        case '[mV/km]/[nT]'
                            fprintf(fid,'%s \n',' -- > Converting Impedance from [mV/km]/[nT] to [V/m]/[T]');
                        case {'[V/m]/[A/m]','Ohm'}
                            fprintf(fid,'%s \n',' -- > Converting Impedance from [V/m]/[A/m] or Ohm to [V/m]/[T]');
                        otherwise
                            fprintf(fid,'%s \n',' -- > Error in specifying units for Impedance');
                            fprintf(fid,'%s \n',' -- > Converting to default [V/m]/[T]');
                    end
            end
        end
        MT.Ndata = t(1)+t(2);    % Total number of independent data points
        fprintf(fid,'%s %d \n',' Total independent data points are : ', MT.Ndata);
        fclose(fid);
    end
end
end