function[ CData, AllData ] = ReadStatData(fname, CData, type)
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

CData.mu0 = 4*pi*1.0E-07;
CData.ncomp = 0;

% Quickly read the data file for finding all unique stations and periods. 
[codes, periods] = Sort_Station_Frequency(fname);
CData.nfreq = length(periods);    % Maximum number of periods
CData.nobs = length(codes);       % Maximum number of Stations
CData.T = periods;
% Initialising AllData, complete data set
AllData.id = ones(CData.nobs,6,CData.nfreq);
AllData.data = zeros(CData.nobs,6,CData.nfreq);
AllData.err = ones(CData.nobs,6,CData.nfreq);

fid = fopen(fname,'r');  % read file
id = 0;
while(1)
    
    id = id + 1;
    % read header
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    tmp = char(tmp{1});
    if isempty(tmp); break; end
    CData.header = tmp;                                % header
    % read the block header
    blockinfo = textscan(fid,'# %s',1,'delimiter','\n');
    header = char(blockinfo{1});
    blockinfo = textscan(fid,'> %s',6,'delimiter','\n');
    blockinfo = char(blockinfo{1});
    dataType = blockinfo(1,:);
    signstr  = blockinfo(2,:);
    if findstr(signstr,'-')
        CData.isign = -1;
        CData.ii = complex(0,-1);
        CData.conjugate = 0;
    else  % Need to replace complex data as complex conjugate
        CData.isign = -1;
        CData.ii = complex(0,-1);        
        CData.conjugate = 1;
        % CData.isign = +1;
        % CData.ii = complex(0,+1);
    end
    CData.signstr = signstr;                            % field type
    % need copy only for impedance units, tipper is dimension less.
    typeUnits = strtrim(blockinfo(3,:));                % units
    CData.orientation = sscanf(blockinfo(4,:),'%f',1);  % orientation
    CData.origin = sscanf(blockinfo(5,:),'%f',3);       % origin
    tmp  = sscanf(blockinfo(6,:),'%d',2);               % nfreq and nobs, discarded now

    % All datatype and their corresponding units stored here, in cells
    CData.Nfreq{id} = tmp(1);
    CData.Nobs{id} = tmp(2);
    CData.datatype{id} = dataType;      
    CData.typeUnits{id} = typeUnits;
    CData.lineinfo{id} = header;
    
    switch strtrim(dataType)
        case 'Off_Diagonal_Impedance'
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            ncomp = 2;  CData.ncomp = CData.ncomp + ncomp;
            comp = ['ZXY';'ZYX'];   CData.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            CData.cord.x(k,1) = data{5}(ind);
                            CData.cord.y(k,1) = data{6}(ind);
                            CData.cord.z(k,1) = data{7}(ind);
                            CData.Stat.code(k,:) = data{2}(ind);
                            CData.Stat.lat(k,1) = data{3}(ind);
                            CData.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                if CData.conjugate == 1   % Data is complex conjugate
                                    AllData.data(k,i+1,j) = complex(data{9}(ind), -data{10}(ind));
                                else
                                   AllData.data(k,i+1,j) = complex(data{9}(ind), data{10}(ind)); 
                                end
                                AllData.err(k,i+1,j) = data{11}(ind);
                                if abs(AllData.data(k,i+1,j)) == 0
                                    AllData.id(k,i+1,j) = 0;
                                end
                            end
                        else
                            AllData.id(k,i+1,j) = 0;
                        end
                    end
                end
            end
        case 'Full_Impedance'                    % Reading Impedance
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            ncomp = 4;  CData.ncomp = CData.ncomp + ncomp;
            comp = ['ZXX';'ZXY';'ZYX';'ZYY'];   CData.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            CData.Stat.lat(k,1) = data{3}(ind);
                            CData.Stat.lon(k,1) = data{4}(ind);
                            CData.Stat.code(k,:) = data{2}(ind);
                            CData.cord.x(k,1) = data{5}(ind);
                            CData.cord.y(k,1) = data{6}(ind);
                            CData.cord.z(k,1) = data{7}(ind);
                            if strcmp(type,'ReadData');
                                if CData.conjugate == 1   % Data is complex conjugate
                                    AllData.data(k,i,j) = complex(data{9}(ind), -data{10}(ind));
                                else
                                   AllData.data(k,i,j) = complex(data{9}(ind), data{10}(ind)); 
                                end
                                AllData.err(k,i,j) = data{11}(ind);
                                if abs(AllData.data(k,i,j)) == 0
                                    AllData.id(k,i,j) = 0;
                                end
                            end
                        else
                            AllData.id(k,i,j) = 0;
                        end
                    end
                end
            end
        case 'Off_Diagonal_Rho_Phase'                    % Apparent Resistivity and phase, all data set is rael
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');
            ncomp = 4;  CData.ncomp = CData.ncomp + ncomp;
            comp = ['RHOXY';'PHSXY';'RHOYX';'PHSYX'];   CData.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            CData.cord.x(k,1) = data{5}(ind);
                            CData.cord.y(k,1) = data{6}(ind);
                            CData.cord.z(k,1) = data{7}(ind);
                            CData.Stat.code(k,:) = data{2}(ind);
                            CData.Stat.lat(k,1) = data{3}(ind);
                            CData.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                AllData.data(k,i,j) = data{9}(ind);
                                AllData.err(k,i,j) = data{10}(ind);
                                if abs(AllData.data(k,i,j)) == 0
                                    AllData.id(k,i,j) = 0;
                                end
                            end
                        else
                            AllData.id(k,i,j) = 0;
                        end
                    end
                end
            end
        case 'Phase_Tensor'                    % Phase Tensor, all data set is rael
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');
            ncomp = 4;  CData.ncomp = CData.ncomp + ncomp;
            comp = ['PTXX';'PTXY';'PTYX';'PTYY'];   CData.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            CData.cord.x(k,1) = data{5}(ind);
                            CData.cord.y(k,1) = data{6}(ind);
                            CData.cord.z(k,1) = data{7}(ind);
                            CData.Stat.code(k,:) = data{2}(ind);
                            CData.Stat.lat(k,1) = data{3}(ind);
                            CData.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                AllData.data(k,i,j) = data{9}(ind);
                                AllData.err(k,i,j) = data{10}(ind);
                                if abs(AllData.data(k,i,j)) == 0
                                    AllData.id(k,i,j) = 0;
                                end
                            end
                        else
                            AllData.id(k,i,j) = 0;
                        end
                    end
                end
            end
        case 'Full_Vertical_Components'          % Reading Tipper                           
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            ncomp = 2;  CData.ncomp = CData.ncomp + ncomp;
            comp = ['TX';'TY'];     CData.comp{id} = comp;
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp	% ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            CData.cord.x(k,1) = data{5}(ind);
                            CData.cord.y(k,1) = data{6}(ind);
                            CData.cord.z(k,1) = data{7}(ind);
                            CData.Stat.code(k,:) = data{2}(ind);
                            CData.Stat.lat(k,1) = data{3}(ind);
                            CData.Stat.lon(k,1) = data{4}(ind);
                            if strcmp(type,'ReadData');
                                if CData.conjugate == 1   % Data is complex conjugate
                                    AllData.data(k,i+4,j) = complex(data{9}(ind), -data{10}(ind));
                                else
                                   AllData.data(k,i+4,j) = complex(data{9}(ind), data{10}(ind)); 
                                end
                                AllData.err(k,i+4,j) = data{11}(ind);
                                if abs(AllData.data(k,i+4,j)) == 0
                                    AllData.id(k,i+4,j) = 0;
                                end
                            end
                        else
                            AllData.id(k,i+4,j) = 0;
                        end
                    end
                end
            end
        otherwise
            disp('Unknown data type');
            break;
    end
end
CData.omega = 2*pi./CData.T;    % Omega's
CData.id = AllData.id;          % Storing component in the CData

% Unit conversion introduced here. we will convert the data to [v/m]/[T]
% as we work in these units and later if required, re convert to the old
% units, unit conversion for data and error as well and only for impedance
% and not for dimension less quantity like VTF and Phase Tensor
CData.NoofDataTypes = length(CData.datatype);  % Number of data types

CData.f1 = -1/CData.ii;
CData.factor = CData.mu0;
CData.factor1 = 1;  CData.factor2 = 1;

for i = 1:CData.NoofDataTypes
    datatype = strtrim(char(CData.datatype{i}));
    switch(datatype)
        case {'Off_Diagonal_Impedance','Full_Impedance'}
            typeUnits = strtrim(char(CData.typeUnits{1}));
            switch(typeUnits)
                case '[V/m]/[T]'
                    CData.factor1 = 1;
                    CData.factor2 = 1;
                case '[mV/km]/[nT]'
                    CData.factor1 = 1000;
                    CData.factor2 = 1/1000;
                case {'[V/m]/[A/m]','Ohm'}
                    CData.factor1 = 1/CData.mu0;
                    CData.factor2 = CData.mu0;
            end
            % Convert data to desired units is required
            if (CData.factor1 ~= 1) && strcmp(type,'ReadData')
                AllData.data(:,1:4,:) = CData.factor1 * AllData.data(:,1:4,:);
                AllData.err(:,1:4,:) = CData.factor1 * AllData.err(:,1:4,:);
            end
        case 'Off_Diagonal_Rho_Phase'
            % Convert rho to log10(rho) if so desired
            if (strcmp(type,'ReadData') && strcmp(CData.rho,'log10'))
                AllData.data(:,1,:) = log10(AllData.data(:,1,:));
                AllData.data(:,3,:) = log10(AllData.data(:,3,:));
                AllData.err(:,1,:) = log10(AllData.err(:,1,:));
                AllData.err(:,3,:) = log10(AllData.err(:,3,:));
            end
    end
end

% checking for elevation if exists;
t1 = min(CData.cord.z);
t2 = max(CData.cord.z);
if t1 ~= t2
    CData.elevation = 'ElevExist';
    CData.refer = t2;           % This becomes the reference coordinate, -ve in case of elevation
else
    CData.elevation = 'NoElevExist';
    CData.refer = 0.0;          % This becomes the reference coordinate, 0 in case of no elevation
end
fclose(fid);

% display information about data and stations
CData = DataSizeDisplay(CData, type, []);
if isfield(CData,'OutputFile')
    CData = DataSizeDisplay(CData, type, CData.OutputFile);
end
end

function [codes, periods] = Sort_Station_Frequency(fname)
% Intermediary step for finding all the stations and time-periods. This
% is necessary while working for real data set with more than 1 data type
% with different stations and different set of frequencies for each data
% set
fid = fopen(fname,'r');  % read file
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
fprintf('%s %d \n','Total number of unique stations are : ',length(codes));
fprintf('%s %d \n','Total number of unique periods are  : ',length(periods));
if id == 2
    fprintf('%s \n',' -- > More than 1 data types exists in the data file');
else
    fprintf('%s \n',' -- > Only 1 data types exists in the data file');
end
fclose(fid);
end

function CData = DataSizeDisplay(CData, type, logfile)
% Display Stations and data size on the command window oe write it in a log
% file, which ever is suited.
if nargin < 3 || isempty(logfile)
%     fprintf('%s %d \n',' Total number of stations are    : ', CData.nobs);
%     fprintf('%s %d \n',' Total number of frequencies are : ', CData.nfreq);
else
    fid = fopen(logfile,'a');
    fprintf(fid,'%s %d \n',' Total number of stations are    : ', CData.nobs);
    fprintf(fid,'%s %d \n',' Total number of frequencies are : ', CData.nfreq);
    fclose(fid);
end

% Print information about the data types
if strcmp(type,'ReadData')
    t = zeros(2,1);
    for i = 1:CData.NoofDataTypes
        datatype = strtrim(char(CData.datatype{i}));
        switch(datatype)
            case 'Off_Diagonal_Impedance'
                temp(:,1:2,:) = CData.id(:,2:3,:);  t(i) = 2*sum(temp(:));  clear temp
            case 'Full_Impedance'
                temp(:,1:4,:) = CData.id(:,1:4,:);  t(i) = 2*sum(temp(:));  clear temp
            case{'Off_Diagonal_Rho_Phase','Phase_Tensor'}
                temp(:,1:4,:) = CData.id(:,1:4,:);  t(i) = sum(temp(:));    clear temp
            case 'Full_Vertical_Components'
                temp(:,1:2,:) = CData.id(:,5:6,:);  t(i) = 2*sum(temp(:));  clear temp
        end
    end
    % Check where to print
    if nargin < 3 || isempty(logfile)
        if CData.conjugate == 1
            fprintf('%s \n',' -- > Data File has time dependency of exp{+iwt}');
            fprintf('%s \n',' -- > Complex data read as complex comjugate');
        end
        for i = 1:CData.NoofDataTypes
            datatype = strtrim(char(CData.datatype{i}));
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
        for i = 1:CData.NoofDataTypes
            datatype = strtrim(char(CData.datatype{i}));
            switch(datatype)
                case {'Off_Diagonal_Impedance','Full_Impedance'}
                    typeUnits = strtrim(char(CData.typeUnits{1}));
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
        CData.Ndata = t(1)+t(2);    % Total number of independent data points
        fprintf('%s %d \n',' Total independent data points are : ', CData.Ndata);
    else
        fid = fopen(logfile,'a');
        if CData.conjugate == 1
            fprintf(fid,'%s \n',' -- > Data File has time dependency of exp{+iwt}');
            fprintf(fid,'%s \n',' -- > Complex data read as complex comjugate');
        end
        for i = 1:CData.NoofDataTypes
            datatype = strtrim(char(CData.datatype{i}));
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
        for i = 1:CData.NoofDataTypes
            datatype = strtrim(char(CData.datatype{i}));
            switch(datatype)
                case {'Off_Diagonal_Impedance','Full_Impedance'}
                    typeUnits = strtrim(char(CData.typeUnits{1}));
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
        CData.Ndata = t(1)+t(2);    % Total number of independent data points
        fprintf(fid,'%s %d \n',' Total independent data points are : ', CData.Ndata);
        fclose(fid);
    end
end
end