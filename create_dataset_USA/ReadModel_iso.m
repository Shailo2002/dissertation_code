function [rho, x, y, z, a, b, c] = ReadModel_iso(file)   

fileID = fopen(file);
Line = fgets(fileID);   % comment  Line for model
Line = fgets(fileID);   Line = strtrim(Line);
[T] = sscanf(Line, '%d',[3 1]);
if findstr(Line,'LOGE')
    type = 'LOGE';
else
    type = 'LINEAR';
end

nx = T(1); ny = T(2); nz = T(3); 
a = fscanf(fileID,'%f',nx);
b = fscanf(fileID,'%f',ny);
c = fscanf(fileID,'%f',nz);
rho = zeros(nx, ny, nz);
for iz = 1:nz
    for iy = 1:ny
        for ix = nx:-1:1
            rho(ix,iy,iz) = fscanf(fileID,'%f',1);
        end
    end
end
origin = fscanf(fileID,'%f',3);

x = [origin(1);  origin(1)+cumsum(a)];
y = [origin(2);  origin(2)+cumsum(b)];
z = [origin(3);  origin(3)+cumsum(c)];
   
fclose(fileID);
end