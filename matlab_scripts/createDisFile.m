function createDisFile(DZ, DX, phanSizeZ, phanSizeX, fname)
%%INPUT:
%%DZ = The axial displacements
%%DX = The lateral displacements
%%phanSizeX = the size of the phantom in the lateral direction
%%phanSizeZ = the size of the phantom in the axial direction
%%
%%NOTE: 1)The units are unimportant so long as you are consistent
%%	2)This function expects square matrices as inputs

if(size(DZ,1) ~= size(DZ,2) || size(DX,1) ~= size(DX,2) )
disp('Error.  Square matrices expected')
end

temp = size(DZ,1);
phanSizeZ = 40;
phanSizeX = 40;

h = fopen(fname, 'wb');
fwrite(h, temp, 'int32');
fwrite(h, DX/phanSizeX, 'double');
fwrite(h, DZ/phanSizeZ, 'double');

fclose(h);
