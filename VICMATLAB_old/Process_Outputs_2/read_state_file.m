% Function for reading in the VIC state file
%
% See VIC classic driver documentation
%
% This could be more generalized. Customize to suit your needs, depending
% on snow, vegetation bands.
%
% Something weird is going on. There are more columns in my state file than
% there are in the documentation. Perhaps it is out of date.
%
% Making this reader is not a fruitful use of my time right now. This VIC
% data assimilation using classic mode project will have to stop. Perhaps I
% can do it in image mode more easily.

function state = read_state_file(statefileprefix)

% statefilename = state_file_name{1};

state_m1_i_1_19990102_00000

y1 = num2str(year(end_date));
m1 = num2str(month(end_date));
d1 = num2str(day(end_date));

if str2double(d1)<10
    d1 = ['0' d1];
end
   
if str2double(m1)<10
    m1 = ['0' m1];
end

statefilename = [statefileprefix '_' y1 m1 d1 '_00000'];

fID = fopen(statefilename, 'r');

% Read the first line
tline = fgetl(fID);
tmp1 = str2num(tline);
y2 = tmp1(1);
m2 = tmp1(2);
d2 = tmp1(3);

% Read the second line
tline = fgetl(fID);
tmp2 = str2num(tline);
n_layers = tmp2(1);
n_nodes = tmp2(2);

% Read the third line
tline = fgetl(fID);
tmp3 = str2num(tline);
cell_number = tmp3(1);
n_veg = tmp3(2);
n_bands = tmp3(3);
dz_node = zeros(n_nodes, 1);
node_depth = zeros(n_nodes, 1);
for k=1:n_nodes
    dz_node(k) = tmp3(3+k);  
    node_depth(k) = tmp3(6+k);
end

% Read all the other lines
for v=1:n_veg

    tline = fgetl(fID);
    tmp4 = str2num(tline);

    veg = tmp4(1);
    band = tmp4(2);

    moist = zeros(n_layers, 1);
    ice = zeros(n_layers, 1);
    for k=1:n_layers
        moist(k) = tmp4(2+k);
        ice(k) = tmp4(5+k);
    end

    wdew = tmp4(5+n_layers+1);
    last_snow = tmp4(5+n_layers+2);
    melting = tmp4(5+n_layers+3);
    coverage = tmp4(5+n_layers+4);
    swq = tmp4(5+n_layers+5);
    surf_temp = tmp4(5+n_layers+6);
    surf_water = tmp4(5+n_layers+7);
    pack_temp = tmp4(5+n_layers+8);
    pack_water = tmp4(5+n_layers+9);
    density = tmp4(5+n_layers+10);
    coldcontent = tmp4(5+n_layers+11);
    snow_canopy = tmp4(5+n_layers+12);

    node_T = zeros(n_nodes, 1);
    for k=1:n_nodes
        node_T(k) = tmp4(5+n_layers+12+k);
    end

end


fclose(fID)



return