function map = mapcomb(siz)

map = [
    0.7 linspace(0,.95,19);
    0.1 linspace(0,.95,19);
    0.1 linspace(0.7,1,19)]';
map = flipud(map);
