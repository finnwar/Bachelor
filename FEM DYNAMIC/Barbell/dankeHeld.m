%% Dr. Held Dankeschön

x = linspace(0,1,61);
d = -60:0;

aoH = 1/(exp(1^2)-1);
amH = 0.1/(exp(1^2)-1);

ohneHeld = aoH*exp(x.^2)-aoH;
mitHeld = amH*exp(x.^2)-amH;

col.blue = [0 100 222]/255;
col.green = [0 140 0]/255;
col.red = [220 33 77]/255;
col.purple = [102 0 102]/255;
col.orange = [255 102 0]/255;
col.teal = [55 200 171]/255;
col.black = [0 0 0];
col.darkGrey = [70 70 70]/255;
col.grey = [150 150 150]/255;
col.lightGrey = [225 225 225]/255;
col.mumBlueLogo = [24 59 101]/255;
col.mumTealLogo = [45 198 214]/255;
col.white = [1 1 1];