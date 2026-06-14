%set(0,'DefaultAxesFontName','SansSerif')
set(0,'DefaultAxesFontName','Ubuntu')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultLineMarkerSize',20)
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
set(groot,'defaultFigureRenderer','painters')
v_order = [7 6 3 1 5 2];
cmap = [0.619608000000000	0.00392200000000000	0.258824000000000;
0.916340000000000	0.366013000000000	0.278431333333333;
0.993464000000000	0.747712333333333	0.435294000000000;
1	1	0.749020000000000;
0.747712333333333	0.898039333333333	0.627450666666667;
0.332026000000000	0.684967000000000	0.678431333333333;
0.368627000000000	0.309804000000000	0.635294000000000];
set(0,'DefaultAxesColorOrder',cmap(v_order,:));
%clear