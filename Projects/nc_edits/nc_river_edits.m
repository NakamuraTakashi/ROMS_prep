NC_FILE1 = './tokyobay_grid/TokyoBay3_river_2d_2018_2020.nc';
ncdisp(NC_FILE1)
%% 
NC_FILE2 = './tokyobay_grid/TokyoBay3_river_NZ30_2018_2020v1_aong.nc';
ncdisp(NC_FILE2)
%% 
NC_FILE3 = './tokyobay_grid/test2.nc';
ncdisp(NC_FILE3)
%% 
NC_VAR = 'river_time';
nctime2 = ncread(NC_FILE2,NC_VAR);
nctime3 = ncread(NC_FILE3,NC_VAR);
ncwrite(NC_FILE1,NC_VAR,nctime2);
nctime1 = ncread(NC_FILE1,NC_VAR);
%% 
NC_VAR = 'river_Xposition';
ncdata2 = ncread(NC_FILE2,NC_VAR);
ncdata3 = ncread(NC_FILE3,NC_VAR);
% ncdata2 = squeeze(ncdata2);
% ncwrite(NC_FILE1,NC_VAR,ncdata2);
% ncdata1 = ncread(NC_FILE1,NC_VAR);
%% 
NC_VAR = 'river_Eposition';
ncdata2 = ncread(NC_FILE2,NC_VAR);
ncdata3 = ncread(NC_FILE3,NC_VAR);
% ncdata2 = squeeze(ncdata2);
% ncwrite(NC_FILE1,NC_VAR,ncdata2);
% ncdata1 = ncread(NC_FILE1,NC_VAR);
%% 
NC_VAR = 'river_direction';
ncdata2 = ncread(NC_FILE2,NC_VAR);
% ncdata3 = ncread(NC_FILE3,NC_VAR);
ncdata2 = squeeze(ncdata2);
ncwrite(NC_FILE1,NC_VAR,ncdata2);
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_Vshape';
ncdata2 = ncread(NC_FILE2,NC_VAR);
% ncdata3 = ncread(NC_FILE3,NC_VAR);
ncdata2 = squeeze(ncdata2);
ncwrite(NC_FILE1,NC_VAR,ncdata2);
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_transport';
ncdata2 = ncread(NC_FILE2,NC_VAR);
ncdata3 = ncread(NC_FILE3,NC_VAR);
% ncdata2 = squeeze(ncdata2);
% ncwrite(NC_FILE1,NC_VAR,ncdata2);
% ncdata1 = ncread(NC_FILE1,NC_VAR);


%% 
NC_VAR = 'river_temp';
ncdata2 = ncread(NC_FILE2,NC_VAR, [1 1 1],[Inf 1 Inf]);
% ncdata3 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
ncwrite(NC_FILE1,NC_VAR,ncdata2);
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_salt';
ncdata2 = ncread(NC_FILE2,NC_VAR, [1 1 1],[Inf 1 Inf]);
% ncdata3 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
ncwrite(NC_FILE1,NC_VAR,ncdata2);
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_Oxyg';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
%% 

NC_VAR = 'river_DO';
% ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncwrite(NC_FILE1,NC_VAR,ncdata1);
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_NO3';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_NO3_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_NH4';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_NH4_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_PO4';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_PO4_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_DOC_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_DOC01_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_DOC_02';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_DOC02_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_DON_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_DON01_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_DON_02';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_DON02_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_DOP_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_DOP01_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_DOP_02';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_DOP02_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);


%% 
NC_VAR = 'river_POC_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_POC01_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_POC_02';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_POC02_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_PON_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_PON01_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_PON_02';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_PON02_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_POP_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_POP01_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

NC_VAR = 'river_POP_02';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2);
NC_VAR = 'river_POP02_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_POC_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2)*0.01;
NC_VAR = 'river_POC03_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);


NC_VAR = 'river_PON_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2)*0.01;
NC_VAR = 'river_PON03_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);


NC_VAR = 'river_POP_01';
ncdata2 = ncread(NC_FILE3,NC_VAR, [1 1 1],[Inf 1 Inf]);
ncdata2 = squeeze(ncdata2)*0.01;
NC_VAR = 'river_POP03_01';
ncwrite(NC_FILE1,NC_VAR,ncdata2(1:9,:));
ncdata1 = ncread(NC_FILE1,NC_VAR);

%% 
NC_VAR = 'river_TA';
ncdata1 = ncread(NC_FILE1,NC_VAR);
ncdata1 = ones(size(ncdata1))*400;
ncwrite(NC_FILE1,NC_VAR,ncdata1);

%% 
NC_VAR = 'river_DIC_01';
ncdata1 = ncread(NC_FILE1,NC_VAR);
ncdata1 = ones(size(ncdata1))*400;
ncwrite(NC_FILE1,NC_VAR,ncdata1);


%% 
NC_VAR = 'river_transport';
ncdata2 = ncread(NC_FILE2,NC_VAR);
ncdata3 = ncread(NC_FILE3,NC_VAR);
%% 

plot(nctime1,ncdata1(4,:))