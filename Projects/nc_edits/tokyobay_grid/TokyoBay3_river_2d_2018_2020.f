      program main
      include 'netcdf.inc'
* error status return
      integer stat
* netCDF ncid
      integer  ncid

* dimension lengths
      integer 's_rho'_len
      parameter ('s_rho'_len = 30)
      integer 'river'_len
      parameter ('river'_len = 9)
      integer 'river_time'_len
      parameter ('river_time'_len = NF_UNLIMITED)
* dimension ids
      integer 's_rho'_dim
      integer 'river'_dim
      integer 'river_time'_dim

* variable ids
      integer 'river'_id;
      integer 'river_Xposition'_id;
      integer 'river_Eposition'_id;
      integer 'river_direction'_id;
      integer 'river_Vshape'_id;
      integer 'river_time'_id;
      integer 'river_transport'_id;
      integer 'river_salt'_id;
      integer 'river_temp'_id;
      integer 'river_DO'_id;
      integer 'river_TA'_id;
      integer 'river_DIC_01'_id;
      integer 'river_NO3_01'_id;
      integer 'river_NH4_01'_id;
      integer 'river_PO4_01'_id;
      integer 'river_DOC01_01'_id;
      integer 'river_DOC02_01'_id;
      integer 'river_DON01_01'_id;
      integer 'river_DON02_01'_id;
      integer 'river_DOP01_01'_id;
      integer 'river_DOP02_01'_id;
      integer 'river_POC01_01'_id;
      integer 'river_POC02_01'_id;
      integer 'river_POC03_01'_id;
      integer 'river_PON01_01'_id;
      integer 'river_PON02_01'_id;
      integer 'river_PON03_01'_id;
      integer 'river_POP01_01'_id;
      integer 'river_POP02_01'_id;
      integer 'river_POP03_01'_id;

* rank (number of dimensions) for each variable
      integer 'river'_rank
      parameter ('river'_rank = 1)
      integer 'river_Xposition'_rank
      parameter ('river_Xposition'_rank = 1)
      integer 'river_Eposition'_rank
      parameter ('river_Eposition'_rank = 1)
      integer 'river_direction'_rank
      parameter ('river_direction'_rank = 1)
      integer 'river_Vshape'_rank
      parameter ('river_Vshape'_rank = 2)
      integer 'river_time'_rank
      parameter ('river_time'_rank = 1)
      integer 'river_transport'_rank
      parameter ('river_transport'_rank = 2)
      integer 'river_salt'_rank
      parameter ('river_salt'_rank = 2)
      integer 'river_temp'_rank
      parameter ('river_temp'_rank = 2)
      integer 'river_DO'_rank
      parameter ('river_DO'_rank = 2)
      integer 'river_TA'_rank
      parameter ('river_TA'_rank = 2)
      integer 'river_DIC_01'_rank
      parameter ('river_DIC_01'_rank = 2)
      integer 'river_NO3_01'_rank
      parameter ('river_NO3_01'_rank = 2)
      integer 'river_NH4_01'_rank
      parameter ('river_NH4_01'_rank = 2)
      integer 'river_PO4_01'_rank
      parameter ('river_PO4_01'_rank = 2)
      integer 'river_DOC01_01'_rank
      parameter ('river_DOC01_01'_rank = 2)
      integer 'river_DOC02_01'_rank
      parameter ('river_DOC02_01'_rank = 2)
      integer 'river_DON01_01'_rank
      parameter ('river_DON01_01'_rank = 2)
      integer 'river_DON02_01'_rank
      parameter ('river_DON02_01'_rank = 2)
      integer 'river_DOP01_01'_rank
      parameter ('river_DOP01_01'_rank = 2)
      integer 'river_DOP02_01'_rank
      parameter ('river_DOP02_01'_rank = 2)
      integer 'river_POC01_01'_rank
      parameter ('river_POC01_01'_rank = 2)
      integer 'river_POC02_01'_rank
      parameter ('river_POC02_01'_rank = 2)
      integer 'river_POC03_01'_rank
      parameter ('river_POC03_01'_rank = 2)
      integer 'river_PON01_01'_rank
      parameter ('river_PON01_01'_rank = 2)
      integer 'river_PON02_01'_rank
      parameter ('river_PON02_01'_rank = 2)
      integer 'river_PON03_01'_rank
      parameter ('river_PON03_01'_rank = 2)
      integer 'river_POP01_01'_rank
      parameter ('river_POP01_01'_rank = 2)
      integer 'river_POP02_01'_rank
      parameter ('river_POP02_01'_rank = 2)
      integer 'river_POP03_01'_rank
      parameter ('river_POP03_01'_rank = 2)

* variable shapes
      integer 'river'_dims('river'_rank)
      integer 'river_Xposition'_dims('river_Xposition'_rank)
      integer 'river_Eposition'_dims('river_Eposition'_rank)
      integer 'river_direction'_dims('river_direction'_rank)
      integer 'river_Vshape'_dims('river_Vshape'_rank)
      integer 'river_time'_dims('river_time'_rank)
      integer 'river_transport'_dims('river_transport'_rank)
      integer 'river_salt'_dims('river_salt'_rank)
      integer 'river_temp'_dims('river_temp'_rank)
      integer 'river_DO'_dims('river_DO'_rank)
      integer 'river_TA'_dims('river_TA'_rank)
      integer 'river_DIC_01'_dims('river_DIC_01'_rank)
      integer 'river_NO3_01'_dims('river_NO3_01'_rank)
      integer 'river_NH4_01'_dims('river_NH4_01'_rank)
      integer 'river_PO4_01'_dims('river_PO4_01'_rank)
      integer 'river_DOC01_01'_dims('river_DOC01_01'_rank)
      integer 'river_DOC02_01'_dims('river_DOC02_01'_rank)
      integer 'river_DON01_01'_dims('river_DON01_01'_rank)
      integer 'river_DON02_01'_dims('river_DON02_01'_rank)
      integer 'river_DOP01_01'_dims('river_DOP01_01'_rank)
      integer 'river_DOP02_01'_dims('river_DOP02_01'_rank)
      integer 'river_POC01_01'_dims('river_POC01_01'_rank)
      integer 'river_POC02_01'_dims('river_POC02_01'_rank)
      integer 'river_POC03_01'_dims('river_POC03_01'_rank)
      integer 'river_PON01_01'_dims('river_PON01_01'_rank)
      integer 'river_PON02_01'_dims('river_PON02_01'_rank)
      integer 'river_PON03_01'_dims('river_PON03_01'_rank)
      integer 'river_POP01_01'_dims('river_POP01_01'_rank)
      integer 'river_POP02_01'_dims('river_POP02_01'_rank)
      integer 'river_POP03_01'_dims('river_POP03_01'_rank)

* variable declarations
* attribute vectors
      integer textval(1)


* enter define mode
      stat = nf_create('TokyoBay3_river_2d_2018_2020.nc', nf_clobber, nc
     1id);
      call check_err(stat)
* define dimensions
      stat = nf_def_dim(ncid, 's_rho', 's_rho'_len, 's_rho'_dim);
      call check_err(stat)
      stat = nf_def_dim(ncid, 'river', 'river'_len, 'river'_dim);
      call check_err(stat)
      stat = nf_def_dim(ncid, 'river_time', 'river_time'_len, 'river_tim
     1e'_dim);
      call check_err(stat)
* define variables

      'river'_dims(1) = 'river'_dim
      stat = nf_def_var(ncid, 'river', nf_double, 'river'_rank, 'river'_
     1dims, 'river'_id);
      call check_err(stat)

      'river_Xposition'_dims(1) = 'river'_dim
      stat = nf_def_var(ncid, 'river_Xposition', nf_double, 'river_Xposi
     1tion'_rank, 'river_Xposition'_dims, 'river_Xposition'_id);
      call check_err(stat)

      'river_Eposition'_dims(1) = 'river'_dim
      stat = nf_def_var(ncid, 'river_Eposition', nf_double, 'river_Eposi
     1tion'_rank, 'river_Eposition'_dims, 'river_Eposition'_id);
      call check_err(stat)

      'river_direction'_dims(1) = 'river'_dim
      stat = nf_def_var(ncid, 'river_direction', nf_double, 'river_direc
     1tion'_rank, 'river_direction'_dims, 'river_direction'_id);
      call check_err(stat)

      'river_Vshape'_dims(1) = 'river'_dim
      'river_Vshape'_dims(2) = 's_rho'_dim
      stat = nf_def_var(ncid, 'river_Vshape', nf_double, 'river_Vshape'_
     1rank, 'river_Vshape'_dims, 'river_Vshape'_id);
      call check_err(stat)

      'river_time'_dims(1) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_time', nf_double, 'river_time'_rank
     1, 'river_time'_dims, 'river_time'_id);
      call check_err(stat)

      'river_transport'_dims(1) = 'river'_dim
      'river_transport'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_transport', nf_double, 'river_trans
     1port'_rank, 'river_transport'_dims, 'river_transport'_id);
      call check_err(stat)

      'river_salt'_dims(1) = 'river'_dim
      'river_salt'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_salt', nf_double, 'river_salt'_rank
     1, 'river_salt'_dims, 'river_salt'_id);
      call check_err(stat)

      'river_temp'_dims(1) = 'river'_dim
      'river_temp'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_temp', nf_double, 'river_temp'_rank
     1, 'river_temp'_dims, 'river_temp'_id);
      call check_err(stat)

      'river_DO'_dims(1) = 'river'_dim
      'river_DO'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DO', nf_double, 'river_DO'_rank, 'r
     1iver_DO'_dims, 'river_DO'_id);
      call check_err(stat)

      'river_TA'_dims(1) = 'river'_dim
      'river_TA'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_TA', nf_double, 'river_TA'_rank, 'r
     1iver_TA'_dims, 'river_TA'_id);
      call check_err(stat)

      'river_DIC_01'_dims(1) = 'river'_dim
      'river_DIC_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DIC_01', nf_double, 'river_DIC_01'_
     1rank, 'river_DIC_01'_dims, 'river_DIC_01'_id);
      call check_err(stat)

      'river_NO3_01'_dims(1) = 'river'_dim
      'river_NO3_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_NO3_01', nf_double, 'river_NO3_01'_
     1rank, 'river_NO3_01'_dims, 'river_NO3_01'_id);
      call check_err(stat)

      'river_NH4_01'_dims(1) = 'river'_dim
      'river_NH4_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_NH4_01', nf_double, 'river_NH4_01'_
     1rank, 'river_NH4_01'_dims, 'river_NH4_01'_id);
      call check_err(stat)

      'river_PO4_01'_dims(1) = 'river'_dim
      'river_PO4_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_PO4_01', nf_double, 'river_PO4_01'_
     1rank, 'river_PO4_01'_dims, 'river_PO4_01'_id);
      call check_err(stat)

      'river_DOC01_01'_dims(1) = 'river'_dim
      'river_DOC01_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DOC01_01', nf_double, 'river_DOC01_
     101'_rank, 'river_DOC01_01'_dims, 'river_DOC01_01'_id);
      call check_err(stat)

      'river_DOC02_01'_dims(1) = 'river'_dim
      'river_DOC02_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DOC02_01', nf_double, 'river_DOC02_
     101'_rank, 'river_DOC02_01'_dims, 'river_DOC02_01'_id);
      call check_err(stat)

      'river_DON01_01'_dims(1) = 'river'_dim
      'river_DON01_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DON01_01', nf_double, 'river_DON01_
     101'_rank, 'river_DON01_01'_dims, 'river_DON01_01'_id);
      call check_err(stat)

      'river_DON02_01'_dims(1) = 'river'_dim
      'river_DON02_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DON02_01', nf_double, 'river_DON02_
     101'_rank, 'river_DON02_01'_dims, 'river_DON02_01'_id);
      call check_err(stat)

      'river_DOP01_01'_dims(1) = 'river'_dim
      'river_DOP01_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DOP01_01', nf_double, 'river_DOP01_
     101'_rank, 'river_DOP01_01'_dims, 'river_DOP01_01'_id);
      call check_err(stat)

      'river_DOP02_01'_dims(1) = 'river'_dim
      'river_DOP02_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_DOP02_01', nf_double, 'river_DOP02_
     101'_rank, 'river_DOP02_01'_dims, 'river_DOP02_01'_id);
      call check_err(stat)

      'river_POC01_01'_dims(1) = 'river'_dim
      'river_POC01_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_POC01_01', nf_double, 'river_POC01_
     101'_rank, 'river_POC01_01'_dims, 'river_POC01_01'_id);
      call check_err(stat)

      'river_POC02_01'_dims(1) = 'river'_dim
      'river_POC02_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_POC02_01', nf_double, 'river_POC02_
     101'_rank, 'river_POC02_01'_dims, 'river_POC02_01'_id);
      call check_err(stat)

      'river_POC03_01'_dims(1) = 'river'_dim
      'river_POC03_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_POC03_01', nf_double, 'river_POC03_
     101'_rank, 'river_POC03_01'_dims, 'river_POC03_01'_id);
      call check_err(stat)

      'river_PON01_01'_dims(1) = 'river'_dim
      'river_PON01_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_PON01_01', nf_double, 'river_PON01_
     101'_rank, 'river_PON01_01'_dims, 'river_PON01_01'_id);
      call check_err(stat)

      'river_PON02_01'_dims(1) = 'river'_dim
      'river_PON02_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_PON02_01', nf_double, 'river_PON02_
     101'_rank, 'river_PON02_01'_dims, 'river_PON02_01'_id);
      call check_err(stat)

      'river_PON03_01'_dims(1) = 'river'_dim
      'river_PON03_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_PON03_01', nf_double, 'river_PON03_
     101'_rank, 'river_PON03_01'_dims, 'river_PON03_01'_id);
      call check_err(stat)

      'river_POP01_01'_dims(1) = 'river'_dim
      'river_POP01_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_POP01_01', nf_double, 'river_POP01_
     101'_rank, 'river_POP01_01'_dims, 'river_POP01_01'_id);
      call check_err(stat)

      'river_POP02_01'_dims(1) = 'river'_dim
      'river_POP02_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_POP02_01', nf_double, 'river_POP02_
     101'_rank, 'river_POP02_01'_dims, 'river_POP02_01'_id);
      call check_err(stat)

      'river_POP03_01'_dims(1) = 'river'_dim
      'river_POP03_01'_dims(2) = 'river_time'_dim
      stat = nf_def_var(ncid, 'river_POP03_01', nf_double, 'river_POP03_
     101'_rank, 'river_POP03_01'_dims, 'river_POP03_01'_id);
      call check_err(stat)
* assign global attributes
* define type
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'type', 17, 'ROMS FORCING 
     1file')
      call check_err(stat)
* define title
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'title', 22, 'TokyoBay Riv
     1er Forcing')
      call check_err(stat)
* define grd_file
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'grd_file', 22, 'TokyoBay2
     1_grid_v3.1.nc')
      call check_err(stat)
* define rivers
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'rivers', 165, '(1)Edogawa
     1 river, (2)Aragawa river (3)Tamagawa river, (4)Nakagawa river,(5)S
     1umidagawa river, (6)Tsurumigawa river, (7)Obitsugawa river, (8)Yor
     1o river, (9)Koito river')
      call check_err(stat)

* assign per-variable attributes
* define long_name
      stat = nf_put_att_text(ncid, 'river'_id, 'long_name', 34, 'river r
     1unoff identification number')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_Xposition'_id, 'long_name', 17
     1, 'river XI-position')
      call check_err(stat)
* define LuvSrc_meaning
      stat = nf_put_att_text(ncid, 'river_Xposition'_id, 'LuvSrc_meaning
     1', 40, 'i point index of U or V face source/sink')
      call check_err(stat)
* define LwSrc_meaning
      stat = nf_put_att_text(ncid, 'river_Xposition'_id, 'LwSrc_meaning'
     1, 39, 'i point index of RHO center source/sink')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_Eposition'_id, 'long_name', 18
     1, 'river ETA-position')
      call check_err(stat)
* define LuvSrc_meaning
      stat = nf_put_att_text(ncid, 'river_Eposition'_id, 'LuvSrc_meaning
     1', 40, 'j point index of U or V face source/sink')
      call check_err(stat)
* define LwSrc_meaning
      stat = nf_put_att_text(ncid, 'river_Eposition'_id, 'LwSrc_meaning'
     1, 39, 'j point index of RHO center source/sink')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_direction'_id, 'long_name', 22
     1, 'river runoff direction')
      call check_err(stat)
* define flag_values
      stat = nf_put_att_text(ncid, 'river_direction'_id, 'flag_values', 
     14, '0, 1')
      call check_err(stat)
* define flag_meanings
      stat = nf_put_att_text(ncid, 'river_direction'_id, 'flag_meanings'
     1, 38, 'flow across u-face, flow across v-face')
      call check_err(stat)
* define LwSrc_True
      stat = nf_put_att_text(ncid, 'river_direction'_id, 'LwSrc_True', 1
     13, 'flag not used')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_Vshape'_id, 'long_name', 44, '
     1river runoff mass transport vertical profile')
      call check_err(stat)
* define requires
      stat = nf_put_att_text(ncid, 'river_Vshape'_id, 'requires', 24, 'm
     1ust sum to 1 over s_rho')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_time'_id, 'long_name', 17, 'ri
     1ver runoff time')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_time'_id, 'units', 30, 'days s
     1ince 2000-01-01 00:00:00')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_transport'_id, 'long_name', 49
     1, 'river runoff vertically integrated mass transport')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_transport'_id, 'units', 15, 'm
     1eter3 second-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_transport'_id, 'time', 10, 'ri
     1ver_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_salt'_id, 'long_name', 21, 'ri
     1ver runoff salinity')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_salt'_id, 'time', 10, 'river_t
     1ime')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_temp'_id, 'long_name', 34, 'ri
     1ver runoff potential temperature')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_temp'_id, 'units', 7, 'Celsius
     1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_temp'_id, 'time', 10, 'river_t
     1ime')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DO'_id, 'long_name', 15, 'rive
     1r runoff DO')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DO'_id, 'units', 8, 'umol L-1'
     1)
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DO'_id, 'time', 10, 'river_tim
     1e')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_TA'_id, 'long_name', 23, 'rive
     1r runoff alkalinity')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_TA'_id, 'units', 9, 'umol kg-1
     1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_TA'_id, 'time', 10, 'river_tim
     1e')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DIC_01'_id, 'long_name', 19, '
     1river runoff DIC_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DIC_01'_id, 'units', 9, 'umol 
     1kg-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DIC_01'_id, 'time', 10, 'river
     1_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_NO3_01'_id, 'long_name', 19, '
     1river runoff NO3_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_NO3_01'_id, 'units', 8, 'umol 
     1L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_NO3_01'_id, 'time', 10, 'river
     1_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_NH4_01'_id, 'long_name', 19, '
     1river runoff NH4_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_NH4_01'_id, 'units', 8, 'umol 
     1L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_NH4_01'_id, 'time', 10, 'river
     1_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_PO4_01'_id, 'long_name', 19, '
     1river runoff PO4_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_PO4_01'_id, 'units', 8, 'umol 
     1L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_PO4_01'_id, 'time', 10, 'river
     1_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DOC01_01'_id, 'long_name', 21,
     1 'river runoff DOC01_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DOC01_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DOC01_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DOC02_01'_id, 'long_name', 21,
     1 'river runoff DOC02_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DOC02_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DOC02_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DON01_01'_id, 'long_name', 21,
     1 'river runoff DON01_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DON01_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DON01_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DON02_01'_id, 'long_name', 21,
     1 'river runoff DON02_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DON02_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DON02_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DOP01_01'_id, 'long_name', 21,
     1 'river runoff DOP01_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DOP01_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DOP01_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_DOP02_01'_id, 'long_name', 21,
     1 'river runoff DOP02_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_DOP02_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_DOP02_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_POC01_01'_id, 'long_name', 21,
     1 'river runoff POC01_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_POC01_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_POC01_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_POC02_01'_id, 'long_name', 21,
     1 'river runoff POC02_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_POC02_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_POC02_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_POC03_01'_id, 'long_name', 21,
     1 'river runoff POC03_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_POC03_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_POC03_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_PON01_01'_id, 'long_name', 21,
     1 'river runoff PON01_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_PON01_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_PON01_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_PON02_01'_id, 'long_name', 21,
     1 'river runoff PON02_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_PON02_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_PON02_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_PON03_01'_id, 'long_name', 21,
     1 'river runoff PON03_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_PON03_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_PON03_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_POP01_01'_id, 'long_name', 21,
     1 'river runoff POP01_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_POP01_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_POP01_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_POP02_01'_id, 'long_name', 21,
     1 'river runoff POP02_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_POP02_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_POP02_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)
* define long_name
      stat = nf_put_att_text(ncid, 'river_POP03_01'_id, 'long_name', 21,
     1 'river runoff POP03_01')
      call check_err(stat)
* define units
      stat = nf_put_att_text(ncid, 'river_POP03_01'_id, 'units', 8, 'umo
     1l L-1')
      call check_err(stat)
* define time
      stat = nf_put_att_text(ncid, 'river_POP03_01'_id, 'time', 10, 'riv
     1er_time')
      call check_err(stat)

* leave define mode
      stat = nf_enddef(ncid);
      call check_err(stat)

* assign scalar and fixed dimension variable data



* perform variable data writes
      stat = nf_close(ncid)
      call check_err(stat)
      end
      subroutine check_err(stat)
      integer stat
      include 'netcdf.inc'
      if (stat .ne. NF_NOERR) then
      print *, nf_strerror(stat)
      stop
      endif
      end
