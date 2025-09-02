********************************************************************************

TE-Japan �񋟃f�[�^�ɂ���


�F���q�󌤋��J���@�\�iJAXA�j
���F���Z�p����@�n���ϑ������Z���^�[�iEORC�j

2020�N03��06�� ����
2021�N11��16�� ����
2022�N05��23�� ����

********************************************************************************

�{�f�B���N�g������́A���{�̈��ΏۂƂ���MATSIRO�ACaMa-Flood�̃V�~�����[�V�������ʂ�
netCDF4�`���Ŏ擾�ł��܂��B

********************************************************************************
# �ڎ�
  1. Versions of Today�fs Earth 
  2. ���p�\�ȕϐ��ꗗ
  3. FTP�f�B���N�g���̍\��
  4. �t�@�C�������K��
  5. �t�H�[�}�b�g
  6. ��͒l�Ɨ\���l�̔��ʂɂ���
  7. �֘A����
  8. �Q�l����

********************************************************************************
# 1. Versions of TE-Japan

[TE-Japan]
    (Experiment)
        MSM/GPV ver.
        Satellite ver.

    (Spatial resolution)
        1/60-deg.

    (Temporal resolution)
        hourly
        daily
        monthly

    (Period)
        MSM/GPV ver. : 2007-present
        Satellite ver. : 2015-present

    (Latency)
        Realtime

    (Forcing)
        MSM/GPV ver. : Surface meteorological parameters by MSM/GPV
        Satellite ver. : Same as JRA-55 ver. except solar radiation from Himawari-8

********************************************************************************
# 2. ���p�\�ȕϐ��ꗗ

[Forcing]
  GPRCT     rainfall
  RPGPRCT   return period of rainfall*8
  GSNWL     snowfall
  GDU       wind speed
  GDT       surface air temperature
  GDQ       specific humidity
  SSRD      surface shortwave radiation (downward)
  SLRD      surface longwave radiation (downward)
  GDPS      surface air pressure

[MATSIRO]
 (Water balance (State))
  GLW       soil moisture (at each level) [Z1-Z6]*1
  GLWtot    soil moisture (total volume)
  GLWC      canopy water
  GLSNW     snow amount

 (Water balance (Flux))
  SNMLT     snow melt
  SNFRZ     snow freeze
  SNSUB     snow sublimation
  ICEMLT    ice melt
  ICESUB    ice sublimation
  SSUB      snow & ice sublimation
  ETFLX     transpiration
  EIFLX     canopy evaporation
  EISUB     canopy sublimation
  EBFLX     soil evaporation
  EBSUB     soil sublimation
  RUNOFF    total runoff (total) [W1-W2]*4
  RUNOFFB   base runoff
  SRUNOF    surface runoff
  RUNOFFA   runoff (lake & land) [W1-W2]*4

 (Heat balance (State))
  GLG       soil temperature [Z1-Z6]*1
  GLTSN     snow temperature [L1-L3]*2
  GLTS      land skin temperature [C1-C2]*3
  GLTC      canopy temperature [C1-C2]*3

 (Heat balance (Flux))
  GFLUXS    soil heat flux
  SNFLXS    snow surface heat flux
  GFLXTL    ground heat flux in total
  SSRU      surface shortwave radiation (upward)
  SLRU      surface longwave radiation (upward)
  SENS      sensible heat flux
  LTNT      latent heat flux
  EVAP      latent heat flux (evaporation)

 (River)
  RFLOW     river flow [W1-W2]*4
  GDRIV*    river water [W1-W2]*4
  GDRIVL    river storage [W1-W2]*4

 (Others)
  SNRAT     snow covered fraction
  ALB       albedo
  GLASN     snow albedo [A1-A3]*5
  GPSI      soil potential [Z1-Z6]*1
  CDSTM     dust density in snow [L1-L3]*2
  WA2L      water flux atmosphere to land
  WL2R      water flux land to river
  GLFRS     soil ice (at each level) [Z1-Z6]*1
  GLFRStot  soil ice (total volume)
  WLND      land water
  BUDIND    inland water sinkbudget
  RBUDIND   distributed water sinkbudget
  WINPT     ground water input
  SHLK      lake sh
  TSIL      lake surface temperature

[CaMa-Flood]
  RIVOUT    river discharge
  RIVSTO*   river water storage
  RIVDPH*   river water depth
  RPRIVDPH  return period of river water depth*8
  RIVVEL    river flow velocity
  FLDOUT    floodplain flow (discharge)
  FLDSTO*   floodplain water storage
  FLDDPH*   floodplain water depth
  FLDARE*   flood area
  FLDFRC*   flood fraction*9
  FLDFRC*   downscaled flood fraction*7*9
  SFCELV*   water surface elevation
  RPOUTFLW  return period of total discharge (RIVOUT + FLDOUT)*8


* : For hourly data, instantaneous values at the simulated time are stored.
    For daily & monthly data, averaged instantaneous values are stored.
    (Other variables are time averaged values for each temporal resolution.)

*1: Z1-Z6 represents the soil layers, the depth (m) of which is
    Z1: 0 - 0.05, Z2: 0.05 - 0.25, Z3: 0.25 - 1, Z4: 1 - 2, Z5: 2 - 4, and Z6: 4 - 14.
*2: L1-L3 represents the snow layers. The number of the effective layers
    and their depth are variable. See Takata et al. (2003) for more details.
*3: C1 and C2 represent the outputs for snow-free canopy and snow-covered canopy, respectively.
*4: W1 and W2 represent the values regarding water and ice, respectively.
*5: A1, A2 and A3 represent the snow albedo of visible, near-infrared and infrared area, respectively.
*6: River width is displayed with enhancement depending on its catchment size.
*7: Only provided in TE-Japan. 
    (Downscaled from 1/60 deg. to 1/240 deg. along with the unit catchment defined in the CaMa-Flood model.)
*8: Return period calculated from the statistical distribution of TE-Japan's simulation results 
    for the period 2007-2020. Only provided in TE-Japan.
*9: TE-Japan updated to FLDFRC with permanent water removed from 2020/09/02, 
    TE-Global updated to FLDFRC with permanent water removed from 2022/01/01.

********************************************************************************
# 3. FTP�f�B���N�g���̍\��

./TE-japan/
    ./MSM/         : Experiment name (MSM-GPV)
        ./monthly/ : Monthly average data (2007-present)
            ./yyyy
                ./mm

        ./daily/   : Daily average data
            ./yyyy
                ./mm
                    ./dd

        ./hourly/  : Houlry data
            ./yyyy
                ./mm
                    ./dd
                        ./hh

    ./SAT/         : Experiment name (Satellite)
       ...

 yyyy: �Z�o�Ώێ����i�^�C�����C���j�̔N�i4���j
 mm: �Z�o�Ώێ����i�^�C�����C���j�̌��i2���j
 dd: �Z�o�Ώێ����i�^�C�����C���j�̓��i2���j
 hh: �Z�o�Ώێ����i�^�C�����C���j�̎��i2���j

********************************************************************************
# 4. �t�@�C�������K��

# TE-japan:
  (monthly)
    TE-JPNrrr_eee_Myyyy01_xxxx.nc

  (daily)
    TE-JPNrrr_eee_Dyyyymm01_xxxx.nc

  (hourly)
    TE-JPNrrr_eee_Hyyyymmdd00_xxxx.nc


 rrr: ��ԕ���\��\��3����
      (e.g. '01M'�͈ܓx�o�x1���i�q��\��)
 eee: ��������\��3����
      (e.g. 'MSM'��MSM/GPV ver.��\��)
 yyyy: �Z�o�Ώێ����i�^�C�����C���j�̔N�i4���j
 mm: �Z�o�Ώێ����i�^�C�����C���j�̌��i2���j
 dd: �Z�o�Ώێ����i�^�C�����C���j�̓��i2���j
 xxxx: �ϐ���

 �t�@�C�����̗�: 
    TE-JPN01M_MSM_H2020030100_GPRCT.nc

********************************************************************************
# 5. �t�H�[�}�b�g
  �S�Ẵf�[�^��NetCDF4�`���ł���Agzip�ň��k���Ċi�[����Ă��܂��B

********************************************************************************
# 6. ��͒l�Ɨ\���l�̔��ʂɂ���(MSM/GPV ver.)

MSM/GPV ver.�̍~���ʂɂ��ẮC���ݎ����܂ł͋C�ے���͉J�ʁiRadar-AMeDAS�j���A
�\��ɂ��Ă͋C�ے����\���l�\�񃂃f��GPV�iMSM�j��39���ԗ\��l���@�B�w�K�ɂ��
�␳�������̂𗘗p���Ă��܂��i�ڍׂ̓v���_�N�g������2.1���Q�Ɓj�B

���҂����ʂ���ɂ́AnetCDF�t�@�C����global_attributes�ɋL�ڂ���Ă���"rain_forcing"
���m�F���Ă�������(�grain_forcing = amedas�h�ł���Ή�͉J�ʂ���͂������ʁA
"rain_forcing = msm"�ł���Ε␳��39���ԗ\��l����͂�������)�B
"rain_forcing"�̍��ڂ́A2022�N5��11���ȍ~�̃v���_�N�g�Ɋ܂܂�܂��B

��2022�N2��28���ȑO�̃f�[�^�́A���ׂċC�ے����\���l�\�񃂃f��GPV�iMSM�j�̉�͒l��p���Ă��܂��B

********************************************************************************
# 7. �֘A����

# TE-japan�T�v�ƌ��؂ɂ���:
  https://www.eorc.jaxa.jp/water/TE-japan/documents/TE-japan_20191129.pdf

********************************************************************************
# 8. �Q�l����
[Today's Earth]
  Yoshimura, K., T. Sakimura, T. Oki, S. Kanae, and S. Seto, Toward
  flood risk prediction: a statistical approach using a 29-year river
  discharge simulation over Japan, Hydrol. Res. Let., 2, 22-26, 2008.

  Ma, W., Ishitsuka, Y., Takeshima, A. et al. Applicability of a nationwide 
  flood forecasting system for Typhoon Hagibis 2019. Sci Rep 11, 10213 (2021). 
  https://doi.org/10.1038/s41598-021-89522-8

  Hibino, K., et al. Development of global hydrological monitoring and
  warning system: Today's Earth (TE) (in preparation)

  Yamamoto, K., et al. Evaluation of Hydrological Acceleration over 15 Years
  Analyzed from Global Terrestrial Hydrological Simulation (in preparation)


[MATSIRO]
  Takata, K, S. Emori, T. Watanabe, Development of the minimal advanced
  treatments of surface interaction and runoff, Global Planet.
  Change, 38, 209-222, 2003.

  Nitta, T, K. Yoshimura, K. Takata, R. O�fishi, T. Sueyoshi, S. Kanae,
  T. Oki, A. Abe-Ouchi, and G. E. Liston, Representing variability in
  subgrid snow cover and snow depth in a global land model: Offline
  validation, J.Clim., 27, 3318-3330, doi: 10.1175/jcli-d-13-00310.1,
  2014.

  Nitta, T., K. Yoshimura, A. Abe-Ouchi, Impact of arctic wetlands on
  the climate system: Model sensitivity simulations with the MIROC5 AGCM
  and a snow-fed wetland scheme, J. Hydrometeor.,
  doi:10.1175/JHM-D-16-0105.1, 2017.


[CaMa-Flood]
  Dai Yamazaki, Shinjiro Kanae, Hyungjun Kim, & Taikan Oki, A physically-based 
  description of floodplain inundation dynamics in a global river routing model,
  Water Resources Research, vol.54, W04501, 2011, DOI: 10.1029/2010WR009726

  Dai Yamazaki, Gustavo A. M. de Almeida, & Paul D. Bates, Improving computational
  efficiency in global river models by implementing the local inertial flow 
  equation and a vector-based river network map, Water Resources Research, 
  vol.49(11), pp.7221-7235, 2013, DOI:10.1002/wrcr.20552

  Dai Yamazaki, Tomoko Sato, Shinjiro Kanae, Yukiko Hirabayashi, & Paul D. Bates,
  Regional flood dynamics in a bifurcating mega delta simulated in a global 
  river model, Geophysical Research Letters, vol.41, pp.3127-3135, 2014, 
  DOI: 10.1002/2014GL059774

  Dai Yamazaki, Daiki Ikeshima, Ryunosuke Tawatari, Tomohiro Yamaguchi,
  Fiachra O'Loughlin, Jeff C. Neal, Christopher C. Sampson, Shinjiro
  Kanae, and Paul D. Bates, A high-accuracy map of global terrain elevations,
  Geophysical Research Letters, vol.44, pp.5844-5853, 2017, DOI:
  10.1002/2017GL072874

  Yamazaki, D., D. Ikeshima, J. Sosa, P.D. Bates, G.H. Allen, T.M. Pavelsky,
  MERIT Hydro: A high-resolution global hydrography map based on latest
  topography datasets, Water Resources Research, vol.55, pp.5053-5073, 2019, 
  DOI: 10.1029/2019WR024873

  �R���, �y�~���, �|����, ���R�h�m, ���{�S�捂�𑜓x�̕\�ʗ����f�[�^����,
  �y�؊w��_���WB1(���H�w), Vol.74(5), I_163-I_168, 2018

********************************************************************************
