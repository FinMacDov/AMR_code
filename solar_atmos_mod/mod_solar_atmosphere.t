
module mod_solar_atmosphere

  use mod_global_parameters, only: std_len
  use mod_physics
  implicit none  
  !> Helium abundance over Hydrogen
  double precision, private    :: He_abundance

  !> resolution of profiles in interpolated tables
  integer, private :: npts

  !> Name of temperture profile chosen 
  character(len=std_len), private  :: Te_profile

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density
  integer, private, protected              :: e_

  !> The adiabatic index
  double precision, private :: sol_gamma

  double precision, allocatable :: z_pos(:), Te_intrpl(:), dLdz_pos(:)
  double precision, allocatable :: Yc(:), invYc(:)
  double precision  :: Tref, Lref, zmin,zmax
  double precision  :: lg_z_pos_min, lg_z_pos_max, lgstep
  integer          :: n_VALCIII, n_McWhirter, n_C7
  double precision :: rho_VALCIII(1:51), n_e_VALCIII(1:51), T_VALCIII(1:51), &
                      z_VALCIII(1:51), n_e_McWhirter(1:29), T_McWhirter(1:29), &
                      p_McWhirter(1:29), z_McWhirter(1:29), z_C7(1:140), & 
                      T_C7(1:140), n_e_C7(1:140), p_C7(1:140)
  double precision :: km_to_cgs 

  data    n_VALCIII / 51 / !< Number of data points
  
  !> gas density [g cm-3]
  data    rho_VALCIII / 3.192d-07, 3.08d-07, 2.949d-07, 2.727d-07, 2.152d-07, & 
                        1.606d-07, 1.15d-07, 5.413d-08, 2.334d-08, 9.327d-09, & 
                        4.902d-09, 3.232d-09, 1.899d-09, 1.121d-09, 6.846d-10, & 
                        4.358d-10, 1.902d-10, 1.297d-10, 7.359d-11, 4d-11, & 
                        1.839d-11, 9.822d-12, 5.315d-12, 2.45d-12, 1.493d-12, & 
                        6.082d-13, 3.227d-13, 2.417d-13, 2.122d-13, 1.802d-13, & 
                        1.628d-13, 1.53d-13, 1.433d-13, 1.225d-13, 1.093d-13, & 
                        9.569d-14, 6.39d-14, 5.619d-14, 5.216d-14, 5.058d-14, & 
                        4.795d-14, 4.517d-14, 4.355d-14, 4.203d-14, 4.017d-14, & 
                        3.665d-14, 3.222d-14, 2.808d-14, 2.113d-14, 1.179d-14, &
                        7.494d-15 /

  !> electron number denisty [cm-3]  
  data    n_e_VALCIII / 1.20400d15, 4.64500d14, 1.54700d14, 6.43300d13, 2.12200d13, & 
                        1.06600d13, 2.12200d13, 6.47600d12, 2.67400d12, 1.11000d12, & 
                        4.51600d11, 1.73300d11, 1.11200d11, 8.08500d10, 7.66400d10, & 
                        8.83800d10, 1.06400d11, 1.04900d11, 1.04100d11, 9.34900d10, & 
                        8.10800d10, 7.48600d10, 7.60000d10, 6.45600d10, 6.00500d10, & 
                        4.77100d10, 4.02800d10, 3.85800d10, 3.81100d10, 3.79200d10, & 
                        3.78300d10, 3.78000d10, 3.79900d10, 3.70500d10, 3.53500d10, & 
                        3.30600d10, 2.62000d10, 2.40200d10, 2.27600d10, 2.21900d10, & 
                        2.12000d10, 2.00900d10, 1.94300d10, 1.88100d10, 1.81200d10, & 
                        1.67700d10, 1.49800d10, 1.31800d10, 9.99300d10, 5.96100d09, &
                        3.83900d9 / 
  !> Temperature [K]
  data  T_VALCIII     / 8320, 7610, 6910, 6420, 5840, & 
                        5455, 5180, 4780, 4465, 4220, & 
                        4170, 4230, 4420, 4730, 5030, & 
                        5280, 5650, 5755, 5925, 6040, & 
                        6150, 6220, 6280, 6370, 6440, & 
                        6630, 6940, 7160, 7360, 7660, & 
                        7940, 8180, 8440, 9500, 10700, & 
                        12300, 18500, 21000, 22500, 23000, & 
                        23500, 24000, 24200, 24500, 25500, & 
                        28000, 32000, 37000, 50000, 89100, 141000 /
  
  !> height with respect to photosphere [cm]
  data  z_VALCIII      / -7500000.0, -5000000.0, -2500000.0, 0.0, 5000000.0, & 
                         10000000.0, 15000000.0, 25000000.0, 35000000.0, 45000000.0, & 
                         51500000.0, 55500000.0, 60500000.0, 65500000.0, 70500000.0, & 
                         75500000.0, 85500000.0, 90500000.0, 98000000.0, 106500000.0, & 
                         118000000.0, 128000000.0, 138000000.0, 151500000.0, 160500000.0, & 
                         178500000.0, 192500000.0, 199000000.0, 201600000.0, 205000000.0, & 
                         207000000.0, 208000000.0, 209000000.0, 210400000.0, 210700000.0, & 
                         210900000.0, 211300000.0, 211500000.0, 212000000.0, 212900000.0, & 
                         216000000.0, 220000000.0, 223000000.0, 225500000.0, 226300000.0, & 
                         226700000.0, 227100000.0, 227400000.0, 228000000.0, 229000000.0, &
                         229800000.0 /  

 data    n_McWhirter / 29 / !< Number of data points

!> electron number denisty [cm-3]  
 data    n_e_McWhirter / 4.13d6, 4.04d9, 2.11d10, 1.56d10, 1.05d10, 7.07d9, &
                         4.89d9, 3.24d9, 2.13d9, 1.37d9, 8.76d8, 6.38d8, &
                         4.82d8, 3.65d8, 2.74d8, 2.01d8, 1.43d8, 9.75d7, &
                         6.32d7, 3.89d7, 2.26d7, 1.26d7, 6.73d6, 3.56d6, &
                         1.92d6, 1.10d6, 6.88d5, 4.74d5, 3.57d5 /
  !> Temperature [K]
  data  T_McWhirter     / 8.94d3, 1.34d4, 2.01d4, 3.02d4, 4.53d4, 6.80d4, &
                          1.02d5, 1.53d5, 2.30d5, 3.45d5, 5.06d5, 6.47d5, &
                          7.89d5, 9.37d5, 1.10d6, 1.27d6, 1.45d6, 1.64d6, &
                          1.84d6, 2.02d6, 2.17d6, 2.28d6, 2.34d6, 2.35d6, &
                          2.36d6, 2.36d6, 2.36d6, 2.36d6, 2.36d6 /
  !> Hydrostatic Pressure [dyne cm-2]
  data  p_McWhirter     / 1.48d-1, 1.32d-1, 1.32d-1, 1.32d-1, 1.32d-1, 1.31d-1, &
                          1.30d-1, 1.29d-1, 1.27d-1, 1.23d-1, 1.15d-1, 1.07d-1, &
                          9.88d-2, 8.91d-2, 7.82d-2, 6.64d-2, 5.40d-2, 4.17d-2, &
                          3.02d-2, 2.04d-2, 1.28d-2, 7.45d-3, 4.09d-3, 2.18d-3, &
                          1.18d-3, 6.75d-4, 4.22d-4, 2.95d-4, 2.19d-4 /  
  !> height with respect to photosphere [cm]
  data  z_McWhirter     / 1.85d8, 1.88d8, 1.88d8, 1.88d8, 1.88d8, 1.89d8, 1.91d8, &
                          1.97d8, 2.13d8, 2.60d8, 3.90d8, 5.86d8, 8.79d8, 1.32d9, &
                          1.98d9, 2.97d9, 4.45d9, 6.68d9, 1.00d10, 1.51d10, 2.26d10, &
                          3.39d10, 5.08d10, 7.63d10, 1.14d11, 1.72d11, 2.58d11, 3.87d11, &
                          5.80d11 / 


  data    n_C7 / 140 / !< Number of data points

  data  z_C7            / -1.00000d7, -9.00000d6, -8.00000d6, -7.00000d6, -6.00000d6, & 
                          -5.00000d6, -4.00000d6, -3.00000d6, -2.00000d6, -1.00000d6, & 
                          0.00000d0, 1.00000d6, 2.00000d6, 3.50000d6, 5.00000d6, & 
                          7.50000d6, 1.00000d7, 1.25000d7, 1.50000d7, 1.75000d7, & 
                          2.00000d7, 2.50000d7, 3.00000d7, 3.50000d7, 4.00000d7, & 
                          4.50000d7, 4.90000d7, 5.25000d7, 5.60000d7, 6.15000d7, & 
                          6.60000d7, 7.00000d7, 7.50000d7, 8.00000d7, 8.54000d7, & 
                          9.00000d7, 9.46000d7, 9.71000d7, 1.00300d8, 1.03200d8, & 
                          1.06500d8, 1.10100d8, 1.14300d8, 1.21400d8, 1.29900d8, & 
                          1.39800d8, 1.52000d8, 1.61700d8, 1.72200d8, 1.82000d8, & 
                          1.89400d8, 1.94600d8, 1.98900d8, 2.02400d8, 2.05500d8, & 
                          2.08300d8, 2.09800d8, 2.11000d8, 2.12000d8, 2.12600d8, & 
                          2.13000d8, 2.13200d8, 2.13400d8, 2.13600d8, 2.13800d8, & 
                          2.14193d8, 2.14516d8, 2.14705d8, 2.14865d8, 2.15012d8, & 
                          2.15116d8, 2.15203d8, 2.15287d8, 2.15348d8, 2.15407d8, & 
                          2.15473d8, 2.15550d8, 2.15616d8, 2.15682d8, 2.15757d8, & 
                          2.15821d8, 2.15884d8, 2.15942d8, 2.16001d8, 2.16063d8, & 
                          2.16126d8, 2.16229d8, 2.16329d8, 2.16456d8, 2.16603d8, & 
                          2.16769d8, 2.16999d8, 2.17232d8, 2.17528d8, 2.17874d8, & 
                          2.18244d8, 2.18608d8, 2.19146d8, 2.19652d8, 2.20192d8, & 
                          2.20766d8, 2.21374d8, 2.22004d8, 2.22774d8, 2.23529d8, & 
                          2.24332d8, 2.25164d8, 2.25943d8, 2.26969d8, 2.28238d8, & 
                          2.29437d8, 2.30855d8, 2.32578d8, 2.34610d8, 2.37819d8, & 
                          2.41866d8, 2.45790d8, 2.49556d8, 2.53547d8, 2.57823d8, & 
                          2.62832d8, 2.68810d8, 2.76149d8, 2.84440d8, 2.97076d8, & 
                          3.13379d8, 3.42499d8, 3.75271d8, 4.29586d8, 4.96864d8, & 
                          5.76283d8, 7.36082d8, 8.97424d8, 1.15962d9, 1.53920d9, & 
                          2.11331d9, 2.66766d9, 3.60795d9, 4.70093d9, 6.80844d9 /
  !> Temperature [K]                         
  data  T_C7            /9.38000d3, 9.12000d3, 8.85000d3, 8.54000d3, 8.22000d3, & 
                          7.90000d3, 7.59000d3, 7.28000d3, 7.02000d3, 6.78000d3, & 
                          6.58300d3, 6.39700d3, 6.23100d3, 6.00600d3, 5.82600d3, & 
                          5.60700d3, 5.43100d3, 5.28800d3, 5.16500d3, 5.08000d3, & 
                          5.01000d3, 4.90700d3, 4.80500d3, 4.70000d3, 4.59000d3, & 
                          4.48500d3, 4.43500d3, 4.41000d3, 4.40000d3, 4.43500d3, & 
                          4.51000d3, 4.64000d3, 4.84000d3, 5.09000d3, 5.43000d3, & 
                          5.72000d3, 5.96900d3, 6.10000d3, 6.22500d3, 6.31500d3, & 
                          6.40000d3, 6.47400d3, 6.53100d3, 6.57600d3, 6.59800d3, & 
                          6.61000d3, 6.62300d3, 6.63300d3, 6.64300d3, 6.65200d3, & 
                          6.66000d3, 6.66700d3, 6.67400d3, 6.68000d3, 6.68600d3, & 
                          6.69400d3, 6.70000d3, 6.70600d3, 6.71800d3, 6.74000d3, & 
                          6.76800d3, 6.80000d3, 6.87000d3, 6.99200d3, 7.24800d3, & 
                          7.95000d3, 9.11500d3, 1.09800d4, 1.32000d4, 1.57600d4, & 
                          1.81400d4, 2.05100d4, 2.31000d4, 2.51200d4, 2.71300d4, & 
                          2.95000d4, 3.22600d4, 3.45800d4, 3.68700d4, 3.94000d4, & 
                          4.14500d4, 4.34000d4, 4.51400d4, 4.68000d4, 4.84900d4, & 
                          5.01400d4, 5.26900d4, 5.50200d4, 5.77900d4, 6.07900d4, & 
                          6.39500d4, 6.80000d4, 7.18100d4, 7.63300d4, 8.12200d4, & 
                          8.61200d4, 9.06400d4, 9.68600d4, 1.02300d5, 1.07800d5, & 
                          1.13200d5, 1.18700d5, 1.24000d5, 1.30200d5, 1.35800d5, & 
                          1.41600d5, 1.47200d5, 1.52300d5, 1.58600d5, 1.66000d5, & 
                          1.72600d5, 1.80100d5, 1.88600d5, 1.98100d5, 2.12000d5, & 
                          2.27800d5, 2.41900d5, 2.54400d5, 2.66700d5, 2.79000d5, & 
                          2.92500d5, 3.07300d5, 3.24000d5, 3.41300d5, 3.65000d5, & 
                          3.92000d5, 4.32900d5, 4.71600d5, 5.24300d5, 5.77400d5, & 
                          6.29300d5, 7.11700d5, 7.78300d5, 8.65000d5, 9.63700d5, & 
                          1.08000d6, 1.17000d6, 1.29400d6, 1.41000d6, 1.58600d6 /
!> electron number denisty [cm-3]  
  data  n_e_C7          / 3.80800d15, 2.92900d15, 2.19600d15, 1.54700d15, 1.05000d15, & 
                          6.91000d14, 4.46500d14, 2.79400d14, 1.82800d14, 1.20600d14, & 
                          8.39700d13, 5.89300d13, 4.25800d13, 2.72900d13, 1.90900d13, & 
                          1.23300d13, 8.74800d12, 6.59600d12, 5.13200d12, 4.10200d12, & 
                          3.30800d12, 2.18200d12, 1.43300d12, 9.33700d11, 6.03400d11, & 
                          3.87100d11, 2.71800d11, 1.99600d11, 1.46700d11, 9.16300d10, & 
                          6.42400d10, 5.15900d10, 4.76200d10, 5.45800d10, 8.05500d10, & 
                          1.16200d11, 1.51300d11, 1.72700d11, 1.90300d11, 2.02100d11, & 
                          2.09100d11, 2.10400d11, 1.99500d11, 1.76000d11, 1.48900d11, & 
                          1.43800d11, 1.42300d11, 1.26700d11, 1.02700d11, 7.96900d10, & 
                          6.45000d10, 5.52600d10, 4.82600d10, 4.30700d10, 3.88300d10, & 
                          3.49800d10, 3.29300d10, 3.12800d10, 2.98500d10, 2.89100d10, & 
                          2.82600d10, 2.78700d10, 2.74100d10, 2.68600d10, 2.60800d10, & 
                          2.43200d10, 2.32500d10, 2.27400d10, 2.04600d10, 1.76100d10, & 
                          1.55000d10, 1.38300d10, 1.23700d10, 1.14400d10, 1.06300d10, & 
                          9.81500d9, 9.01300d9, 8.43400d9, 7.93300d9, 7.45400d9, & 
                          7.10700d9, 6.80900d9, 6.56500d9, 6.34800d9, 6.14200d9, & 
                          5.95500d9, 5.68700d9, 5.46400d9, 5.22200d9, 4.98400d9, & 
                          4.75700d9, 4.49700d9, 4.28000d9, 4.04900d9, 3.82700d9, & 
                          3.62800d9, 3.46300d9, 3.25900d9, 3.09900d9, 2.95300d9, & 
                          2.82100d9, 2.70000d9, 2.59200d9, 2.47700d9, 2.37900d9, & 
                          2.28800d9, 2.20500d9, 2.13500d9, 2.05400d9, 1.96600d9, & 
                          1.89400d9, 1.81900d9, 1.73900d9, 1.65900d9, 1.55400d9, & 
                          1.44800d9, 1.36600d9, 1.30000d9, 1.24000d9, 1.18700d9, & 
                          1.13600d9, 1.08400d9, 1.03000d9, 9.79400d8, 9.16500d8, & 
                          8.52900d8, 7.69500d8, 7.02200d8, 6.24200d8, 5.57800d8, & 
                          5.02100d8, 4.27400d8, 3.76800d8, 3.20600d8, 2.66900d8, & 
                          2.14700d8, 1.80500d8, 1.41100d8, 1.10800d8, 7.49100d7 /
  !> gas Pressure [dyne cm-2]
  data  P_C7         / 1.99500d5, 1.90900d5, 1.82300d5, 1.73800d5, 1.65400d5, & 
                          1.57000d5, 1.48900d5, 1.41100d5, 1.33600d5, 1.26100d5, & 
                          1.18800d5, 1.11600d5, 1.04700d5, 9.48900d4, 8.57000d4, & 
                          7.18600d4, 5.97900d4, 4.94200d4, 4.06000d4, 3.31900d4, & 
                          2.70200d4, 1.77400d4, 1.15100d4, 7.38700d3, 4.69000d3, & 
                          2.94500d3, 2.01700d3, 1.44400d3, 1.03200d3, 6.10200d2, & 
                          3.99000d2, 2.75900d2, 1.76800d2, 1.15600d2, 7.50400d1, & 
                          5.30600d1, 3.81300d1, 3.20600d1, 2.58000d1, 2.12600d1, & 
                          1.71200d1, 1.35600d1, 1.03800d1, 6.64000d0, 3.92300d0, & 
                          2.15100d0, 1.06100d0, 6.37000d-01, 3.83000d-01, 2.48000d-01, & 
                          1.82000d-01, 1.48000d-01, 1.25000d-01, 1.09000d-01, 9.72000d-02, & 
                          8.76000d-02, 8.29000d-02, 7.93000d-02, 7.65000d-02, 7.49000d-02, & 
                          7.38000d-02, 7.33000d-02, 7.29000d-02, 7.25000d-02, 7.24000d-02, & 
                          7.24000d-02, 7.34000d-02, 7.57000d-02, 7.76000d-02, 7.81000d-02, & 
                          7.84000d-02, 7.87000d-02, 7.89000d-02, 7.91000d-02, 7.92000d-02, & 
                          7.92000d-02, 7.93000d-02, 7.93000d-02, 7.94000d-02, 7.95000d-02, & 
                          7.96000d-02, 7.97000d-02, 7.99000d-02, 8.00000d-02, 8.01000d-02, & 
                          8.03000d-02, 8.05000d-02, 8.07000d-02, 8.09000d-02, 8.11000d-02, & 
                          8.14000d-02, 8.17000d-02, 8.20000d-02, 8.23000d-02, 8.27000d-02, & 
                          8.30000d-02, 8.33000d-02, 8.37000d-02, 8.40000d-02, 8.43000d-02, & 
                          8.46000d-02, 8.48000d-02, 8.51000d-02, 8.53000d-02, 8.55000d-02, & 
                          8.57000d-02, 8.59000d-02, 8.60000d-02, 8.62000d-02, 8.63000d-02, & 
                          8.65000d-02, 8.66000d-02, 8.68000d-02, 8.69000d-02, 8.71000d-02, & 
                          8.73000d-02, 8.74000d-02, 8.74000d-02, 8.75000d-02, 8.76000d-02, & 
                          8.79000d-02, 8.81000d-02, 8.83000d-02, 8.84000d-02, 8.85000d-02, & 
                          8.84000d-02, 8.81000d-02, 8.76000d-02, 8.65000d-02, 8.52000d-02, & 
                          8.36000d-02, 8.04000d-02, 7.76000d-02, 7.33000d-02, 6.80000d-02, & 
                          6.13000d-02, 5.58000d-02, 4.83000d-02, 4.13000d-02, 3.14000d-02 / 
  contains

    !> Read this module"s parameters from a file
    subroutine sol_params_read(files)
      use mod_global_parameters, only: unitpar
      character(len=*), intent(in) :: files(:)
      integer                      :: n
  
      namelist /atmos_list/ Te_profile, npts
  
      do n = 1, size(files)
        open(unitpar, file=trim(files(n)), status="old")
        read(unitpar, atmos_list, end=111)
111     close(unitpar)
      end do

    end subroutine sol_params_read

    subroutine solar_atmosphere_init(phys_gamma,He_abund)
      use mod_global_parameters

      double precision, intent(in) :: phys_gamma,He_abund
      
      double precision, dimension(:), allocatable :: z_table
      double precision, dimension(:), allocatable :: Te_table
      double precision :: ratt, Lerror
      double precision ::fact1, fact2, fact3, dL1, dL2, z_shift
      double precision :: tstep, Lstep
      integer :: ntable, i, j, ic, nwx,idir, ix, clip1, clip2
      logical :: jump

      sol_gamma=phys_gamma
      He_abundance=He_abund
      Te_profile = 'T_VALCIII'

      call sol_params_read(par_files)

      ! Determine flux variables
      nwx = 1                  ! rho (density)

      allocate(mom(ndir))
      do idir = 1, ndir
         nwx    = nwx + 1
         mom(idir) = nwx       ! momentum density
      end do

      nwx = nwx + 1
      e_     = nwx          ! energy density

      allocate(z_pos(1:npts), Te_intrpl(1:npts), dLdz_pos(1:npts))
      allocate(Yc(1:npts), invYc(1:npts))
      
      z_pos(1:npts)    = zero
      Te_intrpl(1:npts)    = zero
      dLdz_pos(1:npts) = zero
      
      ! Read in the selected profiles
      select case(Te_profile)
      
      case('VALC3')
         if(mype ==0) &
         print *,'VALCIII table'

         ntable = n_VALCIII
      
         allocate(z_table(1:ntable))
         allocate(Te_table(1:ntable))
         z_table(1:ntable) = z_VALCIII(1:ntable) !< height
         Te_table(1:ntable) = T_VALCIII(1:ntable) !< temperature

      case('C7')
         if(mype ==0) &
         print *,'C7 table'      
         ntable = n_C7

         allocate(z_table(1:ntable))
         allocate(Te_table(1:ntable))
         z_table(1:ntable) = z_C7(1:ntable) !< height
         Te_table(1:ntable) = T_C7(1:ntable) !< temperature
      case('McWhirter')
         if(mype ==0) &
         print *,'McWhirter'      
         ntable = n_McWhirter

         allocate(z_table(1:ntable))
         allocate(Te_table(1:ntable))
         z_table(1:ntable) = z_McWhirter(1:ntable) !< height
         Te_table(1:ntable) = T_McWhirter(1:ntable) !< temperature
      !> This case needs to improved by manipulating data. Better to use C7. 
      case('VALMc')
         if(mype==0)&
         print*, 'VALCIII and McWhirter'
         clip1 = 15
         clip2 = 1
         z_shift = 2.2975e7
         ntable = n_VALCIII+n_McWhirter-clip1-clip2
         allocate(z_table(1:ntable))
         allocate(Te_table(1:ntable))
         z_table(1:n_VALCIII-clip1) = z_VALCIII(1:n_VALCIII-clip1) !< height
         z_table(n_VALCIII-clip1+1:ntable) = z_McWhirter(1+clip2:n_McWhirter)+z_shift !< height
         Te_table(1:n_VALCIII-clip1) = T_VALCIII(1:n_VALCIII-clip1) !< height
         Te_table(n_VALCIII-clip1+1:ntable) = T_McWhirter(1+clip2:n_McWhirter) !< height
      case default
         call mpistop("This T profile is unknown")
      end select
      
      ! create table(s) for use in MPI-AMRVAC
      zmax = z_table(ntable) !<physical location
      zmin = z_table(1)
      ratt = (zmax-zmin)/( dble(npts-1) + smalldouble) !<Step size
      !> starting points
      z_pos(1) = zmin 
      Te_intrpl(1) = Te_table(1)
      !> end points
      z_pos(npts) = zmax 
      Te_intrpl(npts) = Te_table(ntable)
      
      do i=2,npts        !< loop to create one table
        z_pos(i) = z_pos(i-1)+ratt
        do j=1,ntable-1   !< loop to create one spot on a table
        !> Second order polynomial interpolation, except at the outer edge, 
        !> or in case of a large jump.
          if(z_pos(i) < z_table(j+1)) then
             if(j.eq. ntable-1 )then
               fact1 = (z_pos(i)-z_table(j+1))     &
                     /(z_table(j)-z_table(j+1)) 
      
               fact2 = (z_pos(i)-z_table(j))       &
                     /(z_table(j+1)-z_table(j)) 
      
               Te_intrpl(i) = Te_table(j)*fact1 + Te_table(j+1)*fact2 
               exit
             else 
               dL1 = Te_table(j+1)-Te_table(j)
               dL2 = Te_table(j+2)-Te_table(j+1)
               jump =(max(dabs(dL1),dabs(dL2)) > 2*min(dabs(dL1),dabs(dL2)))
             endif
               
             if( jump ) then
               fact1 = (z_pos(i)-z_table(j+1))     &
                     /(z_table(j)-z_table(j+1)) 
      
               fact2 = (z_pos(i)-z_table(j))       &
                     /(z_table(j+1)-z_table(j)) 
                      
               Te_intrpl(i) = Te_table(j)*fact1 + Te_table(j+1)*fact2
               exit          
             else
               fact1 = ((z_pos(i)-z_table(j+1))     &
                     * (z_pos(i)-z_table(j+2)))   &
                     / ((z_table(j)-z_table(j+1)) &
                     * (z_table(j)-z_table(j+2)))
      
               fact2 = ((z_pos(i)-z_table(j))       &
                     * (z_pos(i)-z_table(j+2)))   &
                     / ((z_table(j+1)-z_table(j)) &
                     * (z_table(j+1)-z_table(j+2)))
      
               fact3 = ((z_pos(i)-z_table(j))       &
                     * (z_pos(i)-z_table(j+1)))   &
                     / ((z_table(j+2)-z_table(j)) &
                     * (z_table(j+2)-z_table(j+1)))
      
               Te_intrpl(i) = Te_table(j)*fact1 + Te_table(j+1)*fact2 &
                        + Te_table(j+2)*fact3
               exit
             endif
          endif
        enddo  ! end loop to find create one spot on a table
      enddo    ! end loop to create one table
      
      ! Scale both Z and Te
      z_pos(1:npts) = z_pos(1:npts)/unit_length
      Te_intrpl(1:npts) = Te_intrpl(1:npts)/unit_temperature  

    if (mype==0) then
      open(123,file='pruns',form='formatted')
!      write(123,*) npts
      do ix=1,npts
      write(123,"(i7,7es12.4)") ix, z_pos(ix)*unit_length, Te_intrpl(ix)*unit_temperature 
      enddo
      close(123)
    endif


!      deallocate(z_table)
!      deallocate(Te_table)           
    end subroutine solar_atmosphere_init

end module mod_solar_atmosphere
