c   onrdp.inc
c   Dimensions for onrdp.f
      integer MAXZ,MAXT                  ! Max. vertical points, time steps
      parameter(MAXZ=2047,MAXT=200000)
      integer NICEMAX, NSEDMAX           ! Max. no. ice and sediment classes
      parameter(NICEMAX=2,NSEDMAX=2)
c   ...note that the format statements for printing ice and sediment data
c   are not very flexible
c
