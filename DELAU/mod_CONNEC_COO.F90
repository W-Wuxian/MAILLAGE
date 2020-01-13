#define D_STCO = 6
#define D_STCN = 3
!TRIANGLE T3
#define D_NBPELE = 3
!EDGE 2 PTS
#define D_NBPEDG = 2

module mod_CONNEC_COO
  USE mod_parameters, ONLY : NT, PR, STCO, STCN, NBPELE, NBPEDG, quality, tab_2D
  IMPLICIT NONE

CONTAINS
  !COMPUTE QUALITY ELEMENT:
  real(PR) function QK(CO)
    Implicit None
    !Integer,  Dimension(NBPELE), Intent(IN), optional :: TRI
    Real(PR), Dimension(STCO),   Intent(IN) :: CO
    !Integer  :: P1, P2, P3
    Real(PR) :: x1, y1, x2, y2, x3, y3
    Real(PR) :: l12x, l23x, l31x
    Real(PR) :: l12y, l23y, l31y
    Real(PR) :: surf, long
    !  P1 = TRI(1); P2 = TRI(2); P3 = TRI(3)
    x1 = CO(1); y1 = CO(2)
    x2 = CO(3); y2 = CO(4)
    x3 = CO(5); y3 = CO(6)
    l12x = x2 - x1; l12y = y2 - y1
    l23x = x3 - x2; l23y = y3 - y2
    l31x = x1 - x3; l31y = y1 - y3

    surf = (0.5_PR)*abs(( l12y*l31x - l12x*l31y ))
    long = l12x**2+l12y**2 + l23x**2+l23y**2 + l31x**2+l31y**2
    !print*,'surf',surf
    QK = (long/surf)*(sqrt(3.0_PR)/12)

  end function QK
  !COMPUTE QK OVER ALL ELEMENT
  subroutine QK_INFO_CNCO(T1CN, T1CO, info)
    Implicit None
    Integer,  Dimension(:), Intent(IN) :: T1CN
    Real(PR), Dimension(:), Intent(IN) :: T1CO
    type(quality), Intent(INOUT) :: info
    Integer  :: A, i
    Real(PR) :: tmp, switch_mi, switch_ma
    Real(PR), Dimension(:), Allocatable :: tmp_co
    A = ubound(T1CN,1)
    info%mo = 0.0_PR; info%ma = 0.0_PR; info%mi = 0.0_PR
    info%Qpc(:) = 0.0_PR
    switch_mi = 10000000.0_PR; switch_ma = -switch_mi
    Allocate(tmp_co(STCO))
    DO i = 1, A, NBPELE
      !tmp_cn(1) = TRI(i);tmp_cn(2) = TRI(i+1);tmp_cn(3) = TRI(i+2)
      tmp_co(1) = T1CO(2*T1CN(i)-1);   tmp_co(2) = T1CO(2*T1CN(i))
      tmp_co(3) = T1CO(2*T1CN(i+1)-1); tmp_co(4) = T1CO(2*T1CN(i+1))
      tmp_co(5) = T1CO(2*T1CN(i+2)-1); tmp_co(6) = T1CO(2*T1CN(i+2))
      tmp = QK(tmp_co)!QK(T1CN(i:i+2),tmp_co)
      info%mi = min(switch_mi, tmp)
      info%ma = max(switch_ma, tmp)
      switch_mi = tmp
      switch_ma = tmp
      info%mo = info%mo + tmp
      if(tmp>= 1.0_PR .AND. tmp<=2.0_PR )then
        info%Qpc(1) = info%Qpc(1) + 1
      else if ( tmp> 2.0_PR .AND. tmp<=3.0_PR ) then
        info%Qpc(2) = info%Qpc(2) + 1
      else if ( tmp> 3.0_PR .AND. tmp<=5.0_PR ) then
        info%Qpc(3) = info%Qpc(3) + 1
      else if ( tmp> 5.0_PR .AND. tmp<=10.0_PR ) then
        info%Qpc(4) = info%Qpc(4) + 1
      else if ( tmp> 10.0_PR .AND. tmp<=50.0_PR ) then
        info%Qpc(5) = info%Qpc(5) + 1
      else if ( tmp> 50.0_PR .AND. tmp<=100.0_PR ) then
        info%Qpc(6) = info%Qpc(6) + 1
      else if ( tmp> 100.0_PR ) then
        info%Qpc(7) = info%Qpc(7) + 1
      end if
    END DO
    info%mo = NBPELE * info%mo / A
    info%Qpc(:) = (info%Qpc(:) * 100 * NBPELE) / A
    Deallocate(tmp_co)
    PRINT*,'ubound(T1CN,1)',A,' NBPELE',NBPELE,' NT', A/NBPELE
  end subroutine QK_INFO_CNCO
  ! 4 T2D
  subroutine QK_INFO_T2DCO(TT2D, T1CO, info)
    Implicit None
    Type(tab_2D), Dimension(:), Intent(IN) :: TT2D
    Real(PR),     Dimension(:), Intent(IN) :: T1CO
    Type(quality), Intent(INOUT) :: info
    Integer  :: i
    Real(PR) :: tmp, switch_mi, switch_ma
    Real(PR), Dimension(:), Allocatable :: tmp_co

    info%mo = 0.0_PR; info%ma = 0.0_PR; info%mi = 0.0_PR
    info%Qpc(:) = 0.0_PR
    switch_mi = 10000000.0_PR; switch_ma = -switch_mi
    Allocate(tmp_co(STCO))
    DO i = 1, NT
      !tmp_cn(1) = TRI(i);tmp_cn(2) = TRI(i+1);tmp_cn(3) = TRI(i+2)
      tmp_co(1) = T1CO(2*TT2D(i)%TCN(1)-1);   tmp_co(2) = T1CO(2*TT2D(i)%TCN(1))
      tmp_co(3) = T1CO(2*TT2D(i)%TCN(2)-1); tmp_co(4) = T1CO(2*TT2D(i)%TCN(2))
      tmp_co(5) = T1CO(2*TT2D(i)%TCN(3)-1); tmp_co(6) = T1CO(2*TT2D(i)%TCN(3))
      tmp = QK(tmp_co)!QK(T1CN(i:i+2),tmp_co)
      info%mi = min(switch_mi, tmp)
      info%ma = max(switch_ma, tmp)
      switch_mi = tmp
      switch_ma = tmp
      info%mo = info%mo + tmp
      if(tmp>= 1.0_PR .AND. tmp<=2.0_PR )then
        info%Qpc(1) = info%Qpc(1) + 1
      else if ( tmp> 2.0_PR .AND. tmp<=3.0_PR ) then
        info%Qpc(2) = info%Qpc(2) + 1
      else if ( tmp> 3.0_PR .AND. tmp<=5.0_PR ) then
        info%Qpc(3) = info%Qpc(3) + 1
      else if ( tmp> 5.0_PR .AND. tmp<=10.0_PR ) then
        info%Qpc(4) = info%Qpc(4) + 1
      else if ( tmp> 10.0_PR .AND. tmp<=50.0_PR ) then
        info%Qpc(5) = info%Qpc(5) + 1
      else if ( tmp> 50.0_PR .AND. tmp<=100.0_PR ) then
        info%Qpc(6) = info%Qpc(6) + 1
      else if ( tmp> 100.0_PR ) then
        info%Qpc(7) = info%Qpc(7) + 1
      end if
    END DO
    info%mo = info%mo / NT
    info%Qpc(:) = (info%Qpc(:) * 100 ) / NT
    Deallocate(tmp_co)
  end subroutine QK_INFO_T2DCO


  !READ MESH FILE COMPUTE COO TAB & CONNECTIVITY TAB
  SUBROUTINE RD_MF_CONNEC_COO(MFILE,DIM,NV,NE,NT,T1CO,T1CN,T1ED,T1ED_L,T2D)
    Implicit None
    CHARACTER(LEN=*), INTENT(IN) :: MFILE           !mesh file name
    INTEGER, INTENT(INOUT)       :: DIM, NV, NE, NT          !nbr vertex, nbr triangles
    REAL(PR),     DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: T1CO    !TABLE DE COORDONNEES
    INTEGER,      DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: T1CN    !TABLE DE CONNECTIVITE
    INTEGER,      DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: T1ED    !TABLE EDGES
    INTEGER,      DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: T1ED_L    !TABLE EDGES AVEC LES LABELS
    TYPE(tab_2D), DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: T2D !IDEM+PHYSICAL TAG



    REAL(PR) :: X, Y, Z
    INTEGER  :: LABEL

    INTEGER          :: TAG
    CHARACTER(LEN=100) :: MFORMAT
    INTEGER          :: MDIM

    INTEGER :: i, i1, i2, e1, e3!, j
    INTEGER :: P1, P2, P3

    TAG = 10
    OPEN(TAG, file = MFILE, ACTION = "READ")
    READ(TAG,*)MFORMAT
    PRINT*,TRIM(ADJUSTL(MFORMAT))
    READ(TAG,*)MFORMAT
    PRINT*,TRIM(ADJUSTL(MFORMAT))
    READ(TAG,*)MDIM
    PRINT*,'DIM=', MDIM-1 !2D
    MDIM = MDIM -1; DIM = MDIM
    READ(TAG,*)MFORMAT
    PRINT*,TRIM(ADJUSTL(MFORMAT))
    READ(TAG,*)NV
    ALLOCATE( T1CO( (MDIM)*NV ) )
    PRINT*,'NBR VERTICES = ',NV
    i1 = 1
    vertex:DO i = 1, NV
       READ(TAG,*) X, Y, Z, LABEL
       PRINT*, X, Y, Z, LABEL
       i2 = i1 + 1
       T1CO(i1) = X
       T1CO(i2) = Y
       i1 = i2 + 1
    END DO vertex
    READ(TAG,*)MFORMAT
    PRINT*,MFORMAT
    READ(TAG,*)NE
    ALLOCATE( T1ED( NE*(NBPEDG) ), T1ED_L( NE*(NBPEDG+1) ) )
    PRINT*,'NBR EDGES =',NE
    i1 = 1; e1= 1
    edges: DO i=1, NE
    READ(TAG,*)P1, P2, LABEL
    PRINT*, P1, P2, LABEL
    i2 = i1 + 1
    T1ED(i1) = P1
    T1ED(i2) = P2
    i1 = i2 + 1
    e3 = e1 +2
    T1ED_L(e1)   = P1
    T1ED_L(e1+1) = P2
    T1ED_L(e3)   = LABEL
    e1 = e3 +1
    END DO edges
    READ(TAG,*)MFORMAT
    PRINT*,MFORMAT
    READ(TAG,*)NT
    ALLOCATE( T1CN( NT*(NBPELE) ) )
    PRINT*,'NBR TRIANGLES =',NT
    i1 = 1
    triangles: DO i=1, NT
    READ(TAG,*)P1, P2, P3, LABEL
    PRINT*, P1, P2, P3, LABEL
    i2 = i1 + 2
    T1CN(i1)   = P1
    T1CN(i1+1) = P2
    T1CN(i2)   = P3
    i1 = i2 + 1
  END DO triangles
  READ(TAG,*)MFORMAT
  PRINT*,MFORMAT
  CLOSE(TAG)

  !FORMATION DE TD:
  ALLOCATE(T2D(NT))
  i1 = 1
  DO i = 1, NT*NBPELE, NBPELE
    T2D(i1)%TCN(1) = T1CN(i)
    !T2D(i1)%TCO(1) = T1CO(2 * (T2D(i1)%TCN(1)) - 1); T2D(i1)%TCO(2) = T1CO(2 * (T2D(i1)%TCN(1)))
    T2D(i1)%TCN(2) = T1CN(i+1)
    !T2D(i1)%TCO(3) = T1CO(2 * (T2D(i1)%TCN(2)) - 1); T2D(i1)%TCO(4) = T1CO(2 * (T2D(i1)%TCN(2)))
    T2D(i1)%TCN(3) = T1CN(i+2)
    !T2D(i1)%TCO(5) = T1CO(2 * (T2D(i1)%TCN(3)) - 1); T2D(i1)%TCO(6) = T1CO(2 * (T2D(i1)%TCN(3)))
    i1 = i1 + 1
  END DO
  END SUBROUTINE RD_MF_CONNEC_COO

end module mod_CONNEC_COO
