program main
  USE mod_parameters
  USE mod_VERIF
  USE mod_VOI
  USE mod_CONNEC_COO
  Implicit None
  logical  :: test
  Real(PR),dimension(2) :: TXY
  Real(PR) :: WIZ
  Integer  :: i, j, UB
  Real(PR), Dimension(:), Allocatable :: CO3P

  !----------------------------------------------------------------------
  PRINT*,'DEFINE:', STCO,STCN,NBPELE,NBPEDG
  CALL RD_MF_CONNEC_COO(MESH_FILE,DIM,NV,NE,NT,T1CO,T1CN,T1ED,T1ED_L,T2D)

  CALL QK_INFO_CNCO(T1CN, T1CO, MESH_QK_CNCO)
  CALL BUILD_HEAD_VOI_CN(T1CN,head,voi)
  FCRTL = WR_verif_CNCO(MESH_QK_CNCO)
  DEALLOCATE(head, voi)

  CALL QK_INFO_T2DCO(T2D, T1CO, MESH_QK_T2DCO)
  CALL  BUILD_HEAD_VOI_T2D(T2D,head,voi)
  FCRTL = WR_verif_T2DCO(MESH_QK_T2DCO)
  DEALLOCATE(head, voi)

  TXY(1) = 0.1_PR; TXY(2) = 0.01_PR
  WIZ =0.000000000001_PR; ALLOCATE(CO3P(6))
  UB = ubound(T1CN,1)
  j = 0
  OPEN(30,file='geo.plt',position='append')
  WRITE(30,*)'set key'
  WRITE(30,*)'set size ratio -1'
  WRITE(30,*)'set style arrow 1 nohead lc rgb "forest-green"'
  WRITE(30,*)'set xrange [',-maxval(abs(T1CO))-1.0_PR,':',maxval(abs(T1CO))+1.0_PR,']'
  WRITE(30,*)'set yrange [',-maxval(abs(T1CO))-1.0_PR,':',maxval(abs(T1CO))+1.0_PR,']'
  CLOSE(30)
  PRINT*,'TARGET TXY',TXY(1),TXY(2),':'
  DO i = 1, UB, NBPELE
     j = j + 1
     CO3P(1) = T1CO(2*T1CN(i)-1);   CO3P(2) = T1CO(2*T1CN(i))
     CO3P(3) = T1CO(2*T1CN(i+1)-1); CO3P(4) = T1CO(2*T1CN(i+1))
     CO3P(5) = T1CO(2*T1CN(i+2)-1); CO3P(6) = T1CO(2*T1CN(i+2))
     test = in_c_inscr(TXY,CO3P,WIZ)
     SELECT CASE(test)
     CASE(.TRUE.)
        PRINT*,'EST DANS LE CERCLE CIRCONSCRIT AU TRIANGLE NUM=',j,'test',test
     CASE(.FALSE.)
        PRINT*,'N''EST DANS LE CERCLE CIRCONSCRIT AU TRIANGLE NUM=',j,'test',test
     END SELECT
  END DO
  OPEN(10,file='in_c_inscr.txt',position='append')
  WRITE(10,*)' '
  WRITE(10,*)' '
  WRITE(10,*)TXY(1),TXY(2)
  CLOSE(10)
  OPEN(30,file='geo.plt',position='append')
  WRITE(30,*)'plot "in_c_inscr.txt" index 0 u 1:2 ps 3 pt 2 title "pt tri" , "in_c_inscr.txt" index 1 u 1:2 ps 4 pt 3 title "target"  '
  WRITE(30,*)'replot "edg.txt" u 1:2:($3-$1):($4-$2) with vectors arrowstyle 1 title "edge" '
  WRITE(30,*)'replot "cercle.txt" with circles'
  WRITE(30,*)'replot "cercle.txt" u 1:2 ps 2 pt 4 title "center" '
  CLOSE(30)
  DEALLOCATE(CO3P)
  DEALLOCATE(T2D, T1CO, T1CN, T1ED, T1ED_L)

CONTAINS

  logical function in_c_inscr(Ptarget,CO3P,WIZ)
    implicit none
    Real(PR), Dimension(2), Intent(IN) :: Ptarget !X Y
    Real(PR), Dimension(6), Intent(IN) :: CO3P   !X1 Y1 X2 Y2 X3 Y3
    Real(PR), Intent(IN) :: WIZ ! What IS ZERO i.e la tol pour Dist
    Real(PR), Dimension(4) :: mid!, med!Point Milieu au segment, M�diatrice du segment
    Real(PR) :: P12X, P12Y, P23X, P23Y, P31X, P31Y
    Real(PR) :: Xc, Yc !CENTRE DU CERCLE CIRCONSCRIT
    Real(PR) :: Rayon, Dist

    Real(PR) :: a, aa, b, bb !a, aa:pente mediatrice et b, bb ordonnées à l'origine

    Integer  :: S1, S2, S3, S4
    Integer  :: k1, k2, k3, k4, k5, k6, k7, k8
    Integer  :: S5, S6

    Integer, save  :: cpt = 1


    !Calcul des composantes vectorielles du vecteur PiPj, (i,j=1,3):
    P12X = CO3P(3)-CO3P(1); P12Y = CO3P(4)-CO3P(2)
    P23X = CO3P(5)-CO3P(3); P23Y = CO3P(6)-CO3P(4)
    P31X = CO3P(1)-CO3P(5); P31Y = CO3P(2)-CO3P(6)
    !TEST POUR SAVOIR SI UN COTE EST PARALLELE A L'AXE Ox ET DONC EVITEE UNE DIVISION PAR 0:
    S5 = 0; S6 = 0
    IF(abs(P12Y) <= WIZ )THEN!P12 COLINEAIRE A L'AXE Ox
       ! DONC UTILISER P23 ET P31
       S1 = 2; k2 = 2*S1; k1 = k2-1
       S2 = 3; k4 = 2*S2; k3 = k4-1
       S3 = 3; k6 = 2*S3; k5 = k6-1
       S4 = 1; k8 = 2*S4; k7 = k8-1
       S5 = 10; S6 = 1
       !print*,'P12Y',P12Y
       !!$mid(1) = (CO3P(3) + CO3P(5))/2.0_PR; mid(2) = (CO3P(4) + CO3P(6))/2.0_PR !mid [P2P3].x, mid [P2P3].y
       !!$mid(3) = (CO3P(1) + CO3P(5))/2.0_PR; mid(4) = (CO3P(2) + CO3P(6))/2.0_PR !mid [P1P3].x, mid [P1P3].y
    ELSE IF(abs(P23Y) <= WIZ)THEN!P23 COLINEAIRE A L'AXE Ox
       ! DONC UTILISER P12 ET P31
       S1 = 1; k2 = 2*S1; k1 = k2-1
       S2 = 2; k4 = 2*S2; k3 = k4-1
       S3 = 3; k6 = 2*S3; k5 = k6-1
       S4 = 1; k8 = 2*S4; k7 = k8-1
       S5 = 10; S6 = 2
       !print*,'P23Y',P23Y
       !!$mid(1) = (CO3P(1) + CO3P(3))/2.0_PR; mid(2) = (CO3P(2) + CO3P(4))/2.0_PR !mid [P1P2].x, mid [P1P2].y
       !!$mid(3) = (CO3P(1) + CO3P(5))/2.0_PR; mid(4) = (CO3P(2) + CO3P(6))/2.0_PR !mid [P1P3].x, mid [P1P3].y
    ELSE IF(abs(P31Y) <= WIZ)THEN!P31 COLINEAIRE A L'AXE Ox
       ! DONC UTILISER P12 ET P23
       S1 = 1; k2 = 2*S1; k1 = k2-1
       S2 = 2; k4 = 2*S2; k3 = k4-1
       S3 = 2; k6 = 2*S3; k5 = k6-1
       S4 = 3; k8 = 2*S4; k7 = k8-1
       S5 = 10; S6 = 3
       !print*,'P31Y',P31Y
       !!$mid(1) = (CO3P(1) + CO3P(3))/2.0_PR; mid(2) = (CO3P(2) + CO3P(4))/2.0_PR !mid [P1P2].x, mid [P1P2].y
       !!$mid(3) = (CO3P(3) + CO3P(5))/2.0_PR; mid(4) = (CO3P(4) + CO3P(6))/2.0_PR !mid [P2P3].x, mid [P2P3].y
    END IF

    !Test Pour savoir si un second côtés est parallèle Oy:
    IF(abs(P12X) <= WIZ .OR. abs(P23Y) <= WIZ .OR. abs(P31Y) <= WIZ)THEN
      if(S5 == 10)then
         !TEST POUR SAVOIR SI UN COTE EST PARALLELE A L'AXE Oy:
         IF(abs(P12X) <= WIZ)THEN!P12 COLINEAIRE A L'AXE Oy
            ! DONC UTILISER P31 OU P23 POUR Xc (mean abscisses non egales des pts)
            Xc = (CO3P(1)+CO3P(5))*0.5_PR
            !Pour Yc comme on ne traite ques des points distincts,
            !Le vecteur PiPj est forcement non nul donc on peut
            !Directement prendre P12Y
            Yc = (CO3P(2)+CO3P(4))*0.5_PR !or (CO3P(2)+CO3P(6))*0.5_PR
            !SI hit :
            S5 = 15
         ELSE IF(abs(P23Y) <= WIZ)THEN!P23 COLINEAIRE A L'AXE Oy
            ! DONC UTILISER P12 OU P31 POUR Xc
            Xc = (CO3P(1)+CO3P(3))*0.5_PR
            !Pour Yc comme on ne traite ques des points distincts,
            !Le vecteur PiPj est forcement non nul donc on peut
            !Directement prendre P23Y
            Yc = (CO3P(4)+CO3P(6))*0.5_PR !or (CO3P(2)+CO3P(4))*0.5_PR
            !SI hit :
            S5 = 15
         ELSE IF(abs(P31Y) <= WIZ)THEN!P31 COLINEAIRE A L'AXE Oy
            ! DONC UTILISER P12 OU P23 POUR Xc
            Xc = (CO3P(1)+CO3P(3))*0.5_PR
            !Pour Yc comme on ne traite ques des points distincts,
            !Le vecteur PiPj est forcement non nul donc on peut
            !Directement prendre P31Y
            Yc = (CO3P(6)+CO3P(2))*0.5_PR !or (CO3P(2)+CO3P(4))*0.5_PR
            !SI hit :
            S5 = 15
         END IF
      end if
  END IF
    !ET DONC DANS LE CAS S5 == 10:
  IF(S5 == 10)THEN
    SELECT CASE(S6)
    CASE(1)
      mid(1) = (CO3P(k1) + CO3P(k3))/2.0_PR; mid(2) = (CO3P(k2) + CO3P(k4))/2.0_PR !mid [P2P3].x, mid [P2P3].y
      mid(3) = (CO3P(k7) + CO3P(k5))/2.0_PR; mid(4) = (CO3P(k8) + CO3P(k6))/2.0_PR !mid [P1P3].x, mid [P1P3].y
      !print*,'S6',S6,'k1',k1,'k3',k3,'k2',k2,'k4',k4
      !print*,CO3P(k1), CO3P(k3),' ',CO3P(k2),CO3P(k4)
    CASE(2)
      mid(1) = (CO3P(k1) + CO3P(k3))/2.0_PR; mid(2) = (CO3P(k2) + CO3P(k4))/2.0_PR !mid [P1P2].x, mid [P1P2].y
      mid(3) = (CO3P(k7) + CO3P(k5))/2.0_PR; mid(4) = (CO3P(k8) + CO3P(k6))/2.0_PR !mid [P1P3].x, mid [P1P3].y
      !print*,'S6',S6
    CASE(3)
      mid(1) = (CO3P(k1) + CO3P(k3))/2.0_PR; mid(2) = (CO3P(k2) + CO3P(k4))/2.0_PR !mid [P1P2].x, mid [P1P2].y
      mid(3) = (CO3P(k7) + CO3P(k5))/2.0_PR; mid(4) = (CO3P(k8) + CO3P(k6))/2.0_PR !mid [P2P3].x, mid [P2P3].y
      !print*,'S6',S6
    END SELECT
 ELSE IF(S5 == 0)THEN !peut importe dans celui-ci
    S1 = 1; k2 = 2*S1; k1 = k2-1
    S2 = 2; k4 = 2*S2; k3 = k4-1
    S3 = 3; k6 = 2*S3; k5 = k6-1
    S4 = 1; k8 = 2*S4; k7 = k8-1
    mid(1) = (CO3P(1) + CO3P(3))/2.0_PR; mid(2) = (CO3P(2) + CO3P(4))/2.0_PR !mid [P1P2].x, mid [P1P2].y
    mid(3) = (CO3P(3) + CO3P(5))/2.0_PR; mid(4) = (CO3P(4) + CO3P(6))/2.0_PR !mid [P2P3].x, mid [P2P3].y
  END IF

  IF(S5/=15)THEN
    a = -(CO3P(k3)-CO3P(k1))/(CO3P(k4)-CO3P(k2))
    b = (CO3P(k3)*CO3P(k3)-CO3P(k1)*CO3P(k1)+CO3P(k4)*CO3P(k4)-CO3P(k2)*CO3P(k2))*0.50_PR
    b = b/(CO3P(k4)-CO3P(k2))
    !PRINT*,'a,b',a,b
    aa = -(CO3P(k7)-CO3P(k5))/(CO3P(k8)-CO3P(k6))
    bb = (CO3P(k7)*CO3P(k7)-CO3P(k5)*CO3P(k5)+CO3P(k8)*CO3P(k8)-CO3P(k6)*CO3P(k6))*0.50_PR
    bb = bb/(CO3P(k8)-CO3P(k6))
    !PRINT*,'aa,bb',aa,bb
    Xc = (b-bb)/(aa-a); Yc = a * Xc + b
  END IF

  Rayon = SQRT( (CO3P(k1)-Xc)**2 + (CO3P(k2)-Yc)**2 )
  !print*,"k1 k2",k1,k2
  Dist  = SQRT( (Ptarget(1)-Xc)**2 + (Ptarget(2)-Yc)**2 )
  IF(Dist<=Rayon)THEN!ABS(Dist-Rayon)<WIZ)THEN
    in_c_inscr = .TRUE.
  ELSE
    in_c_inscr = .FALSE.
  END IF

  OPEN(10,file='in_c_inscr.txt',position='append')
  !WRITE(10,*)CO3P(:),Xc,Yc,Dist,Rayon
  WRITE(10,*)CO3P(1),CO3P(2)
  WRITE(10,*)CO3P(3),CO3P(4)
  WRITE(10,*)CO3P(5),CO3P(6)
  CLOSE(10)
  OPEN(50,file='edg.txt',position='append')
  WRITE(50,*)CO3P(1),CO3P(2),CO3P(3),CO3P(4)
  WRITE(50,*)CO3P(3),CO3P(4),CO3P(5),CO3P(6)
  WRITE(50,*)CO3P(5),CO3P(6),CO3P(1),CO3P(2)
  CLOSE(50)
  OPEN(20,file='cercle.txt',position='append')
  WRITE(20,*) Xc,Yc,Rayon
  CLOSE(20)
  !OPEN(30,file='geo.plt',position='append')
  !WRITE(30,*)'set object ',cpt,' circle at ',Xc,',',Yc,' radius ',Rayon
  !CLOSE(30)
  cpt = cpt +1
  !Calcul des point milieux:
  !!$mid(1) = (CO3P(1) + CO3P(3))/2.0_PR; mid(2) = (CO3P(2) + CO3P(4))/2.0_PR !mid [P1P2].x, mid [P1P2].y
  !!$mid(3) = (CO3P(3) + CO3P(5))/2.0_PR; mid(4) = (CO3P(4) + CO3P(6))/2.0_PR !mid [P2P3].x, mid [P2P3].y
  !!$mid(5) = (CO3P(1) + CO3P(5))/2.0_PR; mid(6) = (CO3P(2) + CO3P(6))/2.0_PR !mid [P1P3].x, mid [P1P3].y
  end function in_c_inscr
end program main
