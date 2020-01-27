module mod_SEARCH_IN
  USE mod_parameters, ONLY : PR, STCO, STCN, NBPELE, NBPEDG, NV, NT
  implicit none


CONTAINS
  !SEARCH IF A TARGET POINT IS IN CERCLE INSCRIT AU TRIANGLE
  !INPUT:--> Ptarget le point target
  !      --> LE TRIANGLE CO3P
  !      --> LA TOL POUR LE CRITERE WIZ
  !OUTPUT:--> LOGICAL: TRUE IF IN ; FALSE IF NO
  !OTHER OUTPUT: 3 FICHIERS: 1/in_c_inscr.txt save: X1,Y1
  !                                                 X2,Y2
  !                                                 X3,Y3
  !2/edg.txt save les edges du triangle: CO3P(1),CO3P(2),CO3P(3),CO3P(4)
  !                                      CO3P(3),CO3P(4),CO3P(5),CO3P(6)
  !                                      CO3P(5),CO3P(6),CO3P(1),CO3P(2)
  !3/cercle.txt save : Xc,Yc,Rayon
  logical function in_c_inscr(Ptarget,CO3P,WIZ)
    implicit none
    Real(PR), Dimension(2), Intent(IN) :: Ptarget !X Y
    Real(PR), Dimension(6), Intent(IN) :: CO3P   !X1 Y1 X2 Y2 X3 Y3
    Real(PR), Intent(IN) :: WIZ ! What IS ZERO i.e la tol pour Dist
    Real(PR), Dimension(4) :: mid!, med!Point Milieu au segment, Mediatrice du segment
    Real(PR) :: P12X, P12Y, P23X, P23Y, P31X, P31Y
    Real(PR) :: Xc, Yc !CENTRE DU CERCLE CIRCONSCRIT
    Real(PR) :: Rayon, Dist

    Real(PR) :: a, aa, b, bb !a, aa:pente mediatrice et b, bb ordonnees a l'origine

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

    !Test Pour savoir si un second c�t�s est parall�le Oy:
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

!$$$ZEEIAD: SEEMS TO BE THE ROBUST WAY
real(PR) function side(Ptarget,CO2P)
  implicit none
  Real(PR), Dimension(2), Intent(IN) :: Ptarget   !Ptarget.x, Ptarget.y
  Real(PR), Dimension(4), Intent(IN) :: CO2P   !XA YA XB YB
  side = (Ptarget(1)-CO2P(3))*(CO2P(2)-CO2P(4))-(CO2P(1)-CO2P(3))*(Ptarget(2)-CO2P(4))
end function side
! Méthode:  side >0 si p est à gauche de la droite P1-P2
!           side =0 si le point est surla droite
!           side <0 ... droite
logical function IsInsideTriangle(CO3P,Ptarget)
  implicit none
! Teste si point H est intérieur au Triangle ABC.
! Méthode: Test si H et C sont a gauche de AB, ou pas, mais en meme te
  Real(PR), Dimension(6), Intent(IN) :: CO3P   !X1 Y1 X2 Y2 X3 Y3
  Real(PR), Dimension(2), Intent(IN) :: Ptarget   !Ptarget.x, Ptarget.y
  Real(PR) :: d1, d2, d3
  logical  :: neg, pos
  Real(PR), Dimension(4) ::  D1CO2P, D2CO2P, D3CO2P

  D1CO2P(1:4) = CO3P(1:4)
  D2CO2P(1:4) = CO3P(3:6)
  D3CO2P(1:2) = CO3P(5:6); D3CO2P(3:4) = CO3P(1:2)

  d1 = side(Ptarget,D1CO2P)
  d2 = side(Ptarget,D2CO2P)
  d3 = side(Ptarget,D3CO2P)

  neg = (d1 < 0.0_PR).OR.(d2 < 0.0_PR).OR.(d3 < 0.0_PR)
  pos = (d1 > 0.0_PR).OR.(d2 > 0.0_PR).OR.(d3 > 0.0_PR)

  IsInsideTriangle = .NOT.(neg.AND.pos)
  !SI TRUE LE POINT EST DANS LE TRIANGLE
  !SINON, LE POINT N'EST PAS DANS LE TRIANGLE
end function IsInsideTriangle




! not working
!!!!Programme adapt� au prog de base pour connaitre si oui ou non un point
!!!!est dans le triangle gr�ce aux oordonn�es barycentriques w1 w2 w3 du point
!!!!A modifier pour qu'ils renvoient les signes des coordonn�es barycentriques
  logical function in_triangle(Ptarget,COP3)
    implicit none
    Real(PR), Dimension(2), Intent(IN) :: Ptarget !X Y
    Real(PR), Dimension(6), Intent(IN) :: COP3   !X1 Y1 X2 Y2 X3 Y3
    !Real(PR), Intent(IN) :: WIZ ! What IS ZERO i.e la tol pour Dist
    Real(PR) :: Peri
    Real(PR) :: N1, N2, N3, SN
    Real(PR) :: E1, E2, E3, InfE
    Real(PR), Dimension(2) ::  AB, BC, AC
    Real(PR), Dimension(2) ::  MA, MB, MC
    AB(1) = COP3(3)-COP3(1); AB(2) = COP3(4)-COP3(2)
    BC(1) = COP3(5)-COP3(3); BC(2) = COP3(6)-COP3(4)
    AC(1) = COP3(5)-COP3(1); AC(2) = COP3(6)-COP3(2)
    E1 = SQRT(DOT_PRODUCT(AB,AB))
    E2 = SQRT(DOT_PRODUCT(BC,BC))
    E3 = SQRT(DOT_PRODUCT(AC,AC))
    InfE = min(min(E1,E2),E3)
    Peri = E1 + E2 + E3
    MA(1) = COP3(1)-Ptarget(1); MA(2) = COP3(2)-Ptarget(2)
    MB(1) = COP3(3)-Ptarget(1); MB(2) = COP3(4)-Ptarget(2)
    MC(1) = COP3(5)-Ptarget(1); MC(2) = COP3(6)-Ptarget(2)
    N1 = SQRT(DOT_PRODUCT(MA,MA))
    N2 = SQRT(DOT_PRODUCT(MB,MB))
    N3 = SQRT(DOT_PRODUCT(MC,MC))
    SN = N1 + N2 +N3
    in_triangle = .FALSE.
    !PRINT*,'Peri/2,SN,Peri-InfE',0.5_PR*Peri,SN,Peri-InfE
    IF(0.5_PR*Peri-SN<=0.0_PR.AND.SN-Peri+InfE<=0.0_PR)THEN
      in_triangle = .TRUE.
    END IF
   !Calcul de l'aire du triangle
   !      2
   !     /|\
   !    / | \
   !   /  |  \
   !  /P12XP23\
   ! /  /    \  \
   !1_/___P13__\_3

end function in_triangle

end module mod_SEARCH_IN
