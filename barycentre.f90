!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Programme adapté au prog de base pour connaitre si oui ou non un point est dans le triangle grâce aux oordonnées barycentriques w1 w2 w3 du point 
!!!!A modifier pour qu'ils renvoient les signes des coordonnées barycentriques



  logical function in_triangle(Ptarget,CO3P)
    implicit none
    Real(PR), Dimension(2), Intent(IN) :: Ptarget !X Y
    Real(PR), Dimension(6), Intent(IN) :: CO3P   !X1 Y1 X2 Y2 X3 Y3
    !Real(PR), Intent(IN) :: WIZ ! What IS ZERO i.e la tol pour Dist

    Real(PR) :: A123, P12,P23,P31
    Real(PR)::xP2_3,yP2_3,xP3_1,yP2_3,xP1_2,yP1_2,w1,w2,w3




    !Calcul de l'aire du triangle

    P12=COP3(1)*COP3(4)-COP3(2)*COP3(3)
    P23=COP3(3)*COP3(6)-COP3(4)*COP3(5)
    P31=COP3(5)*COP3(2)-COP3(6)*COP3(1)
    A123=abs((P12+P23+P31)/2)


    !Calcul des coordonnées des 'différences' des sommets

    xP2_3=COP3(3)-COP3(5);      yP2_3=COP3(4)-COP3(6);
    xP3_1=COP3(5)-COP3(1);      yP3_1=COP3(6)-COP3(2);
    xP1_2=COP3(1)-COP3(3);      yP1_2=COP3(2)-COP3(4);

    !Calcul des coordonnées barycentriques

    w1=(P23-Ptarget(1)*yP2_3+Ptarget(2)*xP2_3)/A123

    w2=(P31-Ptarget(1)*yP3_1+Ptarget(2)*xP3_1)/A123

    w3=(P12-Ptarget(1)*yP1_2+Ptarget(2)*xP1_2)/A123




  IF(w1>0 .and. w2>0 .and. w3>0)THEN
    in_triangle = .TRUE.
  ELSE
    in_triangle = .FALSE.
  END IF



end function in_triangle

