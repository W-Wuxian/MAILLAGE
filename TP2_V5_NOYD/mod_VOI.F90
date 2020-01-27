module mod_VOI
  USE mod_parameters, ONLY : PR, STCO, STCN, NBPELE, NBPEDG, NV, NT, LB_head, UB_head, UB_voi, T1CO, T1CN, T1ED, T1ED_L
  IMPLICIT NONE

CONTAINS
  !BUILD_HEAD(Shead,ST1CN)
  !IN --> ST1CN must be connectivity table with size of NT*NBPELE
  !INOUT --> Shead size of NV+1 with NV.eqv. nbr vertex & Shead(0) = 0
  !         on stocke à l'indice i (>0) l'indice de FIN de la liste des
  !         voisins (ds ce cas des triangles ayant pour Sommet le Point i dans
  !         le second tableau Svoi).
  !La liste des voisins du point i débute dans le tableau Svoi à : Shead(i-1)+1
  !La liste des voisins du point i termine dans le tableau Svoi à : Shead(i)
  !NBR voisin (triangles ayant pour Sommet le pt i) du point i := Shead(i)-Shead(i-1)
  !INOUT --> Svoi size of (Shead(NV)), stocke de façon contigüe les # des triangles
  !ayant pour sommet le point i (i=1,NV).
  !#.eqv.numero
  subroutine BUILD_HEAD_VOI_CN(ST1CN,Shead,Svoi)
    implicit none
    INTEGER,      DIMENSION(:), intent(IN) :: ST1CN    !TABLE DE CONNECTIVITE
    integer, dimension(:), allocatable, intent(INOUT) :: Shead
    integer, dimension(:), allocatable, intent(INOUT) :: Svoi

    integer :: i, j, UB_ST1CN, hit, cpt

    !construction de Shead:
    UB_ST1CN = size(ST1CN)
    ALLOCATE(Shead(0:NV))
    LB_head = 0
    UB_head = NV
    hit = 0!; cpt = 0
    Shead(0) = 0
    hvertex: DO i = 1, NV
       htriangles:DO j = 1, UB_ST1CN
          IF(ST1CN(j) == i)THEN
             hit = hit + 1
          END IF
       END DO htriangles
       Shead(i) = hit
    END DO hvertex
    !construction de Svoi:
    UB_voi = hit
    ALLOCATE(Svoi(1:hit)) !because here hit == Shead(NV)
    cpt = 1
    vvertex: DO i = 1, NV
       hit = 0; j = 1
       vtriangles: DO WHILE( (j<UB_ST1CN).AND.(hit/=Shead(i)) )
          if(ST1CN(j)==i)then
             Svoi(cpt) = key1(j)
             cpt = cpt + 1; hit = hit + 1
          else if(ST1CN(j+1)==i)then
             Svoi(cpt) = key2(j+1)
             cpt = cpt + 1; hit = hit + 1
          else if (ST1CN(j+2)==i)then
             Svoi(cpt) = key3(j+2)
             cpt = cpt + 1; hit = hit + 1
          end if
          j = j + NBPELE
       END DO vtriangles
    END DO vvertex
  end subroutine BUILD_HEAD_VOI_CN

  !on note X le numero du triangle qui est formé du triplet de ses sommets (P1,P2,P3)
  !donc il existe une fct f de |N dans |N^3 tq:
  !                         f: X   --> (X1,X2,X3)
  !Avec X3 = X*NBPELE tq X2 = X3-1 et X1 = X3-2
  !On peut aussi ecrire f=(f1,f2,f3) tq f1(X)=X1, f2(X)=X2 et f3(X)=X3
  !Comme le maillage est conforme notamment au sens où chaque triangle est défini de
  !manière unique i.e unicité du triplet (P1,P2,P2) et que la fct f est lineraire (f3) et affine (f1 et f2),
  !elle est inversible.
  !On note: key1 l'inverse de f pour la composante P1 := X = (X1+2)/NBPELE
  !         key2 l'inverse de f pour la composante P2 := X = (X2+1)/NBPELE
  !         key3 l'inverse de f pour la composante P3 := X = (X3)/NBPELE
  !key --> input : a := #vertex
  !--> output : key := le # du triangle ayant pour sommet le vertex a
  integer function key1(a)
    implicit none
    integer, intent(IN) :: a
    key1 = (a+2)/NBPELE
  end function key1

  integer function key2(a)
    implicit none
    integer, intent(IN) :: a
    key2 = (a+1)/NBPELE
  end function key2

  integer function key3(a)
    implicit none
    integer, intent(IN) :: a
    key3 = a/NBPELE
  end function key3

  !************************************************************************************
  !INPUT  --> Shead and Svoi Come from computation of subroutine BUILD_HEAD_VOI_CN
  !OUTPUT --> Shead4CN and Svoi4CN
  Subroutine BUILD_CN_VOI_CN(ST1CN, Shead, Svoi,Shead4CN, Svoi4CN)
    Implicit None
    integer, dimension(:), intent(IN)  :: ST1CN    !TABLE DE CONNECTIVITE
    integer, dimension(0:NV), intent(IN)  :: Shead
    integer, dimension(:), intent(IN)  :: Svoi
    integer, dimension(:), allocatable, intent(INOUT)  :: Shead4CN
    integer, dimension(:), allocatable, intent(INOUT)  :: Svoi4CN

    integer :: P1, P2, P3
    integer :: i, j, j1, j2, j3
    integer :: UB_ST1CN
    integer :: P1d1, P1d2, P1nbrvoi
    integer :: P2d1, P2d2, P2nbrvoi
    integer :: P3d1, P3d2, P3nbrvoi

    integer, dimension(:), allocatable :: L1, L2, L3
    integer :: cpt12, cpt13, cpt23, sumcpt, ink, ctri
    integer, dimension(:), allocatable :: tmp_head, tmp_voi

    integer :: num1_voi1, num1_voi2, num1_voi3
    integer :: num2_voi1, num2_voi2, num2_voi3
    integer :: num3_voi1, num3_voi2, num3_voi3
    !AU MAX 3 VOISINS, AU MIN 1 VOISINS: (VOISIN .EQV. EDGE SHARED)
    ALLOCATE(tmp_voi(1:3*NT))
    ALLOCATE(tmp_head(0:NT))
    tmp_voi(:)  = 0
    tmp_head(:) = 0
    num1_voi1=0; num1_voi2=0; num1_voi3=0
    num2_voi1=0; num2_voi2=0; num2_voi3=0
    num3_voi1=0; num3_voi2=0; num3_voi3=0
    UB_ST1CN = size(ST1CN)
    ctri = 0
    TRIANGLE:DO i = 1, UB_ST1CN, NBPELE
       P1 = ST1CN(i); P2 = ST1CN(i+1); P3 = ST1CN(i+2)
       !POUR P1 du triangle i:
       P1d2 = Shead(P1); P1d1 = Shead(P1-1)
       P1nbrvoi = P1d2-P1d1
       ALLOCATE(L1(P1d1+1:P1d2))
       L1(P1d1+1:P1d2) = Svoi(P1d1+1:P1d2)
       WRITE(50,*)L1(P1d1+1:P1d2)
       !POUR P2 du triangle i:
       P2d2 = Shead(P2); P2d1 = Shead(P2-1)
       P2nbrvoi = P2d2-P2d1
       ALLOCATE(L2(P2d1+1:P2d2))
       L2(P2d1+1:P2d2) = Svoi(P2d1+1:P2d2)
       WRITE(50,*)L2(P2d1+1:P2d2)
       !POUR P3 du triangle i:
       P3d2 = Shead(P3); P3d1 = Shead(P3-1)
       P3nbrvoi = P3d2-P3d1
       ALLOCATE(L3(P3d1+1:P3d2))
       L3(P3d1+1:P3d2) = Svoi(P3d1+1:P3d2)
       !Regardons dans les listes de numero de triangle L1,L2etL3 si
       !Un autre numero que le numero du triangle i
       !Apparait Deux fois (=>) une edge commune donc un triangle voisin:
       !***VERIF LISTE L1 ET L2
       sumcpt = 0; ctri = ctri + 1 !ctri current numero du triangle  dans la boucle triangle
       num1_voi1 = 0; num2_voi1 = 0; num3_voi1 = 0
       L1L2:DO j1 = P1d1+1, P1d2
          if(L1(j1)/=ctri)then
             cpt12 = 0; j2 = P2d1+1
             DO WHILE(cpt12<1 .AND. j2<P2d2+1)
                if( (L2(j2)/=ctri).AND.(L2(j2)==L1(j1)) )then
                   cpt12 = cpt12 + 1
                   sumcpt = sumcpt + 1
                   select case(cpt12)
                   case(1)
                      num1_voi1 = L2(j2)
                   end select
                end if
                j2 = j2 +1
             END DO
          end if
       END DO L1L2
       !***VERIF LISTE L1 ET L3
       L1L3:DO j1 = P1d1+1, P1d2
          if(L1(j1)/=ctri)then
             cpt13 = 0; j3 = P3d1+1
             DO WHILE(cpt13<1 .AND. j3<P3d2+1)
                if( (L3(j3)/=ctri).AND.(L3(j3)==L1(j1)) )then
                   cpt13 = cpt13 + 1
                   sumcpt = sumcpt +1
                   select case(cpt13)
                   case(1)
                      num2_voi1 = L3(j3)
                   end select
                end if
                j3 = j3 +1
             END DO
          end if
       END DO L1L3
       !***VERIF LISTE L2 ET L3
       L2L3:DO j2 = P2d1+1, P2d2
          if(L2(j2)/=ctri)then
             cpt23 = 0; j3 = P3d1+1
             DO WHILE(cpt23<1.AND.j3<P3d2+1)
                if( (L3(j3)/=ctri).AND.(L3(j3)==L2(j2)) )then
                   cpt23 = cpt23 + 1
                   sumcpt = sumcpt +1
                   select case(cpt23)
                   case(1)
                      num3_voi1 = L3(j3)
                   end select
                end if
                j3 = j3 +1
             END DO
          end if
       END DO L2L3
       if(sumcpt<1 .OR. sumcpt>3)then
          DEALLOCATE(L1,L2,L3,tmp_head,tmp_voi)
          DEALLOCATE(T1CO, T1CN, T1ED, T1ED_L)
          stop 'sumcpt false, See Subroutine BUILD_CN_VOI_CN in mod_VOI.F90'
       end if
       tmp_head(ctri) = sumcpt
       tmp_voi(i)  = num1_voi1; tmp_voi(i+1)  = num2_voi1;tmp_voi(i+2)  = num3_voi1
       DEALLOCATE(L1,L2,L3)
    END DO TRIANGLE
    !DU COUPS, IL FAUT TRANSPOSER tmp_head dans Shead4CN tel que:
    !Shead4CN(i)-Shead4CN(i-1) = nbr triangle voisin au triangle i
    !ET TRANSPOSER Svoi4CN en "ommettant"/supprimant la ou il ya des zero
    ALLOCATE(Shead4CN(0:NT))
    Shead4CN(0) = 0
    DO i = 1, NT
       Shead4CN(i) = Shead4CN(i-1) + tmp_head(i)
    END DO
    ALLOCATE( Svoi4CN(Shead4CN(NT)) )
    ink = 0
    DO i = 1, 3*NT, 3
       DO j = 0,2
          if(tmp_voi(i+j)/=0)then
             ink = ink + 1
             Svoi4CN(ink) = tmp_voi(i+j)
          end if
       END DO
    END DO
  End Subroutine BUILD_CN_VOI_CN


end module mod_VOI
