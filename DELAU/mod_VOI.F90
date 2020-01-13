module mod_VOI
  USE mod_parameters, ONLY : PR, STCO, STCN, NBPELE, NBPEDG, NV, NT, tab_2D
  IMPLICIT NONE

  CONTAINS
!BUILD_HEAD(Shead,ST1CN)
!IN --> ST1CN must be connectivity table with size of NT*NBPELE
!INOUT --> Shead size of NV+1 with NV.eqv.br vertex & Shead(0) = 0
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

!4 T2D
subroutine BUILD_HEAD_VOI_T2D(ST2D,Shead,Svoi)
  implicit none
  TYPE(tab_2D), DIMENSION(:), intent(IN) :: ST2D    !TABLE DE CONNECTIVITE
  integer, dimension(:), allocatable, intent(INOUT) :: Shead
  integer, dimension(:), allocatable, intent(INOUT) :: Svoi

  integer :: i, j, UB_ST2D, hit, cpt

  !construction de Shead:
  UB_ST2D = size(ST2D)
  PRINT*,'size of T2D',UB_ST2D
  ALLOCATE(Shead(0:NV))
  hit = 0!; cpt = 0
  Shead(0) = 0
  hvertex: DO i = 1, NV
    htriangles:DO j = 1, UB_ST2D
      IF(ST2D(j)%TCN(1) == i)THEN
        hit = hit + 1
      ELSE IF(ST2D(j)%TCN(2) == i)THEN
        hit = hit + 1
      ELSE IF(ST2D(j)%TCN(3) == i)THEN
        hit = hit + 1
      END IF
    END DO htriangles
    Shead(i) = hit!; PRINT*,'i,Shead',i,Shead(i)
  END DO hvertex
  !construction de Svoi:
  ALLOCATE(Svoi(1:hit)) !because here hit == Shead(NV)
  cpt = 0
  vvertex: DO i = 1, NV
    hit = 0; j = 1
    vtriangles: DO WHILE( (j<UB_ST2D+1).AND.(hit/=Shead(i)) )
      !PRINT*,'TCN123',ST2D(j)%TCN(1),ST2D(j)%TCN(2),ST2D(j)%TCN(3),'j',j,'i',i
      if(ST2D(j)%TCN(1)==i)then
        cpt = cpt + 1; hit = hit + 1
        Svoi(cpt) = j!; PRINT*,'i,SVOI,cpt,j',i,Svoi(cpt),cpt,j
      else if(ST2D(j)%TCN(2)==i)then
        cpt = cpt + 1; hit = hit + 1
        Svoi(cpt) = j!; PRINT*,'i,SVOI,cpt,j',i,Svoi(cpt),cpt,j
      else if(ST2D(j)%TCN(3)==i)then
        cpt = cpt + 1; hit = hit + 1
        Svoi(cpt) = j!; PRINT*,'i,SVOI,cpt,j',i,Svoi(cpt),cpt,j
      end if
      j = j + 1
    END DO vtriangles
  END DO vvertex
end subroutine BUILD_HEAD_VOI_T2D
end module mod_VOI
