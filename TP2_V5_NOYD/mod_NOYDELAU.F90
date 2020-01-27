module mod_NOYDELAU
  USE mod_parameters, ONLY : PR, STCN, STCO, NBPEDG, NBPELE, SWIZ, NV, NT, T1CO, T1CN, TAR_GET, ISBASE, CAV, head, voi
  USE mod_VOI
  USE mod_SEARCH_IN
  implicit none

CONTAINS

  !INPUT: LTAR LA LISTE DES POINTS A INSERER
  !EX: DIM(LTAR(:))=6 => LTAR(1:2)<=>X1,Y1
  !                      LTAR(3:4)<=>X2,Y2
  !                      LTAR(5:6)<=>X3,Y3
  !OUTPUT: T1CO ET T1CN
  Integer Function NOYDELAU(LTAR)
    implicit none
    Type(TAR_GET), Dimension(:), Intent(IN) :: LTAR
    REAL(PR), DIMENSION(2) :: Ptarget
    REAL(PR), DIMENSION(6) :: CO3P
    TYPE(ISBASE) :: BASE
    Integer, Dimension(:), Allocatable :: head4CN, voi4CN
    Integer :: S_LTAR ! SIZE OF LTAR
    Integer :: UB_T1CN, UB_head4CN ! SIZE OF T1CN
    integer :: i, j, ctri
    logical :: IN_TRI!, IN_CIRC !in traingle, in cercle circonsrit
    NOYDELAU = 0
    S_LTAR = size(LTAR)!need to be nbr point a inserer

    BLTAR:DO i = 1, S_LTAR
      j = 1; ctri = 0
      UB_T1CN = ubound(T1CN,1)
      Ptarget(1) = LTAR(i)%XY(1);   Ptarget(2) = LTAR(i)%XY(2)
      CALL BUILD_HEAD_VOI_CN(T1CN,head,voi)
      CALL BUILD_CN_VOI_CN(T1CN, head, voi, head4CN, voi4CN)
      UB_head4CN = ubound(head4CN,1)
      !DEALLOCATE(head, voi)
      BTRI:DO WHILE(j<UB_T1CN)! last admissible j = UB_T1CN-NBPELE
        ctri = ctri + 1
        CO3P(1) = T1CO(2*T1CN(j)-1);   CO3P(2) = T1CO(2*T1CN(j))
        CO3P(3) = T1CO(2*T1CN(j+1)-1); CO3P(4) = T1CO(2*T1CN(j+1))
        CO3P(5) = T1CO(2*T1CN(j+2)-1); CO3P(6) = T1CO(2*T1CN(j+2))
        IN_TRI = IsInsideTriangle(CO3P,Ptarget)
        SELECT CASE(IN_TRI)
        CASE(.TRUE.)
          CALL COMP_BASE(BASE,ctri,CO3P)
          !IN_CIRC = in_c_inscr(Ptarget,CO3P,SWIZ)
          CALL search_mod_cav(Ptarget,ctri,j, UB_head4CN, head4CN, voi4CN,CAV)
          j = UB_T1CN
          !DEALLOCATE(head4CN, voi4CN)
        CASE(.FALSE.)
          j = j + NBPELE
        END SELECT
      END DO BTRI
      DEALLOCATE(head, voi)
      DEALLOCATE(head4CN, voi4CN)
    END DO BLTAR
    NOYDELAU = 1

  End Function NOYDELAU

  subroutine COMP_BASE(SBASE,curtri,CO3P)
    implicit none
    TYPE(ISBASE), INTENT(INOUT) :: SBASE
    INTEGER, INTENT(IN) :: curtri
    REAL(PR), DIMENSION(6), INTENT(IN) :: CO3P
    integer :: k, lu
    SBASE%NUPB = curtri
    lu = curtri*3
    SBASE%PBCN(1) = T1CN(lu-2); SBASE%PBCN(2) = T1CN(lu-1); SBASE%PBCN(2) = T1CN(lu)
    DO k = 1, STCO
      SBASE%PBCO(k) = CO3P(k)
    END DO
  end subroutine COMP_BASE
  !0/ we put -1 in T1CN for the BASE, because in this tri means also in cercle circonsrit
  !1/SEARCH is in circle circonsrit
  !2/mod T1CN BY PUTTING -1 for tri found  in 1/
  !3/save cavity
  subroutine search_mod_cav(TXY,curtri,ll, UB_Shead4CN, Shead4CN, Svoi4CN,SCAV)
    implicit none
    REAL(PR), DIMENSION(2) :: TXY
    !ll = j
    INTEGER, INTENT(IN)   :: curtri, ll, UB_Shead4CN
    INTEGER, DIMENSION(0:UB_Shead4CN), INTENT(IN) :: Shead4CN
    INTEGER, DIMENSION(:), INTENT(INOUT) :: Svoi4CN

    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: SCAV
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_SCAV, add_SCAV
    REAL(PR), DIMENSION(6) :: SCO3P
    logical :: IN_CIRC
    integer :: i, j, BnbVOI !BASE NBR VOISIN, j
    integer, dimension(:), allocatable :: BIS
    integer :: a, i1,i2!, alloc
    !-1 dans T1CN POUR LE TRIANGLE BASE:
    ALLOCATE(SCAV(3))
    SCAV(1:3) = T1CN(ll:ll+2)
    T1CN(ll:ll+2) = -1
    DO i=1,size(Svoi4CN)
      if(Svoi4CN(i)==curtri)then
        Svoi4CN(i)=-1
      end if
    END DO
    DO i=1,size(voi)
      if(voi(i)==curtri)then
        voi(i)=-1
      end if
    END DO
    !search
    BnbVOI = Shead4CN(curtri)-Shead4CN(curtri-1)
    ALLOCATE(BIS(3))!(BIS(Shead4CN(curtri-1)+1:Shead4CN(curtri)))
    BIS(1:3) = SCAV(1:3)!BIS = -2
    i1 = Shead4CN(curtri-1)+1; i2 = Shead4CN(curtri)
    DO i = i1, i2
      !if(i/=-1.AND.i/=-2)then
        !BIS(i)=Svoi4CN(i)
        PRINT*,'i',i,'voi',Svoi4CN(i)
        a = (NBPELE*Svoi4CN(i)) - 2 !i.e le triangle numero Svoi4CN(i) est T1CN(a:a+2)
        PRINT*,'noy a',a
        SCO3P(1) = T1CO(2*T1CN(a)-1);   SCO3P(2) = T1CO(2*T1CN(a))
        SCO3P(3) = T1CO(2*T1CN(a+1)-1); SCO3P(4) = T1CO(2*T1CN(a+1))
        SCO3P(5) = T1CO(2*T1CN(a+2)-1); SCO3P(6) = T1CO(2*T1CN(a+2))
        IN_CIRC = in_c_inscr(TXY,SCO3P,SWIZ)
        SELECT CASE(IN_CIRC)
        CASE(.TRUE.)
          ALLOCATE(tmp_SCAV(1:size(SCAV)))
          tmp_SCAV(:) = SCAV(:)
          DEALLOCATE(SCAV)
          ALLOCATE(SCAV(size(tmp_SCAV)+1))
          SCAV(1:size(tmp_SCAV))=tmp_SCAV(:)
          SCAV(size(tmp_SCAV)+1) = compa(tmp_SCAV(1:3),T1CN(a:a+2),a)
          PRINT*,'tmp_SCAV',tmp_SCAV(:)
          PRINT*,'t1CN',T1CN(a:a+2)
          DEALLOCATE(tmp_SCAV)
          T1CN(a:a+2) = -1
            DO j=1,size(Svoi4CN)
              if(Svoi4CN(j)==Svoi4CN(i))then
                Svoi4CN(j)=-1
              end if
            END DO
            DO j=1,size(voi)
              if(voi(j)==Svoi4CN(i))then
                voi(j)=-1
              end if
            END DO
        CASE(.FALSE.)
            DO j=1,size(Svoi4CN)
              if(Svoi4CN(j)==Svoi4CN(i))then
                Svoi4CN(j)=-2
              end if
            END DO
            DO j=1,size(voi)
              if(voi(j)==Svoi4CN(i))then
                voi(j)=-2
              end if
            END DO
        END SELECT
      !end if

      ! if(i==i2)then
      !   i1 = 1; i2 = size(BIS); h=h+1
      !   cycle
      ! end if
    END DO
    CALL tree(BIS, TXY, SCAV, add_SCAV, Svoi4CN)
    DEALLOCATE(SCAV)
    ALLOCATE(SCAV(1:size(add_SCAV)))
    SCAV(:)=add_SCAV(:)
    DEALLOCATE(add_SCAV)
    DEALLOCATE(BIS)
  end subroutine search_mod_cav


  subroutine tree(fBIS, STXY, CHAINEDSCAV, SCAV, Svoi4CN)
    implicit none
    integer, dimension(3), intent(IN) :: fBIS
    REAL(PR), DIMENSION(2), INTENT(IN) :: STXY
    INTEGER, DIMENSION(:), INTENT(IN) :: CHAINEDSCAV
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: SCAV
    !INTEGER, DIMENSION(0:UB_Shead4CN), INTENT(IN) :: Shead4CN
    INTEGER, DIMENSION(:), INTENT(INOUT) :: Svoi4CN
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_SCAV
    REAL(PR), DIMENSION(6) :: SSCO3P
    integer :: i, j, k, di, j1, j2, Sa
    logical :: SIN_CIRC
    ALLOCATE( SCAV(1:size(CHAINEDSCAV)) )
    SCAV(:)=CHAINEDSCAV(:)
    di = 3
    i = 1
    DO WHILE(i<=di)
      j1 = head(fBIS(i)-1)+1; j2 = head(fBIS(i))
      DO j =   j1, j2
        if(voi(j)/=-1.AND.voi(j)/=-2)then
          Sa = (NBPELE*voi(j)) - 2 !i.e le triangle numero Svoi4CN(i) est T1CN(a:a+2)
          if(T1CN(Sa)/=-1)then
              PRINT*,'tree sa',sa,T1CN(Sa),2*T1CN(Sa)-1
              SSCO3P(1) = T1CO(2*T1CN(Sa)-1);   SSCO3P(2) = T1CO(2*T1CN(Sa))
              SSCO3P(3) = T1CO(2*T1CN(Sa+1)-1); SSCO3P(4) = T1CO(2*T1CN(Sa+1))
              SSCO3P(5) = T1CO(2*T1CN(Sa+2)-1); SSCO3P(6) = T1CO(2*T1CN(Sa+2))
              SIN_CIRC = in_c_inscr(STXY,SSCO3P,SWIZ)
            SELECT CASE(SIN_CIRC)
            CASE(.TRUE.)
              ALLOCATE(tmp_SCAV(1:size(SCAV)))
              tmp_SCAV(:) = SCAV(:)
              DEALLOCATE(SCAV)
              ALLOCATE(SCAV(size(tmp_SCAV)+1))
              SCAV(1:size(tmp_SCAV))=tmp_SCAV(:)
              SCAV(size(tmp_SCAV)+1) = compa2(tmp_SCAV(:),T1CN(Sa:Sa+2),Sa)
              PRINT*,'tmp_SCAV',tmp_SCAV(:)
              PRINT*,'t1CN',T1CN(Sa:Sa+2)
              DEALLOCATE(tmp_SCAV)
              T1CN(Sa:Sa+2) = -1
              DO k=1,size(Svoi4CN)
                if(Svoi4CN(k)==voi(j))then
                  Svoi4CN(k)=-1
                end if
              END DO
              DO k=1,size(voi)
                if(voi(k)==voi(j))then
                  voi(k)=-1
                end if
              END DO
            CASE(.FALSE.)
              DO k=1,size(Svoi4CN)
                if(Svoi4CN(k)==voi(j))then
                  Svoi4CN(k)=-2
                end if
              END DO
              DO k=1,size(voi)
                if(voi(k)==voi(j))then
                  voi(k)=-2
                end if
              END DO
            END SELECT
        end if
        end if
      END DO
      i = i +1
    END DO
  end subroutine tree

  integer function compa2(tm,tp,aa)
    implicit none
    integer, dimension(:), intent(in) :: tm
    integer, dimension(3), intent(in) :: tp !tp eqv current t1cn, tm eqv tmp_scav
    integer, intent(in) :: aa
    integer :: flog, clog
    integer :: k, l
    l = 1; flog = 0; clog = 0
    compa2 = 0
    PRINT*,'c2,tp',tp(:)
    PRINT*,'c2,tm',tm(:)
    DO WHILE(l<=size(tp))
      k = 1; flog = 0; clog = 0
      print*,'l',l
      DO WHILE(k<=size(tm))
        print*,'k',k
        if(tp(l)==tm(k))then
          k = k +1; flog = flog + 1
        else
          k = k + 1; clog = l
        end if
      END DO
      l = l +1
      if(flog == 0.AND.clog /=0)then
        compa2 = T1CN(aa+ clog-1)
        PRINT*,'$$'
        l = size(tp)+1
      end if
  END DO
  ! if(tm(k)/=tp(1))then
  !   compa2 = T1CN(aa) ! soit aa, soit aa+1, soit aa+2
  ! else if(tm(k)/=tp(2))then
  !   compa2 = T1CN(aa+1)
  ! else if(tm(k)/=tp(3))then
  !   compa2 = T1CN(aa+2)
  ! end if
  end function compa2

  integer function compa(tm,tp,aa)
    implicit none
    integer, dimension(3), intent(in) :: tm, tp !tp eqv current t1cn, tm eqv tmp_scav
    integer, intent(in) :: aa
    integer :: l
    compa = 0
      l = 1
      do while(l<4)
        if(tm(1)==tp(l).OR.tm(2)==tp(l).OR.tm(3)==tp(l))then
          l = l+1
        else
          compa = T1CN(aa+ l-1) ! soit aa, soit aa+1, soit aa+2
          l = 5
        end if
      end do
  end function compa
end module mod_NOYDELAU
