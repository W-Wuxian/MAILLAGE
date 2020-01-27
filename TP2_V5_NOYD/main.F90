program main
  USE mod_parameters
  USE mod_VERIF
  USE mod_VOI
  USE mod_CONNEC_COO
  USE mod_SEARCH_IN
  USE mod_NOYDELAU
  Implicit None
  logical  :: test_in_c_inscr, test_in_tri
  Real(PR),dimension(2) :: TXY
  Real(PR) :: WIZ
  Integer  :: i, j, UB
  INTEGER :: CRTL_NOY
  Real(PR), Dimension(:), Allocatable :: CO3P
  Integer, Dimension(:), Allocatable :: head4CN, voi4CN

  !----------------------------------------------------------------------
  PRINT*,'DEFINE:', STCO,STCN,NBPELE,NBPEDG
  CALL RD_MF_CONNEC_COO(MESH_FILE,DIM,NV,NE,NT,T1CO,T1CN,T1ED,T1ED_L)
  CALL QK_INFO_CNCO(T1CN, T1CO, MESH_QK_CNCO)
  CALL BUILD_HEAD_VOI_CN(T1CN,head,voi)
  FCRTL = WR_verif_CNCO(MESH_QK_CNCO)
  CALL BUILD_CN_VOI_CN(T1CN, head, voi, head4CN, voi4CN)
  OPEN(30,file='verftri.txt')
  DO i = 1, NT
     WRITE(30,*)head4CN(i), voi4CN(head4CN(i-1)+1:head4CN(i))
  END DO
  CLOSE(30)
  DEALLOCATE(head, voi)

  TXY(1) = 0.3_PR; TXY(2) = 0.22_PR
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
  !TEST IF TARGET IS IN CERCLE CIRCONSCRIT AND IS IN TRIANGLE
  PRINT*,'TARGET TXY',TXY(1),TXY(2),':'
  DO i = 1, UB, NBPELE
     j = j + 1
     CO3P(1) = T1CO(2*T1CN(i)-1);   CO3P(2) = T1CO(2*T1CN(i))
     CO3P(3) = T1CO(2*T1CN(i+1)-1); CO3P(4) = T1CO(2*T1CN(i+1))
     CO3P(5) = T1CO(2*T1CN(i+2)-1); CO3P(6) = T1CO(2*T1CN(i+2))
     test_in_c_inscr = in_c_inscr(TXY,CO3P,WIZ)
     test_in_tri     = IsInsideTriangle(CO3P,TXY)!in_triangle(TXY,CO3P)
     SELECT CASE(test_in_c_inscr)
     CASE(.TRUE.)
        PRINT*,'EST DANS LE CERCLE CIRCONSCRIT AU TRIANGLE NUM=',j,'test',test_in_c_inscr
        OPEN(60,file='CERCLE_INSC_TRUE.txt',POSITION='APPEND')
        WRITE(60,*)CO3P(1:6)
        CLOSE(60)
     END SELECT
     SELECT CASE(test_in_tri)
     CASE(.TRUE.)
        PRINT*,'EST DANS LE TRIANGLE NUM=',j,'test',test_in_tri
     END SELECT
  END DO
  OPEN(10,file='in_c_inscr.txt',position='append')
  WRITE(10,*)' '
  WRITE(10,*)' '
  WRITE(10,*)TXY(1),TXY(2)
  CLOSE(10)
  !ECRITURE DU SCRIPT .plt POUR AFFICHER LE MAILLAGE LU ET LE POINT TARGET ET CERCLE CIRCONSCRIT...
  OPEN(30,file='geo.plt',position='append')
  WRITE(30,*)'plot "in_c_inscr.txt" index 0 u 1:2 ps 3 pt 2 title "pt tri" , "in_c_inscr.txt" index 1 u 1:2 ps 4 pt 3 title "target"  '
  WRITE(30,*)'replot "edg.txt" u 1:2:($3-$1):($4-$2) with vectors arrowstyle 1 title "edge" '
  WRITE(30,*)'replot "cercle.txt" with circles'
  WRITE(30,*)'replot "cercle.txt" u 1:2 ps 2 pt 4 title "center" '
  CLOSE(30)

  DEALLOCATE(CO3P, head4CN, voi4CN)
  ALLOCATE(TARG(1))
  TARG(1)%XY(1)=TXY(1);TARG(1)%XY(2)=TXY(2)
  CRTL_NOY = NOYDELAU(TARG)
  DO i=1,size(CAV)
    PRINT*,'i, CAV',i,CAV(i)
  END DO
  DEALLOCATE(CAV,TARG)
  DEALLOCATE(T1CO, T1CN, T1ED, T1ED_L)


end program main
