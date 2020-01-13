module mod_VERIF
  USE mod_parameters
  Implicit None

  CONTAINS
    Integer function WR_verif_CNCO(FMESH_QK)
      Implicit None
      type(quality), intent(IN) :: FMESH_QK
      Integer :: i, sco, sco2
      sco = size(T1CO);sco2 = size(T1CO)/2
      WR_verif_CNCO = 0
      OPEN(20,file="VERIF_CNCO.txt")
      WRITE(20,*)"QUALITY: moyenne, min, max:"
      WRITE(20,*)FMESH_QK%mo, FMESH_QK%mi, FMESH_QK%ma
      WRITE(20,*)"QUALITY PERCENT: [1:2], ]2:3], ]3:5], ]5:10], ]10:50], ]50:100], ]100:infty[ "
      WRITE(20,*)FMESH_QK%Qpc(1:7)
      WRITE(20,*)"HEAD ET VOI:"
      DO i = 1, NV
         WRITE(20,*)head(i), voi(head(i-1)+1:head(i))
      END DO
      CLOSE(20)
      OPEN(30,file="CNCO.mesh")
      WRITE(30,*)'MeshVersionFormatted'
      WRITE(30,*)"Dimension"
      WRITE(30,*)3
      WRITE(30,*)"Vertices"
      WRITE(30,*)sco2
      do i = 1, (sco-2)+1, 2
        WRITE(30,*) T1CO(i), T1CO(i+1), 0
      end do
      WRITE(30,*)"Edges"
      WRITE(30,*)size(T1ED_L)/3
      do i = 1, (size(T1ED_L)-3)+1, 3
        WRITE(30,*)T1ED_L(i:i+2)
      end do
      WRITE(30,*)"Triangles"
      WRITE(30,*)size(T1CN)/3
      do i = 1, (size(T1CN)-3)+1, 3
        WRITE(30,*) T1CN(i:i+2)
      end do
      WRITE(30,*)"End"
      CLOSE(30)
      WR_verif_CNCO = 1
    End function WR_verif_CNCO



    Integer function WR_verif_T2DCO(FMESH_QK)
      Implicit None
      type(quality), intent(IN) :: FMESH_QK
      Integer :: i, sco, sco2
      sco = size(T1CO);sco2 = size(T1CO)/2
      WR_verif_T2DCO = 0
      OPEN(20,file="VERIF_T2DCO.txt")
      WRITE(20,*)"QUALITY: moyenne, min, max:"
      WRITE(20,*)FMESH_QK%mo, FMESH_QK%mi, FMESH_QK%ma
      WRITE(20,*)"QUALITY PERCENT: [1:2], ]2:3], ]3:5], ]5:10], ]10:50], ]50:100], ]100:infty[ "
      WRITE(20,*)FMESH_QK%Qpc(1:7)
      WRITE(20,*)"HEAD ET VOI:"
      DO i = 1, NV
         WRITE(20,*)head(i), voi(head(i-1)+1:head(i))
      END DO
      CLOSE(20)
      OPEN(30,file="T2DCO.mesh")
      WRITE(30,*)'MeshVersionFormatted'
      WRITE(30,*)"Dimension"
      WRITE(30,*)3
      WRITE(30,*)"Vertices"
      WRITE(30,*)sco2
      do i = 1, (sco-2)+1, 2
        WRITE(30,*) T1CO(i), T1CO(i+1), 0
      end do
      WRITE(30,*)"Edges"
      WRITE(30,*)size(T1ED_L)/3
      do i = 1, (size(T1ED_L)-3)+1, 3
        WRITE(30,*)T1ED_L(i:i+2)
      end do
      WRITE(30,*)"Triangles"
      WRITE(30,*)size(T2D)
      do i = 1, size(T2D)
        WRITE(30,*) T2D(i)%TCN(:)
      end do
      WRITE(30,*)"End"
      CLOSE(30)
      WR_verif_T2DCO = 1
    End function WR_verif_T2DCO

end module mod_VERIF
