module mod_parameters
  Implicit None

  INTEGER, PARAMETER :: PR = 8
  INTEGER, PARAMETER :: STCO = 6, STCN = 3
  INTEGER, PARAMETER :: NBPELE = 3, NBPEDG = 2

  type tab_2D
     !real(PR), dimension(6) :: TCO !X1,Y1,X2,Y2,X3,Y3
     integer,  dimension(3) :: TCN !P1,P2,P3
     !integer,                  :: TTAG
  end type tab_2D

  type quality
    real(PR) :: mo, mi, ma
    real(PR), dimension(7) :: Qpc
  end type quality

  INTEGER :: DIM ,NV, NE, NT
  REAL(PR),     DIMENSION(:),   ALLOCATABLE :: T1CO    !TABLE DE COORDONNEES
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: T1CN    !TABLE DE CONNECTIVITE
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: T1ED    !TABLE EDGES
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: T1ED_L  !TABLE EDGES AVEC LES LABELS


  INTEGER :: FCRTL !CONTROLE VERIF.txt
  CHARACTER(LEN=50), PARAMETER :: MESH_FILE =  "/mnt/d/PROJET_MESH/meshes/OWN/"//"CARRE3.mesh"

  type(tab_2D), dimension(:), allocatable :: T2D
  type(quality) :: MESH_QK_CNCO,MESH_QK_T2DCO

  INTEGER, DIMENSION(:), ALLOCATABLE :: head, voi

  !MESH_FILE = "/mnt/d/PROJET_MESH/meshes/OWN/"//"CARRE2.mesh"

end module mod_parameters
