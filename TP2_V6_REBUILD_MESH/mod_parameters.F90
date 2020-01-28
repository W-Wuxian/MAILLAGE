module mod_parameters
  Implicit None

  INTEGER, PARAMETER :: PR = 8
  INTEGER, PARAMETER :: STCO = 6, STCN = 3
  INTEGER, PARAMETER :: NBPELE = 3, NBPEDG = 2
  REAL(PR), PARAMETER :: SWIZ = 0.000000000001_PR
  INTEGER :: LB_head, UB_head, UB_voi !WILL BE SET WITH Subroutine BUILD_HEAD_VOI_CN
  INTEGER :: imp_ct

  type quality
    real(PR) :: mo, mi, ma
    real(PR), dimension(7) :: Qpc
  end type quality

  type TAR_GET
    real(PR), dimension(2) :: XY
  end type TAR_GET

type ISBASE
  integer :: NUPB
  integer, dimension(3)  :: PBCN
  real(PR), dimension(6) :: PBCO
end type ISBASE

  INTEGER :: DIM ,NV, NE, NT
  REAL(PR),     DIMENSION(:),   ALLOCATABLE :: T1CO    !TABLE DE COORDONNEES
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: T1CN    !TABLE DE CONNECTIVITE
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: T1ED    !TABLE EDGES
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: T1ED_L  !TABLE EDGES AVEC LES LABELS

  INTEGER :: FCRTL !CONTROLE VERIF.txt
  CHARACTER(LEN=50), PARAMETER :: MESH_FILE =  "/mnt/d/PROJET_MESH/meshes/OWN/"//"CARRE2.mesh"

  type(quality) :: MESH_QK_CNCO
  type(TAR_GET), dimension(:), allocatable :: TARG
  INTEGER, DIMENSION(:), ALLOCATABLE :: head, voi

  INTEGER, DIMENSION(:), ALLOCATABLE :: CAV

  !MESH_FILE = "/mnt/d/PROJET_MESH/meshes/OWN/"//"CARRE2.mesh"

end module mod_parameters
