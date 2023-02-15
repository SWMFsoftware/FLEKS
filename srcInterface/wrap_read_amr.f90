!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! Extrernal subroutines to call from other languages
!==============================================================================
subroutine wrapamr_init_mpi() bind(C)
  use FL_BATL_lib,  ONLY: init_mpi

  ! Set the MPI parameters inside BATL
  implicit none

  call init_mpi

end subroutine wrapamr_init_mpi
!==============================================================================
subroutine wrapamr_clean() bind(C)

  ! Deallocate all memory used by READAMR

  use ModReadAmr, ONLY: readamr_clean
  
  call readamr_clean
  
end subroutine wrapamr_clean
!=============================================================================
subroutine wrapamr_read_header(NameFile_I, l, lVerbose) bind(C)

  ! Read header information from .h or .info
  ! This sets number and names of variables, domain size, etc.
  ! l is the maximum length of the filename
  ! If lVerbose is 1, write out verbose information.

  use ModUtilities, ONLY: char_array_to_string
  use ModReadAmr, ONLY:readamr_init

  use iso_c_binding, ONLY: c_char, c_int

  implicit none

  character(c_char), intent(in):: NameFile_I(*)
  integer(c_int), value, intent(in):: l
  integer(c_int), value, intent(in):: lVerbose

  character(len=l):: NameFile
  integer:: i
  !----------------------------------------------------------------------------
  call char_array_to_string(NameFile_I, NameFile)

  ! Cut off extension from the file name (if any)
  i = index(NameFile, ".", BACK=.true.)

  call readamr_init(NameFile(1:i-1), IsVerbose=lVerbose==1_c_int)

end subroutine wrapamr_read_header
  
!==============================================================================
subroutine wrapamr_read_file(NameFile_I, l, lNewGrid, lVerbose) bind(C)

  ! Read data from the file given by the C character array NameFile_I.
  ! The maximum length of the file name is given by.
  ! If lNewGrid == 1 then read all information, otherwise data only
  ! If lVerbose == 1 then write out verbose information

  use ModUtilities, ONLY: char_array_to_string
  use ModReadAmr, ONLY:readamr_read

  use iso_c_binding, ONLY: c_char, c_int

  implicit none

  character(c_char), intent(in):: NameFile_I(*)
  integer(c_int), value, intent(in):: l
  integer(c_int), value, intent(in):: lNewGrid
  integer(c_int), value, intent(in):: lVerbose

  character(len=l):: NameFile
  !----------------------------------------------------------------------------
  call char_array_to_string(NameFile_I, NameFile)

  call readamr_read(NameFile, &
       IsNewGridIn = lNewGrid==1_c_int, IsVerboseIn=lVerbose==1_c_int)

end subroutine wrapamr_read_file
!==============================================================================
subroutine wrapamr_get_ndim(nDimOut) bind(C)

  ! Get number of dimensions of the grid

  use FL_BATL_lib, ONLY: nDim
  use iso_c_binding, ONLY: c_int
  implicit none
  integer(c_int), intent(out):: nDimOut
  
  nDimOut = nDim
end subroutine wrapamr_get_ndim    
!==============================================================================
subroutine wrapamr_get_block_size(nIOut, nJOut, nKOut, nGOut) bind(C)

  ! Get number of cells in the block: there are nI*nJ*nK physical cells 
  ! surrounded by nG ghost cells in the first nDim dimensions.

  use FL_BATL_lib, ONLY: nI, nJ, nK, nG
  use iso_c_binding, ONLY: c_int
  implicit none
  integer(c_int), intent(out):: nIOut, nJOut, nKOut, nGOut
  
  nIOut = nI; nJOut = nJ; nKOut = nK; nGOut = nG
end subroutine wrapamr_get_block_size
!==============================================================================
subroutine wrapamr_get_nvar(nVarOut) bind(C)

  ! Get number of variables per grid cell in the data file

  use ModReadAmr, ONLY: nVarData
  use iso_c_binding, ONLY: c_int
  implicit none
  integer(c_int), intent(out):: nVarOut
  
  nVarOut = nVarData
end subroutine wrapamr_get_nvar    
!==============================================================================
subroutine wrapamr_get_namevar(NameVarOut_I, l) bind(C)

  ! Names of variables that are returned by the interpolation routine 
  ! l is the length of the string (without trailing spaces)

  use ModUtilities, ONLY: string_to_char_array
  use ModReadAmr, ONLY: NameVarData

  use iso_c_binding, ONLY: c_char, c_int
  
  implicit none
  character(c_char), intent(out):: NameVarOut_I(*)
  integer(c_int),    intent(out):: l
  !--------------------------------------------------------------------------
  call string_to_char_array(NameVarData, NameVarOut_I)
  l = len_trim(NameVarData)

end subroutine wrapamr_get_namevar
!==============================================================================
subroutine wrapamr_get_nameunit(NameUnitOut_I, l) bind(C)

  ! units of the coordinates followed by the units of the interpolated
  ! variables
  ! l is the length of the string (without trailing spaces)

  use ModUtilities, ONLY: string_to_char_array
  use ModReadAmr, ONLY: NameUnitData

  use iso_c_binding, ONLY: c_char, c_int
  
  implicit none
  character(c_char), intent(out):: NameUnitOut_I(*)
  integer(c_int),    intent(out):: l

  !--------------------------------------------------------------------------
  call string_to_char_array(NameUnitData, NameUnitOut_I)
  l = len_trim(NameUnitData)
 
end subroutine wrapamr_get_nameunit
!==============================================================================
subroutine wrapamr_get_domain(CoordMinOut_D, CoordMaxOut_D) bind(C)

  ! Boundary of the computational domain in the 
  ! normalized/generalized coordinates

  use FL_BATL_lib, ONLY: nDim
  use ModReadAmr, ONLY: CoordMin_D, CoordMax_D
  use iso_c_binding, ONLY: c_double
  
  implicit none
  real(c_double), intent(out):: CoordMinOut_D(nDim), CoordMaxOut_D(nDim)

  CoordMinOut_D = CoordMin_D(1:nDim)
  CoordMaxOut_D = CoordMax_D(1:nDim)

end subroutine wrapamr_get_domain

!=============================================================================
subroutine wrapamr_get_block(x_D, iProcFound, iBlock, iLevel, &
     State_VG, Xyz_DG, CoordMinBlock_D, CoordMaxBlock_D) bind(C)

  ! Get data for a full grid block. iProcFound is the index of the
  ! processor that contains the point Xyz_D. It is -1 if the point 
  ! is not in the domain. All other variables are only set On the 
  ! processor iProcFound. iBlock is the block index on the processor.
  ! State_VG and Xyz_DG contain the data and the Cartesian coordinates 
  ! for all grid cell centers (including ghost cells) of the block.
  ! CoordMinBlock_D and CoordMaxBlock_D are the (generalized) coordinates
  ! of the block edges. iLevel is the grid level (0 is the root level).

  use FL_BATL_lib, ONLY: nDim, MaxDim, Xyz_DGB, CoordMin_DB, CoordMax_DB, iProc, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, iTree_IA, Level_, find_grid_block
  use ModReadAmr, ONLY: State_VGB, nVar
  use iso_c_binding, ONLY: c_double, c_int

  implicit none

  real(c_double), intent(in) :: x_D(nDim)
  integer(c_int), intent(out):: iProcFound, iBlock, iLevel
  real(c_double), intent(out):: State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
  real(c_double), intent(out):: Xyz_DG(nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
  real(c_double), intent(out):: CoordMinBlock_D(nDim), CoordMaxBlock_D(nDim)

  integer:: iNode
  
  real:: Xyz_D(MaxDim)
  !----------------------------------------------------------------------------
  ! This copy converts real precision and pads extra dimensions with 0
  Xyz_D = 0.0
  Xyz_D(1:nDim) = x_D

  call find_grid_block(Xyz_D, iProcFound, iBlock, iNodeOut=iNode)

  if(iProc == iProcFound)then
     State_VG        = State_VGB(:,:,:,:,iBlock)
     Xyz_DG          = Xyz_DGB(1:nDim,:,:,:,iBlock)
     CoordMinBlock_D = CoordMin_DB(1:nDim,iBlock)
     CoordMaxBlock_D = CoordMax_DB(1:nDim,iBlock)
     iLevel          = iTree_IA(Level_,iNode)
  else
     State_VG        = 0.0
     Xyz_DG          = 0.0
     CoordMinBlock_D = 0.0
     CoordMaxBlock_D = 0.0
     iLevel          = -1
  end if

end subroutine wrapamr_get_block
!=============================================================================
subroutine wrapamr_get_data(x_D, StateOut_V, iFound) bind(C)

  ! Get the interpolated values StateOut_V at the point XyzOut_V
  ! The first index of State_V is the interpolation weight, 
  ! so StateOut_V has nVar+1 elements.
  ! For parallel execution, MPI_SUM is needed and a division by total weight.
  ! iFound is set to 0 if point is not found (outside domain), 1 otherwise.

  use ModReadAmr, ONLY: nVar, readamr_get
  use FL_BATL_lib,   ONLY: nDim, MaxDim
  use iso_c_binding, ONLY: c_double, c_int
  
  implicit none
  real(c_double), intent(in) :: x_D(nDim)
  real(c_double), intent(out):: StateOut_V(0:nVar)
  integer(c_int), intent(out):: iFound

  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  logical:: IsFound
  !----------------------------------------------------------------------------

  ! This copy converts real precision if needed and pads extra dimensions
  Xyz_D = 0.0
  Xyz_D(1:nDim) = x_D
  call readamr_get(Xyz_D, State_V, IsFound)
  
  ! This copy converts real precision if needed
  StateOut_V = State_V

  ! Set integer found flag
  iFound = 0
  if(IsFound) iFound = 1

end subroutine wrapamr_get_data
!=============================================================================
subroutine wrapamr_get_data_cell(x_D, StateOut_V, &
     CellSizeOut_D, iFound) bind(C)

  ! Get the interpolated values StateOut_V at the point XyzOut_V
  ! iFound is set to 0 if point is not found (outside domain), 1 otherwise.
  ! Set CellSizeOut_D to size of the cell containing the point.

  use ModReadAmr, ONLY: nVar, readamr_get
  use FL_BATL_lib, ONLY: nDim, MaxDim
  use iso_c_binding, ONLY: c_double, c_int

  implicit none
  real(c_double), intent(in) :: x_D(nDim)
  real(c_double), intent(out):: StateOut_V(0:nVar)
  real(c_double), intent(out):: CellSizeOut_D(nDim)
  integer(c_int), intent(out):: iFound

  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  real   :: CellSize_D(3)
  logical:: IsFound
  !----------------------------------------------------------------------------

  ! This copy converts real precision if needed
  Xyz_D = 0.0
  Xyz_D(1:nDim) = x_D
  call readamr_get(Xyz_D, State_V, IsFound, CellSize_D)
 
  ! These copies convert real precision if needed
  StateOut_V = State_V
  CellSizeOut_D = CellSize_D(1:nDim)

  ! Set integer found flag
  iFound = 0
  if(IsFound) iFound = 1

end subroutine wrapamr_get_data_cell
!=============================================================================
subroutine wrapamr_get_data_serial(x_D, StateOut_V, iFound) bind(C)

  ! Get the interpolated values StateOut_V at the point XyzOut_V.
  ! Division by sum of weights is done internally, 
  ! so StateOut_V has nVar elements.
  ! iFound is set to 0 if point is not found (outside domain), 1 otherwise.

  use ModReadAmr, ONLY: nVar, readamr_get
  use FL_BATL_lib,   ONLY: nDim, MaxDim
  use iso_c_binding, ONLY: c_double, c_int
  
  implicit none
  real(c_double), intent(in) :: x_D(nDim)         ! point position
  real(c_double), intent(out):: StateOut_V(nVar)  ! interpolated values
  integer(c_int), intent(out):: iFound            ! 1 if point is found, 0 if not

  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  logical:: IsFound
  !----------------------------------------------------------------------------

  ! This copy converts real precision if needed
  Xyz_D = 0.0
  Xyz_D(1:nDim) = x_D
  call readamr_get(Xyz_D, State_V, IsFound)

  if(IsFound)then
     ! Divide by weight.
     ! Also this converts real precision if needed.
     StateOut_V = State_V(1:nVar)/State_V(0)
     iFound = 1
  else
     ! Set integer found flag
     StateOut_V = 0.0
     iFound = 0
  end if

end subroutine wrapamr_get_data_serial

!=============================================================================
subroutine wrapamr_get_array_serial(nPoint, x_DI, State_VI) bind(C)

  ! Get the interpolated values StateOut_V at the point XyzOut_V.
  ! Division by sum of weights is done internally, 
  ! so StateOut_V has nVar elements. 
  ! Points that are not found in the domain are filled with all zero values.

  use ModReadAmr, ONLY: nVar, readamr_get
  use FL_BATL_lib,   ONLY: nDim, MaxDim
  use iso_c_binding, ONLY: c_double, c_int
  
  implicit none

  integer(c_int), intent(in), value :: nPoint         ! number of points
  real(c_double), intent(in) :: x_DI(nDim,nPoint)     ! locations of points
  real(c_double), intent(out):: State_VI(nVar,nPoint) ! interpolated values

  integer:: i 
  real   :: State_V(0:nVar)
  real   :: Xyz_D(MaxDim)
  logical:: IsFound
  !-------------------------------------------------------------------------

  ! Initialize all coordinates
  Xyz_D = 0.0
  do i = 1, nPoint
      Xyz_D(1:nDim) = x_DI(:,i)   
      call readamr_get(Xyz_D, State_V, IsFound)
      if(IsFound)then
         State_VI(:,i) = State_V(1:nVar) / State_V(0)
      else
         State_VI(:,i) = 0.0
      end if
  end do

end subroutine wrapamr_get_array_serial
!=============================================================================

