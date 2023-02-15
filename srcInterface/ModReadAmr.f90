!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModReadAmr

  ! reconstruct AMR grid and read data on this grid
  ! interpolate data

  use FL_BATL_lib, ONLY: MaxDim
  use ModUtilities, ONLY: CON_stop

  implicit none
  save

  private !except

  ! Public methods 
  public:: readamr_init   ! read header information
  public:: readamr_read   ! read AMR data
  public:: readamr_get    ! get data at some point
  public:: readamr_clean  ! clean all variables

  ! Public data (optional)
  real, public, allocatable:: State_VGB(:,:,:,:,:) ! variables stored in blocks

  integer, public:: nVar=0                 ! number of variables in State_VGB

  character(len=20), public:: TypeGeometry = '???'
  real, public:: CoordMin_D(MaxDim) = -0.5 ! lower domain limits in gen coords
  real, public:: CoordMax_D(MaxDim) = +0.5 ! upper domain limits in gen coords

  real,    public:: TimeData = -1.0        ! simulation time of data
  integer, public:: nStepData = 0          ! time step of data
  integer, public:: nVarData = 0           ! number of variables in data file 
  integer, public:: nBlockData = 0         ! number of blocks in data file

  integer, public:: nParamData = 0         ! number of parameters in data file
  real, allocatable, public:: ParamData_I(:) ! paramters in data file

  character(len=500), public:: NameVarData  = '???' ! all variable names
  character(len=500), public:: NameUnitData = '???' ! unit names

  ! Local variables
  character(len=20):: TypeGeometryBatl = '???'
  integer:: nCellData = 0 ! number of cells in data file
  integer:: nProcData ! number of processors that wrote data file
  logical:: IsBinary  ! if the unprocessed IDL files are in binary format
  integer:: nByteReal ! number of bytes for reals in unprocessed IDL file

  character(len=5):: TypeDataFile ! type of processed file (real4/real8/ascii)

contains
  !============================================================================
  subroutine readamr_init(NameFile, IsVerbose)

    use ModIoUnit, ONLY: UnitTmp_
    use FL_BATL_lib,  ONLY: MaxDim, nDim, nIjk, nIjk_D, iProc, nProc, iComm, &
         init_batl
    use BATL_grid, ONLY: create_grid
    use BATL_tree, ONLY: read_tree_file, distribute_tree
    use ModReadParam, ONLY: read_init, read_echo_set, read_file, &
         read_line, read_command, read_var, lStringLine

    character(len=*), intent(in):: NameFile  ! base name
    logical,          intent(in):: IsVerbose ! provide verbose output

    integer:: i, iDim, iError, nDimSim

    character(len=500):: NameFileOrig, NameHeaderFile

    integer:: MaxBlock
    integer:: nRgen=0
    real, allocatable:: Rgen_I(:)

    real:: CellSizePlot_D(MaxDim), CellSizeMin_D(MaxDim)

    integer:: nIjkIn_D(MaxDim),  nRoot_D(nDim)
    logical:: IsPeriodic_D(MaxDim), IsExist

    logical, parameter:: DoDebug = .false.
    
    character (len=lStringLine) :: NameCommand, StringLine

    character(len=*), parameter:: NameSub = 'readamr_init'
    !-------------------------------------------------------------------------
    if(DoDebug) &
         write(*,*) NameSub,' starting on proc=', iProc, ' nProc=', nProc

    NameHeaderFile = trim(NameFile)//'.info'
    inquire(file=NameHeaderFile, exist=IsExist)

    if(DoDebug) write(*,*) NameSub,' IsExist=', IsExist
    if(.not. IsExist) then
       NameHeaderFile = trim(NameFile)//'.h'
       inquire(file=NameHeaderFile, exist=IsExist)
    end if
    if(DoDebug) write(*,*) NameSub,' IsExist=', IsExist
    if(.not. IsExist) call CON_stop(NameSub// &
         ' ERROR: could not open '//trim(NameFile)//'.h or .info')

    if(DoDebug) write(*,*) NameSub,' NameHeaderFile=', trim(NameHeaderFile)
    call read_file(NameHeaderFile, iComm, IsVerbose=IsVerbose)
    if(DoDebug) write(*,*) NameSub,' read file done'
    call read_echo_set(IsVerbose)
    if(DoDebug) write(*,*) NameSub,' read_echo_set done'
    call read_init()
    if(DoDebug) write(*,*) NameSub,' read_init done'

    READPARAM: do
       if(.not.read_line(StringLine) )then
          EXIT READPARAM
       end if

       if(.not.read_command(NameCommand)) CYCLE READPARAM

       select case(NameCommand)
       case('#HEADFILE')
          call read_var('HeadFileName', NameFileOrig)
          call read_var('nProc', nProcData)
          call read_var('IsBinary', IsBinary)
          if(IsBinary) call read_var('nByteReal', nByteReal)
       case('#NDIM')
          call read_var('nDimSim', nDimSim)

       case('#GRIDBLOCKSIZE')
          nIjkIn_D = 1
          do i = 1, nDimSim
             call read_var('BlockSize', nIjkIn_D(i))
          enddo

       case('#ROOTBLOCK')
          nRoot_D = 1
          do i = 1, nDimSim
             call read_var('nRootBlock', nRoot_D(i))
          enddo

       case('#NSTEP')
          call read_var('nStep', nStepData)

       case('#TIMESIMULATION')
          call read_var('TimeSimulation', TimeData)

       case('#PLOTRANGE')
          do i = 1, nDimSim
             call read_var('CoordMin', CoordMin_D(i))
             call read_var('CoordMax', CoordMax_D(i))
          enddo

       case('#PLOTRESOLUTION')
          CellSizePlot_D = 0;
          do i = 1, nDimSim
             call read_var('DxSavePlot', CellSizePlot_D(i))
          enddo
          if(CellSizePlot_D(1) >= 0.0) call CON_stop(NameSub// &
               ': the resolution should be set to -1 for file'//trim(NameFile))

       case('#CELLSIZE')
          CellSizeMin_D = 0;
          do i = 1, nDimSim
             call read_var('CellSizeMin', CellSizeMin_D(i))
          enddo

       case('#NCELL')
          call read_var('nCellPlot',nCellData)

       case('#PLOTVARIABLE')
          call read_var('nPlotVar', nVarData)
          call read_var('NameVar', NameVarData)
          call read_var('NameUnit', NameUnitData)

       case('#SCALARPARAM')
          call read_var('nParam', nParamData)
          allocate(ParamData_I(nParamData))
          do i = 1, nParamData
             call read_var('Param', ParamData_I(i))
          enddo

       case('#GRIDGEOMETRYLIMIT')
          call read_var('TypeGeometry', TypeGeometry)
          if(index(TypeGeometry,'genr') > 0)then
             read(*,*) nRgen
             allocate(Rgen_I(nRgen))
             do i = 1, nRgen
                read(*,*) Rgen_I(i)
             end do
          else
             allocate(Rgen_I(1))
          end if

       case('#PERIODIC')
          do i = 1, nDimSim
             call read_var('IsPeriodic', IsPeriodic_D(i))
          enddo

       case('#OUTPUTFORMAT')
          call read_var('TypeFormat', TypeDataFile)

       case default
          ! write(*,*) 'WARNING: unknow command ', NameCommand
       end select

    enddo READPARAM

    ! Total number of blocks in the data file
    nBlockData = nCellData / nIjk

    ! Number of blocks per processor !!! this is wrong for box cut from sphere
    MaxBlock  = (nBlockData + nProc - 1) / nProc

    if(iProc == 0 .and. any(nIjkIn_D /= nIjk_D))then
       write(*,*) 'ERROR in ',NameSub,' while reading ',trim(NameHeaderFile)
       write(*,*) 'Block size in header file        =', nIjkIn_D
       write(*,*) 'READAMR is configured to nI,nJ,nK=', nIjk_D
       call CON_stop('Read other file or reconfigure and recompile READAMR')
    end if

    ! Initialize BATL (using generalized coordinates and radians)
    if(TypeGeometry(1:9)=='spherical') then
       TypeGeometryBatl = 'rlonlat'//TypeGeometry(10:20)
    else
       TypeGeometryBatl = TypeGeometry
    end if

    call init_batl(CoordMin_D, CoordMax_D, MaxBlock, &
         TypeGeometryBatl, rGenIn_I=rGen_I, nRootIn_D=nRoot_D, &
         IsPeriodicIn_D=IsPeriodic_D, &
         UseRadiusIn=.false., UseDegreeIn=.false.)

    ! Read the full tree information and create grid
    call read_tree_file(trim(NameFile)//'.tree')
    call distribute_tree(.true.)
    call create_grid

    deallocate(Rgen_I)

  end subroutine readamr_init
  !============================================================================
  subroutine readamr_read(NameFile, iVarIn_I, IsNewGridIn, IsVerboseIn, &
       UseCoordTest)

    use ModKind,  ONLY: Real4_, Real8_
    use FL_BATL_lib, ONLY: nDim, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock, nG, iProc, Xyz_DGB, &
         find_grid_block, message_pass_cell, xyz_to_coord
    use ModPlotFile, ONLY: read_plot_file
    use ModIoUnit,   ONLY: UnitTmp_
    use ModConst,    ONLY: cPi

    character(len=*), intent(in):: NameFile     ! data file name
    integer, optional, intent(in):: iVarIn_I(:) ! index of variables to store
    logical, optional, intent(in):: IsNewGridIn ! new grid (read info/tree)
    logical, optional, intent(in):: IsVerboseIn ! provide verbose output
    logical, optional, intent(in):: UseCoordTest! store cos^2(coord) into State

    logical:: IsNewGrid
    logical:: IsVerbose

    ! Allocatable arrays for holding linear file data
    real, allocatable  :: State_V(:), State_VI(:,:), Xyz_DI(:,:)

    real:: Xyz_D(MaxDim) = 0.0, Coord_D(MaxDim)

    real(Real4_)              :: DxCell4, Xyz4_D(3)
    real(Real4_), allocatable :: State4_V(:)
    real(Real8_)              :: DxCell8, Xyz8_D(3)
    real(Real8_), allocatable :: State8_V(:)

    integer:: iCell, iCell_D(MaxDim), i, j, k, l, iBlock, iProcFound, iError

    integer:: nVarLast = -1, MaxBlockLast = -1

    character(len=1):: StringTmp

    character(len=*), parameter:: NameSub = 'readamr_read'
    !--------------------------------------------------------------------------
    IsNewGrid = .true.
    if(present(IsNewGridIn)) IsNewGrid = IsNewGridIn

    IsVerbose = .false.
    if(present(IsVerboseIn)) IsVerbose = IsVerboseIn

    if(IsVerbose)write(*,*) NameSub,&
         ' starting with IsNewGrid, UseCoordTest=', &
         IsNewGrid, present(UseCoordTest)

    ! Find extension in filename
    l = index(NameFile,".",BACK=.true.)

    ! Read grid info if necessary
    if(IsNewGrid) call readamr_init(NameFile(1:l-1), IsVerboseIn)

    if(NameFile(l:len_trim(NameFile)) == '.idl')then
       if(.not.IsBinary)then
          TypeDataFile = 'idl'
       elseif(nByteReal == 4)then
          TypeDataFile = 'idl4'
       else
          TypeDataFile = 'idl8'
       end if
    end if

    if(IsVerbose)write(*,*) NameSub,' TypeDataFile=', TypeDataFile

    nVar = nVarData
    if(present(iVarIn_I))then
       nVar = size(iVarIn_I)
       if(nVar>nVarData .or. any(iVarIn_I<1) .or. any(iVarIn_I>nVarData))then
          write(*,*)'ERROR: nVarData, iVarIn_I=', nVarData, iVarIn_I
          call CON_stop(NameSub//': invalid iVarIn_I array')
       end if
    end if
    if(IsVerbose)write(*,*) NameSub,' nVarData, nVar, present(iVarIn_I)=', &
         nVarData, nVar, present(iVarIn_I)

    if(.not. allocated(State_VGB) .or. &
         nVar /= nVarLast .or. MaxBlock /= MaxBlockLast)then
       if(allocated(State_VGB)) deallocate(State_VGB)
       allocate(&
            State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       !if(IsVerbose)write(*,*) NameSub,' allocated State_VGB(', &
       !     nVar,',',MinI,':',MaxI,',',MinJ,':',MaxJ,',',MinK,':',MaxK, &
       !     ',',MaxBlock,')'

       nVarLast = nVar; MaxBlockLast = MaxBlock
    end if
    State_VGB = 0.0

    ! The tests need at least nDim variables to be stored
    if(present(UseCoordTest)) nVar = max(nVar, MaxDim)

    select case(TypeDataFile)
    case('ascii')
       open(UnitTmp_, file=NameFile, status='old', iostat=iError)
       if(iError /= 0) call CON_stop(NameSub// &
            ' ERROR: could not open ascii file '//trim(NameFile))
       ! Read and discard 4 (nParamData==0) or 5 header lines 
       do i = 1, 4 + min(nParamData, 1)
          read(UnitTmp_,*) StringTmp
       end do
       if(IsVerbose)write(*,*) NameSub,' read header lines from ascii file'
    case('real4', 'real8')
       allocate(State_VI(nVarData,nCellData), Xyz_DI(nDim,nCellData))
       call read_plot_file(NameFile, TypeFileIn=TypeDataFile, &
            CoordOut_DI=Xyz_DI, VarOut_VI = State_VI)
    case('idl')
       open(UnitTmp_, file=NameFile, status='old', iostat=iError)
       if(iError /= 0) call CON_stop(NameSub// &
            ' ERROR: could not open ascii IDL file '//trim(NameFile))
    case('idl4', 'idl8')
       open(UnitTmp_, file=NameFile, form='unformatted', status='old', &
            iostat=iError)
       if(iError /= 0) call CON_stop(NameSub// &
            ' ERROR: could not open binary IDL file '//trim(NameFile))
    case default
       call CON_stop(NameSub//': unknown TypeDataFile='//trim(TypeDataFile))
    end select

    ! State variables for one cell
    allocate(State_V(nVarData), State4_V(nVarData), State8_V(nVarData))

    ! put each data point into the tree
    do iCell = 1, nCellData
       ! get cell data
       select case(TypeDataFile)
       case('ascii')
          read(UnitTmp_,*) Xyz_D(1:nDim), State_V
       case('real4', 'real8')
          Xyz_D(1:nDim) = Xyz_DI(:,iCell)
          State_V = State_VI(:,iCell)
       case('idl')
          read(UnitTmp_,*) DxCell8, Xyz_D, State_V
       case('idl4')
          read(UnitTmp_) DxCell4, Xyz4_D, State4_V
          Xyz_D =  Xyz4_D; State_V = State4_V
       case('idl8')
          read(UnitTmp_) DxCell8, Xyz8_D, State8_V
          Xyz_D =  Xyz8_D; State_V = State8_V
       case default
          call CON_stop(NameSub//': unknown TypeDataFile='//trim(TypeDataFile))
       end select

       ! Find cell
       call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D)

       if(iBlock < 0)then
          write(*,*)'ERROR for iCell, Xyz_D=', iCell, Xyz_D
          call CON_stop(NameSub//': could not find cell on the grid')
       end if

       !check if point belongs on this processor
       if (iProcFound /= iProc) CYCLE

       i = iCell_D(1); j = iCell_D(2); k = iCell_D(3)

       if(any(abs(Xyz_DGB(:,i,j,k,iBlock) - Xyz_D) > 1e-5))then
          write(*,*)NameSub,' ERROR at iCell,i,j,k,iBlock,iProc=', &
               iCell, i, j, k, iBlock, iProc
          write(*,*)NameSub,' Xyz_D  =', Xyz_D
          write(*,*)NameSub,' Xyz_DGB=', Xyz_DGB(:,i,j,k,iBlock)
          call CON_stop(NameSub//': incorrect coordinates')
       end if

       if(present(iVarIn_I))then
          State_VGB(:,i,j,k,iBlock) = State_V(iVarIn_I)
       else
          State_VGB(:,i,j,k,iBlock) = State_V
       end if

       ! For verification tests
       if(present(UseCoordTest))then
          ! Store cos^2 of generalized coordinates into first MaxDim elements
          call xyz_to_coord(Xyz_D, Coord_D)
          State_VGB(1:MaxDim,i,j,k,iBlock) = &
               cos(cPi*(Coord_D - CoordMin_D)/(CoordMax_D - CoordMin_D))**2
       end if

    enddo

    if(TypeDataFile=='ascii' .or. TypeDataFile(1:3)=='idl') close(UnitTmp_)

    if(IsVerbose)write(*,*)NameSub,' read data'

    ! deallocate to save memory
    deallocate(State_V, State4_V, State8_V)
    if(allocated(State_VI)) deallocate (State_VI, Xyz_DI) 

    ! Set ghost cells if any. Note that OUTER ghost cells are not set!
    if(nG > 0) call message_pass_cell(nVar, State_VGB)

    if(IsVerbose)then
       write(*,*)NameSub,' done'
       flush(6)
    endif

  end subroutine readamr_read

  !============================================================================
  subroutine readamr_get(Xyz_D, State_V, IsFound, CellSize_D)

    use FL_BATL_lib, ONLY: nDim, nG, nIJK_D, iProc, Xyz_DGB, CellSize_DB, &
         iTree_IA, Level_, MaxCoord_I, IsCartesianGrid, &
         interpolate_grid, interpolate_grid_amr, find_grid_block

    real,    intent(in)  :: Xyz_D(MaxDim)   ! location on grid
    real,    intent(out) :: State_V(0:nVar) ! weight and variables
    logical, intent(out) :: IsFound         ! true if found on grid

    real, optional, intent(out):: CellSize_D(3) ! cell size

    ! Block and processor index for the point
    integer:: iBlock, iProcOut, iNode, iLevel

    ! Variables for linear interpolation using ghost cells
    integer:: i1, j1=1, k1=1, i2, j2, k2
    real:: Dist_D(MaxDim), Dx1, Dx2, Dy1, Dy2, Dz1, Dz2

    ! Variables for AMR interpolation without ghost cells
    integer:: iCell, nCell, iCell_II(0:nDim,2**nDim), iCell_D(MaxDim), i, j, k
    real:: Weight_I(2**nDim)

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub='readamr_get'
    !-------------------------------------------------------------------------
    if(DoDebug)write(*,*)NameSub,' starting with Xyz_D=', Xyz_D

    State_V = 0.0

    call find_grid_block(Xyz_D, iProcOut, iBlock, iCell_D, Dist_D, iNode)
    if(DoDebug)write(*,*)NameSub,&
         ' found iProcOut, iBlock, iNode, iCell_D, Dist_D=', &
         iProcOut, iBlock, iNode, iCell_D, Dist_D

    IsFound = iBlock > 0

    if(present(CellSize_D))then
       if(IsFound)then
          iLevel = iTree_IA(Level_,iNode)
          CellSize_D = (CoordMax_D - CoordMin_D)/nIjk_D/MaxCoord_I(iLevel)
       else
          CellSize_D = 0.0
       end if
    end if

    if(.not.IsFound) RETURN

    ! Check if all surrounding cells are inside a single block
    if(all(iCell_D(1:nDim) > 0 .and. iCell_D(1:nDim) < nIJK_D(1:nDim)))then

       if(iProcOut /= iProc) RETURN

       ! Set weight to 1.0
       State_V(0) = 1.0

       ! Set indexes and distances for interpolation
       Dx1 = Dist_D(1); Dx2 = 1 - Dx1
       i1  = iCell_D(1); i2 = i1 + 1
       if(nDim > 1)then
          Dy1 = Dist_D(2); Dy2 = 1 - Dy1
          j1 = iCell_D(2); j2 = j1 + 1
       end if
       if(nDim > 2)then
          Dz1 = Dist_D(3); Dz2 = 1 - Dz1
          k1 = iCell_D(3); k2 = k1 + 1
       end if

       ! Interpolate
       if(nDim == 1)then
          State_V(1:nVar) = Dx2*State_VGB(:,i1,j1,k1,iBlock)  &
               +            Dx1*State_VGB(:,i2,j1,k1,iBlock)
       end if
       if(nDim == 2)then
          State_V(1:nVar) = Dy2*(Dx2*State_VGB(:,i1,j1,k1,iBlock)   &
               +                 Dx1*State_VGB(:,i2,j1,k1,iBlock))  &
               +            Dy1*(Dx2*State_VGB(:,i1,j2,k1,iBlock)   &
               +                 Dx1*State_VGB(:,i2,j2,k1,iBlock))
       end if
       if(nDim == 3)then
          State_V(1:nVar) = Dz2*(Dy2*(Dx2*State_VGB(:,i1,j1,k1,iBlock)   &
               +                      Dx1*State_VGB(:,i2,j1,k1,iBlock))  &
               +                 Dy1*(Dx2*State_VGB(:,i1,j2,k1,iBlock)   &
               +                      Dx1*State_VGB(:,i2,j2,k1,iBlock))) &
               +            Dz1*(Dy2*(Dx2*State_VGB(:,i1,j1,k2,iBlock)   &
               +                      Dx1*State_VGB(:,i2,j1,k2,iBlock))  &
               +                 Dy1*(Dx2*State_VGB(:,i1,j2,k2,iBlock)   &
               +                      Dx1*State_VGB(:,i2,j2,k2,iBlock)))
       end if

       if(DoDebug)then
          write(*,*)'!!! i1,j1,k1,i2,j2,k2=',i1,j1,k1,i2,j2,k2
          write(*,*)'!!! Dx1,Dx2,Dy1,Dy2,Dz1,Dz2=',Dx1,Dx2,Dy1,Dy2,Dz1,Dz2
          write(*,*)'!!! State_VGB(1,i1:i2,j1:j2,k1:k2,iBlock) = ',&
               State_VGB(1,i1:i2,j1:j2,k1:k2,iBlock)
       end if

    else
       ! Use interpolation algorithm that does not rely on ghost cells at all
       if(IsCartesianGrid)then
          ! Failed test on spherical grid
          call interpolate_grid_amr(Xyz_D, nCell, iCell_II, Weight_I)
       else
          ! Simple algorithm at resolution changes
          call interpolate_grid(Xyz_D, nCell, iCell_II, Weight_I)
       end if

       if(DoDebug)write(*,*)NameSub,': interpolate iProc, nCell=',iProc, nCell

       do iCell = 1, nCell
          iBlock  = iCell_II(0,iCell)
          iCell_D = 1
          iCell_D(1:nDim) = iCell_II(1:nDim,iCell)
          i      = iCell_D(1)
          j      = iCell_D(2)
          k      = iCell_D(3)
          if(DoDebug)write(*,*)NameSub,': iProc,iBlock,i,j,k,=',&
               iProc, iBlock, i, j, k

          State_V(0) = State_V(0)  + Weight_I(iCell)
          State_V(1:nVar) = State_V(1:nVar) &
               + Weight_I(iCell)*State_VGB(:,i,j,k,iBlock)

          if(DoDebug)write(*,*)NameSub, ': iProc,iBlock,i,j,k,Xyz,State=', &
               iProc, iBlock, i, j, k, &
               Xyz_DGB(:,i,j,k,iBlock), State_VGB(:,i,j,k,iBlock)
       end do
    end if

    if(DoDebug)write(*,*)NameSub,' finished with State_V=', State_V

  end subroutine readamr_get
  !============================================================================
  subroutine readamr_clean
    use FL_BATL_lib, ONLY: clean_batl

    call clean_batl
    if(allocated(State_VGB)) deallocate(State_VGB)
    if(allocated(ParamData_I)) deallocate(ParamData_I)
    nVar       = 0
    nVarData   = 0
    nBlockData = 0

  end subroutine readamr_clean

end module ModReadAmr

