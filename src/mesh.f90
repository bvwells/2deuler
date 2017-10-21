module MeshClass
  implicit none

  integer :: NoOfDimensions=2
  integer :: NoOfNodes
  integer :: NoOfCells
  integer :: MaxNoOfSurroundingCells
  integer :: NoOfNodesInCell
  integer :: NoOfEdges
  integer :: NoOfBoundaryCells
  integer :: MovingMeshChoice

  ! Names of mesh files
  character(len=40) :: MeshFilename
  character(len=40) :: CellFilename

  contains
    subroutine ReadMesh(x,xdot,Cells)
   !----------------------------------------------------------------------------
   !  Routine: ReadMesh
   !
   !  Description: Read in the mesh 
   !
   !----------------------------------------------------------------------------
      implicit none

      integer,dimension(:,:),pointer :: Cells
      double precision,dimension(:,:),pointer :: x,xdot
      integer :: i,j,ios

      ! Read in mesh node co-ordinates
      write(6,*) 'Reading file '//MeshFilename
      open(unit=10,file=MeshFilename,status='old',iostat=ios)
      if( ios.ne.0 ) stop 'Can not open file '//MeshFilename
      read(10,*) NoOfNodes
      allocate(x(NoOfNodes,NoOfDimensions))
      allocate(xdot(NoOfNodes,NoOfDimensions))
      xdot=0.0d0
      do i=1,NoOfNodes
         read(10,*) (x(i,j),j=1,NoOfDimensions)
      enddo
      close(10)

      ! Read in mesh cell data
      write(6,*) 'Reading file '//CellFilename
      open(unit=10,file=CellFilename,status='old',iostat=ios)
      if( ios.ne.0 ) stop 'Can not open file '//CellFilename
      read(10,*) NoOfCells,NoOfNodesInCell
      allocate(cells(NoOfCells,NoOfNodesInCell))
      do i=1,NoOfCells
         read(10,*) (Cells(i,j),j=1,NoOfNodesInCell)
      enddo
      close(10)

      return
    end subroutine ReadMesh

    subroutine WriteMesh(x,Cells)
   !----------------------------------------------------------------------------
   !  Routine: WriteMesh
   !
   !  Description: Write the mesh 
   !
   !----------------------------------------------------------------------------
      implicit none
      double precision,intent(in) :: x(NoOfNodes,NoOfDimensions)
      integer,intent(in) :: Cells(NoOfCells,NoOfNodesInCell)
      integer :: i,j,ios

      write(6,*) 'Writing mesh files...'
      open(unit=10,file='nodes.m',iostat=ios)
      if( ios.ne.0 ) stop 'Can not open file nodes.m'
      do i=1,NoOfNodes
         write(10,'(6f16.6)') (x(i,j),j=1,NoOfDimensions)
      enddo
      close(10)

      open(unit=10,file='cells.m',iostat=ios)
      if( ios.ne.0 ) stop 'Can not open file cells.m'
      do i=1,NoOfCells
         write(10,*) (Cells(i,j),j=1,NoOfNodesInCell)
      enddo
      close(10)

      return
    end subroutine WriteMesh

    subroutine ProcessMesh(Cells,Con,Surnodes,Edges)
      implicit none

      integer,dimension(NoOfCells,NoOfNodesInCell) :: Cells
      integer,dimension(:,:),pointer :: Con,Surnodes,Edges

      call max_no_of_cells(Cells)
      if( MovingMeshChoice.gt.1 ) then
         allocate(Con(NoOfNodes,MaxNoOfSurroundingCells+1))
         call connectivity(Con,Cells)
      endif
      allocate(Surnodes(NoOfNodes,MaxNoOfSurroundingCells))
      call connect(Surnodes,Cells)

      call CalculateNumberOfEdges(Cells,Edges)

      return
    end subroutine ProcessMesh

    function CalculateCellCentre(Cells,x,CellNo) result(CellCentre)
   !----------------------------------------------------------------------------
   !  Routine: CalculateCellCentre
   !
   !  Description: Calculate the centroid of a cell 
   !
   !----------------------------------------------------------------------------
      implicit none
    !---------------------------------------------------------------------------
      integer,intent(in) :: Cells(NoOfCells,NoOfNodesInCell)
      double precision,intent(in) :: x(NoOfNodes,NoOfDimensions)
      integer,intent(in) :: CellNo
    !---------------------------------------------------------------------------
      double precision :: CellCentre(NoOfDimensions)
      integer :: i,j,Node
    !---------------------------------------------------------------------------

      CellCentre = 0.0d0
      do i=1,NoOfNodesInCell
         Node = Cells(CellNo,i)
         do j=1,NoOfDimensions
            CellCentre(j) = CellCentre(j) + x(Node,j)
         enddo
      enddo
      CellCentre = CellCentre / dble(NoOfNodesInCell)

      return
    end function CalculateCellCentre

    function CalculateCellArea(Cells,x,CellNo) result(CellArea)
   !----------------------------------------------------------------------------
   !  Routine: CalculateCellAreas
   !
   !  Description: Calculate the cell area using the Surveyor's Formula
   !
   !----------------------------------------------------------------------------
      implicit none
    !---------------------------------------------------------------------------
      integer,intent(in) :: Cells(NoOfCells,NoOfNodesInCell)
      double precision,intent(in) :: x(NoOfNodes,NoOfDimensions)
      integer,intent(in) :: CellNo
    !---------------------------------------------------------------------------
      double precision :: CellArea
      integer :: i,Node1,Node2
    !---------------------------------------------------------------------------
      ! Surveyor's Formula
      ! Area = 0.5 * sum_{i=1}^{N} (x(i)*y(i+1) - x(i+1)*y(i) )
      CellArea = 0.0d0
      do i=1,NoOfNodesInCell
         ! x(i) * y(i+1) - x(i+1) * y(i)
         Node1 = Cells(CellNo,mod(i-1,NoOfNodesInCell)+1)
         Node2 = Cells(CellNo,mod(i  ,NoOfNodesInCell)+1)
         CellArea = CellArea + 0.5d0*(x(Node1,1)*x(Node2,2)-x(Node2,1)*x(Node1,2))
      enddo

      return
    end function CalculateCellArea

    function CalculateLength(x,y,n) result(Length)
   !----------------------------------------------------------------------------
   !  Routine: CalculateLength
   !
   !  Description: Calculate the length between to vectors
   !
   !----------------------------------------------------------------------------
      implicit none
      integer,intent(in) :: n
      double precision,intent(in),dimension(1:n) :: x,y
      double precision :: Length
      integer :: i
      Length = 0.0d0
      do i=1,n
         Length = Length + (x(i)-y(i))*(x(i)-y(i))
      enddo
      Length = sqrt(Length)
      return
    end function CalculateLength

    subroutine CalculateNumberOfEdges(Cells,Edges)
      implicit none
      integer,intent(in) :: Cells(NoOfCells,NoOfNodesInCell)
      integer,pointer :: Edges(:,:)
      integer :: i,j,k,iEdge,nEdge,Swap

      write(6,*) 'Generating edges in mesh...'

      ! First guess at number of edges in mesh
      NoOfEdges = NoOfNodesInCell*NoOfCells
      allocate(Edges(NoOfEdges,4))
      Edges=0
      iEdge=0
      do i=1,NoOfCells
         do j=1,NoOfNodesInCell
            iEdge=iEdge+1
            Edges(iEdge,1) = Cells(i,mod(j-1,NoOfNodesInCell)+1)
            Edges(iEdge,2) = Cells(i,mod(j  ,NoOfNodesInCell)+1)
            Edges(iEdge,3) = i
         enddo
      enddo

      ! Remove duplicate edges
      nEdge = 0
      do i=1,NoOfEdges
         if( Edges(i,1).eq.0 .and. Edges(i,2).eq.0 ) cycle
         do j=i+1,NoOfEdges
            if( (Edges(i,1).eq.Edges(j,1) .and. Edges(i,2).eq.Edges(j,2)) .or.&
               &(Edges(i,1).eq.Edges(j,2) .and. Edges(i,2).eq.Edges(j,1)) ) then
               Edges(i,4) = Edges(j,3)
               do k=1,4
                  Edges(j,k) = 0
               enddo
               exit
           endif
         enddo
      enddo

      ! Compress the edges
      nEdge = 1
      do i=1,NoOfEdges
         if( Edges(i,1).ne.0 .and. Edges(i,2).ne.0 ) then
            do j=1,4
               Swap = Edges(i,j)
               Edges(i,j)=0
               Edges(nEdge,j)=Swap
            enddo
            nEdge = nEdge + 1
         endif
      enddo

      return
    end subroutine CalculateNumberOfEdges

    
    subroutine max_no_of_cells(cells)
      implicit none
      integer,intent(in),dimension(NoOfCells,NoOfNodesInCell) :: cells
      integer,allocatable :: array(:)
      integer :: i, j

      allocate(array(NoOfNodes))
      array=0
      do i=1,NoOfCells
         do j=1,NoOfNodesInCell
            array(cells(i,j))=array(cells(i,j)) + 1
         enddo
      enddo
      MaxNoOfSurroundingCells=maxval(array)
      deallocate(array)
      return
    end subroutine max_no_of_cells



    subroutine connectivity(con,cells)
      implicit none

      integer,dimension(NoOfCells,NoOfNodesInCell),intent(in) :: cells
      integer,dimension(NoOfNodes,MaxNoOfSurroundingCells+1),intent(inout) :: con
      integer,allocatable :: count(:)
      logical :: location
      integer :: i, j, k, l

      allocate(count(NoOfNodes))

      count = 1
      con=0

      do i=1,NoOfCells
         do L=1,NoOfNodesInCell
            do j=1,NoOfNodesInCell
               location=.false.
               do k=1,MaxNoOfSurroundingCells+1
                  if( cells(i,j)==con(cells(i,L),k) ) then
                     location=.true.
                  endif
               enddo

               if( .not.location ) then
                  con(cells(i,L),count(cells(i,L))) = cells(i,j)
                  count(cells(i,L)) = count(cells(i,L)) +1
               endif
            enddo
         enddo
      enddo

      deallocate(count)

      return
    end subroutine connectivity

    subroutine find_cells(found,cells,i)

      implicit none
      integer,intent(in) :: i
      integer,intent(in),dimension(NoOfCells,NoOfNodesInCell) :: cells
      integer,intent(inout),dimension(NoOfNodesInCell,NoOfNodesInCell+1) :: found
      integer,dimension(NoOfNodesInCell+1) :: swap
      integer :: j, k, L, count, foundcount

      foundcount=1; count=1; found=0
      
      do j=1,NoOfCells
         if( j/=i ) then
            do k=1,NoOfNodesInCell
               count=1
               do L=1,NoOfNodesInCell
                  if( (cells(j,L)==cells(i,mod(k-1,NoOfNodesInCell)+1)).or.&
                     &(cells(j,L)==cells(i,mod(k,NoOfNodesInCell)+1)) ) then
                     count=count+1
                  endif
               enddo
               if( count>2 ) then
                  found(foundcount,1)=j
                  found(foundcount,2:NoOfNodesInCell+1)=cells(j,1:NoOfNodesInCell)
                  foundcount=foundcount+1
               endif
            enddo
         endif
      enddo

      ! check that the cells are ordered in an anti-clockwise direction.
      do L=1,NoOfNodesInCell
         do j=1,NoOfNodesInCell
            do k=1,NoOfNodesInCell
               if( (found(j,mod(k-1,NoOfNodesInCell)+2)==cells(i,mod(L,NoOfNodesInCell)+1)).and.&
               &   (found(j,mod(k,NoOfNodesInCell)+2)==cells(i,mod(L-1,NoOfNodesInCell)+1)) ) then
                  swap=found(L,:)
                  found(L,:)=found(j,:)
                  found(j,:)=swap
               endif
            enddo
         enddo
      enddo

      return
   end subroutine find_cells

   
   subroutine connect(surnodes,cells)

     implicit none

      integer,intent(in),dimension(NoOfCells,NoOfNodesInCell) :: cells
      integer,intent(inout),dimension(NoOfNodes,MaxNoOfSurroundingCells) :: surnodes
      integer,dimension(NoOfNodes) :: count
      integer :: i, j
      count=1; surnodes=0

      do i=1,NoOfCells
         do j=1,NoOfNodesInCell
            surnodes(cells(i,j),count(cells(i,j)))=i
            count(cells(i,j))=count(cells(i,j))+1
         enddo
      enddo

      return
   end subroutine connect

end module MeshClass
