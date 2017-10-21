module SolutionClass
   implicit none

   integer, parameter :: NoOfEquations = 4
   integer :: InitialConditionChoice
   double precision, parameter :: MinimumTimeStep = 1d-6
   double precision, parameter :: Gamma = 1.4d0

contains

   subroutine InitialCondition(ud, un, uvertex, x, cells, surnodes, cellsize)
      use MeshClass
      implicit none
      !==============================================================================
      integer, intent(inout), dimension(NoOfCells, NoOfNodesInCell) :: cells
      integer, intent(in), dimension(NoOfNodes, MaxNoOfSurroundingCells) :: surnodes
      double precision, dimension(:, :), pointer :: ud, un, uvertex, CellSize
      double precision, intent(inout), dimension(NoOfNodes, NoOfDimensions) :: x
      !==============================================================================
      double precision, dimension(1:NoOfDimensions) :: xcent, xcentre
      double precision :: r
      integer :: i
      !==============================================================================
      xcentre = 0.0d0

      ! Allocate memory
      allocate (ud(NoOfCells, NoOfEquations))
      allocate (un(NoOfCells, NoOfEquations))
      allocate (CellSize(NoOfCells, 2))
      allocate (uvertex(NoOfNodes, NoOfEquations))

      write (6, *) 'Calculating initial condition...'
      if (InitialConditionChoice == 1) then
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            r = CalculateLength(xcent, xcentre, NoOfDimensions)
            if (r < 0.5d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 1.0d0/(gamma - 1.0d0)
            else
               ud(i, 1) = 0.125d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 0.1d0/(gamma - 1.0d0)
            endif
            CellSize(i, 1) = CalculateCellArea(Cells, x, i)
            CellSize(i, 2) = CellSize(i, 1)
         enddo
      elseif (InitialConditionChoice == 2) then
         ! Sedov blast wave
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            r = CalculateLength(xcent, xcentre, NoOfDimensions)
            if (r < 0.05d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 8.0d0
            else
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 0.0d0
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      elseif (InitialConditionChoice == 3) then
         ! Sedov blast wave
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            r = CalculateLength(xcent, xcentre, NoOfDimensions)
            if (r < 0.4d0) then
               ud(i, 1) = 10.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 1000.0d0/(gamma - 1.0d0)
            else
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 0.1d0/(gamma - 1.0d0)
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      elseif (InitialConditionChoice == 4) then
         ! Shock tube problem
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            if (xcent(1) < 0.5d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 1.0d0/(gamma - 1.0d0)
            else
               ud(i, 1) = 0.125d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 0.1d0/(gamma - 1.0d0)
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      elseif (InitialConditionChoice == 5) then
         ! Blast wave problem
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            if (xcent(1) <= -0.80d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 1000.0d0/(gamma - 1.0d0)
            elseif (xcent(1) > -0.8d0 .and. xcent(1) < 0.8d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 100.0d0/(gamma - 1.0d0)
            else
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 0.01d0/(gamma - 1.0d0)
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      elseif (InitialConditionChoice == 6) then
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            if (xcent(1) <= 0.0d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = 0.75d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 1.0d0/(gamma - 1.0d0) + 0.5d0*ud(i, 1)*((ud(i, 2)/ud(i, 1))**2 + (ud(i, 3)/ud(i, 1))**2)
            else
               ud(i, 1) = 0.125d0
               ud(i, 2) = 0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = 0.1d0/(gamma - 1.0d0) + 0.5d0*ud(i, 1)*((ud(i, 2)/ud(i, 1))**2 + (ud(i, 3)/ud(i, 1))**2)
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      elseif (InitialConditionChoice == 7) then
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            if (xcent(1) <= 0.0d0) then
               ud(i, 1) = 1.0d0
               ud(i, 2) = -2.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = (0.4d0/(gamma - 1.0d0)) + 0.5d0*ud(i, 1)*((ud(i, 2)/ud(i, 1))**2 + (ud(i, 3)/ud(i, 1))**2)
            else
               ud(i, 1) = 1.0d0
               ud(i, 2) = 2.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = (0.4d0/(gamma - 1.0d0)) + 0.5d0*ud(i, 1)*((ud(i, 2)/ud(i, 1))**2 + (ud(i, 3)/ud(i, 1))**2)
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      elseif (InitialConditionChoice == 8) then
         do i = 1, NoOfCells
            xcent = CalculateCellCentre(Cells, x, i)
            if (xcent(1) <= 0.0d0) then
               ud(i, 1) = 0.445d0
               ud(i, 2) = ud(i, 1)*0.698d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = (3.528d0/(gamma - 1.0d0)) + 0.5d0*ud(i, 1)*((ud(i, 2)/ud(i, 1))**2 + (ud(i, 3)/ud(i, 1))**2)
            else
               ud(i, 1) = 0.5d0
               ud(i, 2) = ud(i, 1)*0.0d0
               ud(i, 3) = 0.0d0
               ud(i, 4) = (0.571d0/(gamma - 1.0d0)) + 0.5d0*ud(i, 1)*((ud(i, 2)/ud(i, 1))**2 + (ud(i, 3)/ud(i, 1))**2)
            endif
            cellsize(i, 1) = CalculateCellArea(Cells, x, i)
            cellsize(i, 2) = cellsize(i, 1)
         enddo
      endif

      call InterpolateSolutionToVertex(ud, uvertex, Cellsize, surnodes)

      return
   end subroutine InitialCondition

   subroutine WriteSolution(u, x, CurrentReportStep)
      use MeshClass
      implicit none
      integer, intent(in) :: CurrentReportStep
      double precision, intent(in) :: u(NoOfNodes, NoOfEquations)
      double precision, intent(in) :: x(NoOfNodes, NoOfDimensions)
      ! Variables for writing output files
      integer :: test_number, hundreds, tens, units, i, j
      character(LEN=32) :: filename
      character(LEN=10) :: numbers
      ! Format statements
      character(LEN=16) :: F
      parameter(F='(6F16.6)')

      numbers = "0123456789"; filename = "solution"

      test_number = CurrentReportStep
      hundreds = test_number/100
      test_number = test_number - 100*hundreds
      tens = test_number/10
      test_number = test_number - 10*tens
      units = test_number

      open (unit=10, file=trim(filename)//numbers(hundreds + 1:hundreds + 1)//&
                                       &numbers(tens + 1:tens + 1)//&
                                       &numbers(units + 1:units + 1)//".m")
      do i = 1, NoOfNodes
         write (10, F) (x(i, j), j=1, NoOfDimensions), (u(i, j), j=1, NoOfEquations)
      enddo
      close (10)

      return
   end subroutine WriteSolution

   subroutine InterpolateSolutionToVertex(u, uvertex, cellsize, surnodes)
      use MeshClass
      implicit none
      integer, intent(in), dimension(NoOfNodes, MaxNoOfSurroundingCells) :: surnodes
      double precision, intent(inout), dimension(NoOfNodes, NoOfEquations) :: uvertex
      double precision, intent(in), dimension(NoOfCells, NoOfEquations) :: u
      double precision, intent(in), dimension(NoOfCells, 2) :: cellsize
      double precision :: total_area
      integer :: i, j, k

      do i = 1, NoOfNodes
         total_area = 0.0d0
         do j = 1, NoOfEquations
            uvertex(i, j) = 0.0d0
         enddo
         do j = 1, MaxNoOfSurroundingCells
            if (surnodes(i, j) /= 0) then
               do k = 1, NoOfEquations
                  uvertex(i, k) = uvertex(i, k) + cellsize(surnodes(i, j), 1)*u(surnodes(i, j), k)
               enddo
               total_area = total_area + cellsize(surnodes(i, j), 1)
            endif
         enddo
         do j = 1, NoOfEquations
            uvertex(i, j) = uvertex(i, j)/total_area
         enddo
      enddo

      return
   end subroutine InterpolateSolutionToVertex

   function Pressure(u) result(P)
      implicit none
      double precision :: u(NoOfEquations)
      double precision :: P

      P = (gamma - 1.0d0)*(u(4) - 0.5d0*u(1)*((u(2)/u(1))**2 + (u(3)/u(1))**2))

      return
   end function Pressure

   subroutine CalculateEigenvalues(EigenValues, RoeAverage, SoundSpeed, ReferenceSpeed)
      implicit none
      double precision, intent(out) :: EigenValues(NoOfEquations)
      double precision, intent(in)  :: RoeAverage(NoOfEquations)
      double precision, intent(in)  :: SoundSpeed, ReferenceSpeed

      EigenValues(1) = RoeAverage(2) - SoundSpeed - ReferenceSpeed
      EigenValues(2) = RoeAverage(2) - ReferenceSpeed
      EigenValues(3) = RoeAverage(2) - ReferenceSpeed
      EigenValues(4) = RoeAverage(2) + SoundSpeed - ReferenceSpeed

      return
   end subroutine CalculateEigenvalues

   subroutine CalculateEigenVectors(EigenVectors, RoeAverage, SoundSpeed)
      implicit none
      double precision, intent(out) :: EigenVectors(NoOfEquations, NoOfEquations)
      double precision, intent(in)  :: RoeAverage(NoOfEquations), SoundSpeed

      EigenVectors(1, 1) = 1.0d0
      EigenVectors(1, 2) = 1.0d0
      EigenVectors(1, 3) = 0.0d0
      EigenVectors(1, 4) = 1.0d0

      EigenVectors(2, 1) = RoeAverage(2) - SoundSpeed
      EigenVectors(2, 2) = RoeAverage(2)
      EigenVectors(2, 3) = 0.0d0
      EigenVectors(2, 4) = RoeAverage(2) + SoundSpeed

      EigenVectors(3, 1) = RoeAverage(3)
      EigenVectors(3, 2) = RoeAverage(3)
      EigenVectors(3, 3) = 1.0d0
      EigenVectors(3, 4) = RoeAverage(3)

      EigenVectors(4, 1) = RoeAverage(4) - RoeAverage(2)*SoundSpeed
      EigenVectors(4, 2) = 0.5d0*(RoeAverage(2)**2 + RoeAverage(3)**2)
      EigenVectors(4, 3) = RoeAverage(3)
      EigenVectors(4, 4) = RoeAverage(4) + RoeAverage(2)*SoundSpeed

      return
   end subroutine CalculateEigenVectors

   subroutine CalculateRoeAverages(RoeAverage, u_l, u_r)
      implicit none
      double precision, intent(out) :: RoeAverage(NoOfEquations)
      double precision, intent(in)  :: u_l(NoOfEquations)
      double precision, intent(in)  :: u_r(NoOfEquations)
      double precision :: rootrho_l, rootrho_r, denom, h_l, h_r

      rootrho_l = sqrt(u_l(1))
      rootrho_r = sqrt(u_r(1))

      RoeAverage(1) = rootrho_l*rootrho_r

      denom = rootrho_l + rootrho_r
      RoeAverage(2) = (rootrho_l*(u_l(2)/u_l(1)) + rootrho_r*(u_r(2)/u_r(1)))/denom
      RoeAverage(3) = (rootrho_l*(u_l(3)/u_l(1)) + rootrho_r*(u_r(3)/u_r(1)))/denom

      H_l = (u_l(4) + Pressure(u_l))/u_l(1)
      H_r = (u_r(4) + Pressure(u_r))/u_r(1)

      RoeAverage(4) = (rootrho_l*H_l + rootrho_r*H_r)/denom

      return
   end subroutine CalculateRoeAverages

   subroutine CalculateSoundSpeed(SoundSpeed, RoeAverage)
      implicit none
      double precision, intent(out) :: SoundSpeed
      double precision, intent(in)  :: RoeAverage(NoOfEquations)

      SoundSpeed = sqrt((gamma - 1.0d0)*(RoeAverage(4) - 0.5d0*(RoeAverage(2)**2 + RoeAverage(3)**2)))

      return
   end subroutine CalculateSoundSpeed

   subroutine CalculateWaveStrengths(WaveStrengths, DeltaU, Eigenvectors, SoundSpeed, RoeAverage)
      implicit none

      double precision, intent(out) :: WaveStrengths(NoOfEquations)
      double precision, intent(in)  :: EigenVectors(NoOfEquations, NoOfEquations)
      double precision, intent(in)  :: DeltaU(NoOfEquations)
      double precision, intent(in)  :: RoeAverage(NoOfEquations), SoundSpeed

      ! Solution to the matrix system
      ! DeltaU = EigenVectors * WaveStrengths
      ! Note: Could be made more generic by replacing with 4x4 matrix solver

      WaveStrengths(3) = DeltaU(3) - RoeAverage(3)*DeltaU(1)
      WaveStrengths(2) = ((gamma - 1.0d0)/(SoundSpeed**2))*(DeltaU(1)*(RoeAverage(4) - (RoeAverage(2)**2)) +&
           & RoeAverage(2)*DeltaU(2) - DeltaU(4) + (DeltaU(3) - RoeAverage(3)*DeltaU(1))*RoeAverage(3))
      WaveStrengths(1) = (1.0d0/(2.0d0*SoundSpeed))*(DeltaU(1)*(RoeAverage(2) + SoundSpeed) - DeltaU(2) -&
           & SoundSpeed*WaveStrengths(2))
      WaveStrengths(4) = DeltaU(1) - (WaveStrengths(1) + WaveStrengths(2))

!      WaveStrengths=0.0d0
!      call SetLinearSolverParameters(1,1,0,50,1d-10)
!      call SolveLinearSystem(EigenVectors,WaveStrengths,DeltaU,NoOfEquations)

      return
   end subroutine CalculateWaveStrengths

   subroutine FluxFunction(Flux, U, ReferenceSpeed)
      implicit none
      double precision, intent(out) :: Flux(NoOfEquations)
      double precision, intent(in)  :: U(NoOfEquations), ReferenceSpeed

      Flux(1) = U(2) - ReferenceSpeed*U(1)
      Flux(2) = U(2)*(U(2)/U(1)) + Pressure(U) - ReferenceSpeed*U(2)
      Flux(3) = U(2)*(U(3)/U(1)) - ReferenceSpeed*U(3)
      Flux(4) = (U(2)/U(1))*(U(4) + Pressure(U)) - ReferenceSpeed*U(4)

      return
   end subroutine FluxFunction

   subroutine RoeRiemannSolver(Flux, u_l, u_r, ReferenceSpeed)
      implicit none
      double precision, intent(inout), dimension(NoOfEquations) :: Flux
      double precision, intent(in), dimension(NoOfEquations) :: u_l, u_r
      double precision, intent(in) :: ReferenceSpeed
      double precision, dimension(NoOfEquations, NoOfEquations) :: EigenVectors
      double precision, dimension(NoOfEquations) :: EigenValues, WaveStrengths, DeltaU
      double precision, dimension(NoOfEquations) :: FluxL, FluxR, RoeAverage
      double precision :: SoundSpeed
      integer :: i, j

      call CalculateRoeAverages(RoeAverage, u_l, u_r)

      call CalculateSoundSpeed(SoundSpeed, RoeAverage)

      call CalculateEigenvalues(EigenValues, RoeAverage, SoundSpeed, ReferenceSpeed)

      call CalculateEigenVectors(EigenVectors, RoeAverage, SoundSpeed)

      DeltaU = u_r - u_l

      call CalculateWaveStrengths(WaveStrengths, DeltaU, Eigenvectors, SoundSpeed, RoeAverage)

      call FluxFunction(FluxL, u_l, ReferenceSpeed)
      call FluxFunction(FluxR, u_r, ReferenceSpeed)

      Flux = 0.5d0*(FluxL + FluxR)
      do i = 1, NoOfEquations
         do j = 1, NoOfEquations
            Flux(j) = Flux(j) - 0.5d0*WaveStrengths(i)*dabs(EigenValues(i))*EigenVectors(j, i)
         enddo
      enddo

      return

   end subroutine RoeRiemannSolver

   subroutine solver(ud, un, x, xdot, TimeStep, cells, cellsize, Edges)

      use MeshClass

      implicit none
      !==============================================================================
      integer, intent(in), dimension(NoOfCells, NoOfNodesInCell) :: cells
      integer, intent(in), dimension(NoOfEdges, 4) :: Edges
      double precision, intent(inout), dimension(NoOfCells, NoOfEquations) :: ud, un
      double precision, intent(in), dimension(NoOfNodes, NoOfDimensions) :: x, xdot
      double precision, intent(inout), dimension(NoOfCells, 1:2) :: cellsize
      double precision, intent(in) :: TimeStep
      !==============================================================================
      double precision, dimension(NoOfEquations) :: Flux, ul_edge, ur_edge
      double precision :: vol1, vol2, Factor, Speed, r
      integer :: Node1, Node2, Cell1, Cell2
      integer :: i, j, k
      !==============================================================================

      ! Mass accumulation terms
      do i = 1, NoOfCells
         cellsize(i, 2) = CalculateCellArea(Cells, x, i)
         Factor = cellsize(i, 1)/cellsize(i, 2)
         do j = 1, NoOfEquations
            un(i, j) = Factor*ud(i, j)
         enddo
      enddo

      ! Flux terms
      do i = 1, NoOfEdges

         Node1 = Edges(i, 1)
         Node2 = Edges(i, 2)
         Cell1 = Edges(i, 3)
         Cell2 = Edges(i, 4)

         if (Node1 .eq. 0 .or. Node2 .eq. 0) cycle
         if (Cell1 .eq. 0) Cell1 = Cell2
         if (Cell2 .eq. 0) Cell2 = Cell1

         ! Interpolate the solution to the common edge
         do j = 1, NoOfEquations
            ul_edge(j) = ud(Cell1, j)
            ur_edge(j) = ud(Cell2, j)
            do k = 1, NoOfDimensions
               ul_edge(j) = ul_edge(j)
               ur_edge(j) = ur_edge(j)
            enddo
         enddo

         ! Calculate the speed of the moving edge
         vol1 = 0.5d0*((x(Node2, 1) - x(Node1, 1) - TimeStep*(xdot(Node2, 1) - xdot(Node1, 1)))*xdot(Node1, 2) - &
                       (x(Node2, 2) - x(Node1, 2) - TimeStep*(xdot(Node2, 2) - xdot(Node1, 2)))*xdot(Node1, 1))
         vol2 = -0.5d0*((x(Node2, 2) - x(Node1, 2) - TimeStep*xdot(Node2, 2))*xdot(Node2, 1) - &
                        (x(Node2, 1) - x(Node1, 1) - TimeStep*xdot(Node2, 1))*xdot(Node2, 2))
         r = CalculateLength(x(Node1, :), x(Node2, :), NoOfDimensions)
         Speed = -(vol1 + vol2)/r

         ! Calculate the Euler flux
         call CalculateEulerFlux(Flux, ul_edge, ur_edge, x, Node1, Node2, Speed)

         ! Update the solution
         do j = 1, NoOfEquations
            if (Edges(i, 3) .gt. 0) un(Cell1, j) = un(Cell1, j) - (TimeStep/cellsize(Cell1, 2))*Flux(j)
            if (Edges(i, 4) .gt. 0) un(Cell2, j) = un(Cell2, j) + (TimeStep/cellsize(Cell2, 2))*Flux(j)
         enddo

      enddo

      ! Re-set the solution
      call ResetSolution(ud, un, Cellsize)

      return

   end subroutine solver

   subroutine ResetSolution(ud, un, CellSize)
      use MeshClass
      implicit none
      double precision, intent(inout) :: ud(NoOfCells, NoOfEquations)
      double precision, intent(inout) :: un(NoOfCells, NoOfEquations)
      double precision, intent(inout) :: cellsize(NoOfCells, NoOfEquations)
      integer :: i, j

      ! Re-set the solution
      do i = 1, NoOfCells
         cellsize(i, 1) = cellsize(i, 2)
         do j = 1, NoOfEquations
            ud(i, j) = un(i, j)
         enddo
      enddo

      return
   end subroutine ResetSolution

   subroutine CalculateEulerFlux(Flux, ul_edge, ur_edge, x, Node1, Node2, Speed)
      use MeshClass
      implicit none

      integer, intent(in) :: Node1, Node2
      double precision, intent(out) :: Flux(NoOfEquations)
      double precision, intent(in) :: ul_edge(NoOfEquations), ur_edge(NoOfEquations)
      double precision, intent(in) :: x(NoOfNodes, NoOfDimensions), Speed

      double precision :: r, normal(NoOfDimensions), u_l(NoOfEquations), u_r(NoOfEquations), FluxRot(NoOfEquations)

      r = CalculateLength(x(Node1, :), x(Node2, :), NoOfDimensions)

      ! Calculate the outward normals to the common edge.
      normal(1) = (x(Node1, 2) - x(Node2, 2))/r
      normal(2) = (x(Node2, 1) - x(Node1, 1))/r

      ! Rotate the riemann data states ready for calculation.
      u_l(1) = ul_edge(1)
      u_l(2) = Normal(1)*ul_edge(2) + Normal(2)*ul_edge(3)
      u_l(3) = -Normal(2)*ul_edge(2) + Normal(1)*ul_edge(3)
      u_l(4) = ul_edge(4)

      u_r(1) = ur_edge(1)
      u_r(2) = Normal(1)*ur_edge(2) + Normal(2)*ur_edge(3)
      u_r(3) = -Normal(2)*ur_edge(2) + Normal(1)*ur_edge(3)
      u_r(4) = ur_edge(4)

      call hllc_ale(FluxRot, u_l, u_r, Speed)
      !call RoeRiemannSolver(FluxRot,u_l,u_r,Speed)

      Flux(1) = r*Fluxrot(1)
      Flux(2) = r*(Normal(1)*Fluxrot(2) - Normal(2)*Fluxrot(3))
      Flux(3) = r*(Normal(2)*Fluxrot(2) + Normal(1)*Fluxrot(3))
      Flux(4) = r*Fluxrot(4)

      return
   end subroutine CalculateEulerFlux

   subroutine CalculateTimestep(Timestep, u, x, xdot, cells, CFLCondition, Edges)

      use MeshClass

      implicit none
      !==============================================================================
      integer, intent(in), dimension(NoOfCells, NoOfNodesInCell) :: cells
      integer, intent(in), dimension(NoOfEdges, 4) :: Edges
      double precision, intent(in), dimension(NoOfCells, NoOfEquations) :: u
      double precision, intent(in), dimension(NoOfNodes, NoOfDimensions) :: x, xdot
      double precision, intent(inout) :: Timestep
      double precision, intent(in) :: CFLCondition
      !==============================================================================
      double precision, dimension(NoOfDimensions) :: xmidl, xmidr, xdotl, xdotr
      double precision :: fluid_speed, mesh_speed, speed, dx, al, ar, pl, pr
      integer :: Cell1, Cell2, Node1, Node2
      integer :: i, j
      !==============================================================================
      speed = 0.0d0

      do i = 1, NoOfEdges

         Node1 = Edges(i, 1)
         Node2 = Edges(i, 2)
         Cell1 = Edges(i, 3)
         Cell2 = Edges(i, 4)

         if (Node1 .eq. 0 .or. Node2 .eq. 0) cycle
         if (Cell1 .eq. 0) cycle
         if (Cell2 .eq. 0) cycle

         do j = 1, NoOfDimensions
            xmidl(j) = (x(cells(Cell1, 1), j) + x(cells(Cell1, 2), j) + x(cells(Cell1, 3), j))/3.0d0
            xmidr(j) = (x(cells(Cell2, 1), j) + x(cells(Cell2, 2), j) + x(cells(Cell2, 3), j))/3.0d0
         enddo
         dx = CalculateLength(xmidl, xmidr, NoOfDimensions)

         pl = Pressure(u(Cell1, :))
         pr = Pressure(u(Cell2, :))

         al = sqrt(gamma*pl/u(Cell1, 1))
         ar = sqrt(gamma*pr/u(Cell2, 1))

         xdotl(:) = (xdot(cells(Cell1, 1), :) + xdot(cells(Cell1, 2), :) + xdot(cells(Cell1, 3), :))/3.0d0
         xdotr(:) = (xdot(cells(Cell2, 1), :) + xdot(cells(Cell2, 2), :) + xdot(cells(Cell2, 3), :))/3.0d0

         mesh_speed = sqrt(dabs(xdotl(1)**2 + xdotl(2)**2))
         fluid_speed = sqrt(dabs((u(Cell1, 2)/u(Cell1, 1))**2 + (u(Cell1, 3)/u(Cell1, 1))**2))
         speed = dmax1(dabs(fluid_speed + al - mesh_speed)/dx, speed)

         mesh_speed = sqrt(dabs(xdotr(1)**2 + xdotr(2)**2))
         fluid_speed = sqrt(dabs((u(Cell2, 2)/u(Cell2, 1))**2 + (u(Cell2, 3)/u(Cell2, 1))**2))
         speed = dmax1(dabs(fluid_speed + ar - mesh_speed)/dx, speed)
      enddo

      Timestep = max(CFLCondition/speed, MinimumTimeStep)

      return
   end subroutine CalculateTimestep

end module SolutionClass
