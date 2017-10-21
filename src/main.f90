program Euler2D

   use MeshClass
   use SolutionClass

   implicit none
!==============================================================================

   ! Main computing variables
   integer, dimension(:, :), pointer :: Cells, Con, Surnodes, Edges
   double precision, dimension(:, :), pointer :: x, uvertex, cellsize, xdot
   double precision, dimension(:, :), pointer :: ud, un ! old, new

   ! Time variables
   double precision :: OutputTime, Time, TimeStep, CFLCondition, ReportTime

   ! Variables for writing output files
   integer :: NoOfReportSteps, CurrentReportStep = 0
   logical :: WriteSolutionFlag

   ! Computational timing variables
   integer :: System_Time_Start, System_Time_Stop, System_Time_Rate
   real :: CPU_Time_Start, CPU_Time_Stop

!  integer :: status
!==============================================================================

   ! Write program banner
   write (6, *)
   write (6, *) '---------------------------------------------------------------------'
   write (6, *) '2D Euler - Version 0.1'
   write (6, *) '---------------------------------------------------------------------'
   write (6, *)

   ! Initialise the variables
   Time = 0.0d0

   ! Start timing procedure
   call system_clock(System_Time_Start, System_Time_Rate)
   call cpu_time(CPU_Time_Start)

   ! Read in the problem data from the file 'variables.data'
   call ReadInput(OutputTime, CFLCondition, NoOfReportSteps)

   ! Read in mesh data from file
   call ReadMesh(x, xdot, Cells)

   ! Set the mesh
   call ProcessMesh(Cells, Con, Surnodes, Edges)

   ! Set the initial conditions
   call InitialCondition(ud, un, uvertex, x, Cells, Surnodes, CellSize)

   ! Write the initial condition
   call WriteSolution(uvertex, x, CurrentReportStep)
   CurrentReportStep = CurrentReportStep + 1

   ! Write the computational mesh to file
   call WriteMesh(x, Cells)

   ! Write time-stepping banner
   write (6, *)
   write (6, *) 'Advancing solution to time ', OutputTime
   write (6, *) 'Courant-Friedrichs-Lewy (CFL) value is ', CFLCondition
   write (6, *)
   write (6, *) '---------------------------------------------------------------------'
   write (6, *) '         Time                Time-step                               '
   write (6, *) '---------------------------------------------------------------------'

   do while (Time < OutputTime)

      ! Calculate the timestep
      call CalculateTimestep(TimeStep, ud, x, xdot, Cells, CFLCondition, Edges)

      ReportTime = dble(CurrentReportStep)*OutputTime/dble(NoOfReportSteps)
      if (Time + TimeStep > ReportTime) then
         TimeStep = ReportTime - Time
         WriteSolutionFlag = .true.
      endif

      ! Move the mesh
      if (MovingMeshChoice .gt. 1) then
         call MoveMesh(x, xdot, uvertex, ud, cellsize, Con, Cells, TimeStep)
      endif

      ! Solve the Euler equations
      call Solver(ud, un, x, xdot, TimeStep, Cells, cellsize, Edges)

      ! Increment the time
      Time = Time + TimeStep

      ! Write out the current time and time-step size
      write (6, '(2f20.8)') Time, TimeStep

      ! Write the solution to file
      if (WriteSolutionFlag) then
         call InterpolateSolutionToVertex(ud, uvertex, Cellsize, Surnodes)
         call WriteSolution(uvertex, x, CurrentReportStep)
         CurrentReportStep = CurrentReportStep + 1
         WriteSolutionFlag = .false.
      endif

   enddo

   ! Stop the timing
   call system_clock(System_Time_Stop)
   call cpu_time(CPU_Time_Stop)

   write (*, '(/)')
   write (*, '(1x,a)') 'Timing report:'
   write (6, *) '---------------------------------------------------------------------'
   write (*, '(1x,a)') '    Elapsed Time   CPU Time', &
      '        (s)           (s)'
   write (6, *) '---------------------------------------------------------------------'
   write (*, '(1x,2e15.4)') dble(System_Time_Stop - System_Time_Start)/dble(System_Time_Rate)&
                        &, CPU_Time_Stop - CPU_Time_Start
   write (6, *) '---------------------------------------------------------------------'
   write (*, '(/)')

end program Euler2D
