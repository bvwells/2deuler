subroutine ReadInput(OutputTime,CFLCondition,NoOfReportSteps)
  !----------------------------------------------------------------------------
  !  Routine: ReadInput
  !
  !  Description: Read input variables 
  !
  !----------------------------------------------------------------------------
  use MeshClass
  use SolutionClass
  implicit none

  ! Time variables   
  double precision :: OutputTime,CFLCondition
  ! Variables for writing output files
  integer :: NoOfReportSteps

  open(unit=10,file='variables.data',status='old')
  read(10,*) OutputTime
  read(10,*) CFLCondition
  read(10,*) MeshFilename
  read(10,*) CellFilename
  read(10,*) NoOfReportSteps
  read(10,*) MovingMeshChoice
  read(10,*) InitialConditionChoice
  close(10)

  return
end subroutine ReadInput
