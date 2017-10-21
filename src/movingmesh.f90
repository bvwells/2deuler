subroutine MoveMesh(x, xdot, uvertex, ud, cellsize, con, Cells, TimeStep)
   use MeshClass
   use SolutionClass
   implicit none

   double precision, intent(in) :: TimeStep
   double precision, intent(inout), dimension(NoOfNodes, NoOfDimensions) :: x, xdot
   double precision, intent(in), dimension(NoOfNodes, NoOfEquations) :: uvertex
   double precision, intent(in), dimension(NoOfCells, NoOfEquations) :: ud
   double precision, intent(inout), dimension(1:NoOfCells, 2) :: cellsize
   integer, intent(in), dimension(NoOfCells, NoOfNodesInCell) :: cells
   integer, intent(in), dimension(NoOfNodes, MaxNoOfSurroundingCells + 1) :: con

   if (MovingMeshChoice == 2) then
      call lagrangian(x, xdot, uvertex, cellsize, Cells, con)
      call TimeStepMesh(x, xdot, Timestep)
   elseif (MovingMeshChoice == 3) then
      call ARC(x, xdot, uvertex, cellsize, cells, con)
      call TimeStepMesh(x, xdot, Timestep)
   endif

   return
end subroutine MoveMesh

subroutine TimeStepMesh(x, xdot, Timestep)
   use MeshClass
   implicit none
   double precision, intent(in) :: TimeStep
   double precision, intent(in), dimension(NoOfNodes, NoOfDimensions) :: xdot
   double precision, intent(inout), dimension(NoOfNodes, NoOfDimensions) :: x
   integer :: i, j

   do i = 1, NoOfDimensions
      do j = 1, NoOfNodes
         x(j, i) = x(j, i) + TimeStep*xdot(j, i)
      enddo
   enddo

   return
end subroutine TimeStepMesh

subroutine ARC(x, xdot, uvertex, cellsize, cells, con)

   use MeshClass
   use SolutionClass

   use routines
   use matrix_solvers

   implicit none
!==============================================================================
   integer, intent(in), dimension(1:NoOfCells, 1:NoOfNodesInCell) :: cells
   integer, intent(in), dimension(1:NoOfNodes, 1:MaxNoOfSurroundingCells + 1) :: con
   double precision, intent(inout), dimension(1:NoOfNodes, 1:NoOfDimensions) :: x, xdot
   double precision, intent(in), dimension(1:NoOfNodes, 1:NoOfEquations) :: uvertex
   double precision, intent(in), dimension(1:NoOfCells, 1:2) :: cellsize
!==============================================================================
   double precision, dimension(1:NoOfNodes, 1:MaxNoOfSurroundingCells + 1) :: B
   double precision, dimension(1:NoOfNodes) :: psi, B_2, rhs
   double precision, dimension(1:NoOfCells, 1:NoOfDimensions, 1:NoOfNodesInCell) :: grad
   double precision, dimension(1:NoOfNodes, 1:NoOfDimensions) :: xdotstore
   double precision, dimension(1:NoOfNodesInCell, 1:NoOfNodesInCell) :: int
   double precision :: integral, alpha, beta, TOL
   integer :: i, j, k, L, loc, count
!==============================================================================
   B = 0; psi = 0; i = 0; j = 0; k = 0; L = 0; B_2 = 0; rhs = 0; grad = 0
   alpha = 0; beta = 20.0d0; TOL = 0.0000001d0; xdotstore = 0.0d0 !2.5 worked   0k
   int(1, 1) = 1.0d0/12.0d0; int(1, 2) = 1.0d0/24.0d0; int(1, 3) = 1.0d0/24.0d0
   int(2, 1) = 1.0d0/24.0d0; int(2, 2) = 1.0d0/12.0d0; int(2, 3) = 1.0d0/24.0d0
   int(3, 1) = 1.0d0/24.0d0; int(3, 2) = 1.0d0/24.0d0; int(3, 3) = 1.0d0/12.0d0

   do i = 1, NoOfCells

      call gradient(grad(i, :, :), x(cells(i, 1), 1), x(cells(i, 2), 1), x(cells(i, 3), 1),&
                               &x(cells(i, 1), 2), x(cells(i, 2), 2), x(cells(i, 3), 2), 2.0d0*cellsize(i, 1))

      alpha = dmax1(dabs((uvertex(cells(i, 1), 1)**2)*DOT_PRODUCT(grad(i, :, 1), grad(i, :, 1)) + &
           &2.0d0*(uvertex(cells(i, 1), 1)*uvertex(cells(i, 2), 1))*DOT_PRODUCT(grad(i, :, 1), grad(i, :, 2)) + &
           &2.0d0*(uvertex(cells(i, 1), 1)*uvertex(cells(i, 3), 1))*DOT_PRODUCT(grad(i, :, 1), grad(i, :, 3)) + &
           &(uvertex(cells(i, 2), 1)**2)*DOT_PRODUCT(grad(i, :, 2), grad(i, :, 2)) + &
           &2.0d0*(uvertex(cells(i, 2), 1)*uvertex(cells(i, 3), 1))*DOT_PRODUCT(grad(i, :, 2), grad(i, :, 3)) + &
           &(uvertex(cells(i, 3), 1)**2)*DOT_PRODUCT(grad(i, :, 3), grad(i, :, 3))), alpha)
   enddo

   if (dabs(alpha) <= TOL) then
      alpha = 0.0d0
   else
      alpha = beta/alpha
   endif

   do i = 1, NoOfCells
      integral = SQRT(1.0d0 + alpha*((uvertex(cells(i, 1), 1)**2)*DOT_PRODUCT(grad(i, :, 1), grad(i, :, 1)) + &
           &2.0d0*(uvertex(cells(i, 1), 1)*uvertex(cells(i, 2), 1))*DOT_PRODUCT(grad(i, :, 1), grad(i, :, 2)) + &
           &2.0d0*(uvertex(cells(i, 1), 1)*uvertex(cells(i, 3), 1))*DOT_PRODUCT(grad(i, :, 1), grad(i, :, 3)) + &
           &(uvertex(cells(i, 2), 1)**2)*DOT_PRODUCT(grad(i, :, 2), grad(i, :, 2)) + &
           &2.0d0*(uvertex(cells(i, 2), 1)*uvertex(cells(i, 3), 1))*DOT_PRODUCT(grad(i, :, 2), grad(i, :, 3)) + &
           &(uvertex(cells(i, 3), 1)**2)*DOT_PRODUCT(grad(i, :, 3), grad(i, :, 3))))

      do j = 1, NoOfNodesInCell
         do k = 1, NoOfNodesInCell
            loc = 0
            do L = 1, MaxNoOfSurroundingCells + 1
               if (cells(i, k) == con(cells(i, j), L)) then
                  loc = L
                  exit
               endif
            enddo
            B(cells(i, j), loc) = B(cells(i, j), loc) + cellsize(i, 1)*integral*DOT_PRODUCT(grad(i, :, j), grad(i, :, k))
         enddo

         rhs(cells(i, j)) = rhs(cells(i, j)) - alpha*cellsize(i, 1)*&
              &(uvertex(cells(i, 1), 1)*DOT_PRODUCT(grad(i, :, 1), grad(i, :, j)) +&
              &  uvertex(cells(i, 2), 1)*DOT_PRODUCT(grad(i, :, 2), grad(i, :, j)) +&
              &  uvertex(cells(i, 3), 1)*DOT_PRODUCT(grad(i, :, 3), grad(i, :, j)))*&
              & (uvertex(cells(i, 1), 2)*grad(i, 1, 1) + uvertex(cells(i, 1), 3)*grad(i, 2, 1) +&
              &  uvertex(cells(i, 2), 2)*grad(i, 1, 2) + uvertex(cells(i, 2), 3)*grad(i, 2, 2) +&
              &  uvertex(cells(i, 3), 2)*grad(i, 1, 3) + uvertex(cells(i, 3), 3)*grad(i, 2, 3))/integral
      enddo
   enddo

   call sparse_cg(B, Con, psi, rhs, NoOfNodes, MaxNoOfSurroundingCells, 1)
   B = 0; rhs = 0

   do i = 1, NoOfCells
      do j = 1, NoOfNodesInCell
         do k = 1, NoOfNodesInCell
            loc = 0
            do L = 1, MaxNoOfSurroundingCells + 1
               if (cells(i, k) == con(cells(i, j), L)) then
                  loc = L
                  exit
               endif
            enddo
            B(cells(i, j), loc) = B(cells(i, j), loc) + 2.0d0*cellsize(i, 1)*int(j, k)
         enddo

         rhs(cells(i, j)) = rhs(cells(i, j)) + 2.0d0*cellsize(i, 1)*(psi(cells(i, 1))*grad(i, 1, 1) + &
              &psi(cells(i, 2))*grad(i, 1, 2) + psi(cells(i, 3))*grad(i, 1, 3))/6.0d0
         B_2(cells(i, j)) = B_2(cells(i, j)) + 2.0d0*cellsize(i, 1)*(psi(cells(i, 1))*grad(i, 2, 1) + &
              &psi(cells(i, 2))*grad(i, 2, 2) + psi(cells(i, 3))*grad(i, 2, 3))/6.0d0

      enddo
   enddo

   call sparse_cg(B, Con, xdot(:, 1), rhs, NoOfNodes, MaxNoOfSurroundingCells, 2)
   call sparse_cg(B, Con, xdot(:, 2), B_2, NoOfNodes, MaxNoOfSurroundingCells, 2)

   do k = 1, 4
      xdotstore(:, :) = xdot(:, :)
      xdot(:, :) = 0.0d0
      do i = 1, NoOfNodes
         count = 0
         do j = 1, MaxNoOfSurroundingCells + 1
            if (con(i, j) /= 0) then
               xdot(i, :) = xdot(i, :) + xdotstore(con(i, j), :)
               count = count + 1
            endif
         enddo
         xdot(i, :) = xdot(i, :)/abs(count)
      enddo
   enddo

   return

end subroutine ARC

subroutine lagrangian(x, xdot, uvertex, cellsize, cells, Con)

   use MeshClass
   use SolutionClass

   use routines
   use matrix_solvers

   implicit none
!==============================================================================
   integer, intent(in), dimension(1:NoOfCells, 1:NoOfNodesInCell) :: cells
   integer, intent(in), dimension(1:NoOfNodes, 1:MaxNoOfSurroundingCells + 1) :: con
   double precision, intent(inout), dimension(1:NoOfNodes, 1:NoOfDimensions) :: x, xdot
   double precision, intent(in), dimension(1:NoOfCells, 1:2) :: cellsize
   double precision, intent(in), dimension(1:NoOfNodes, 1:NoOfEquations) :: uvertex
!==============================================================================
   double precision, dimension(1:NoOfNodes, 1:MaxNoOfSurroundingCells + 1) :: B, mass
   double precision, dimension(1:NoOfNodes) :: psi, B_1, B_2, rhs
   double precision, dimension(1:NoOfCells, 1:NoOfDimensions, 1:NoOfNodesInCell) :: grad
   double precision, dimension(1:NoOfNodesInCell, 1:NoOfNodesInCell) :: int
   double precision :: integral
   integer :: i, j, k, L, loc, count
!==============================================================================
   B = 0; psi = 0; i = 0; j = 0; k = 0; L = 0; mass = 0; B_1 = 0; B_2 = 0; rhs = 0; grad = 0; 
   int(1, 1) = 1.0d0/12.0d0; int(1, 2) = 1.0d0/24.0d0; int(1, 3) = 1.0d0/24.0d0
   int(2, 1) = 1.0d0/24.0d0; int(2, 2) = 1.0d0/12.0d0; int(2, 3) = 1.0d0/24.0d0
   int(3, 1) = 1.0d0/24.0d0; int(3, 2) = 1.0d0/24.0d0; int(3, 3) = 1.0d0/12.0d0

   do i = 1, NoOfCells

      call gradient(grad(i, :, :), x(cells(i, 1), 1), x(cells(i, 2), 1), x(cells(i, 3), 1),&
                               &x(cells(i, 1), 2), x(cells(i, 2), 2), x(cells(i, 3), 2), 2.0d0*cellsize(i, 1))

      integral = cellsize(i, 1)*(uvertex(cells(i, 1), 1) + uvertex(cells(i, 2), 1) + uvertex(cells(i, 3), 1))/3.0d0

      do j = 1, NoOfNodesInCell
         do k = 1, NoOfNodesInCell
            loc = 0
            do L = 1, MaxNoOfSurroundingCells + 1
               if (cells(i, k) == con(cells(i, j), L)) then
                  loc = L
                  exit
               endif
            enddo
            B(cells(i, j), loc) = B(cells(i, j), loc) + integral*DOT_PRODUCT(grad(i, :, j), grad(i, :, k))
            mass(cells(i, j), loc) = mass(cells(i, j), loc) + 2.0d0*cellsize(i, 1)*int(j, k)
         enddo

         rhs(cells(i, j)) = rhs(cells(i, j)) - cellsize(i, 1)*&
              &(grad(i, 1, j)*(uvertex(cells(i, 1), 2) + uvertex(cells(i, 2), 2) + uvertex(cells(i, 3), 2)) +&
              &grad(i, 2, j)*(uvertex(cells(i, 1), 3) + uvertex(cells(i, 2), 3) + uvertex(cells(i, 3), 3)))/3.0d0
      enddo
   enddo

   call sparse_cg(B, Con, psi, rhs, NoOfNodes, MaxNoOfSurroundingCells, 1)

   do i = 1, NoOfCells
      do j = 1, NoOfNodesInCell
         B_1(cells(i, j)) = B_1(cells(i, j)) + cellsize(i, 1)*(psi(cells(i, 1))*grad(i, 1, 1) +  &
              &psi(cells(i, 2))*grad(i, 1, 2) + psi(cells(i, 3))*grad(i, 1, 3))/3.0d0
         B_2(cells(i, j)) = B_2(cells(i, j)) + cellsize(i, 1)*(psi(cells(i, 1))*grad(i, 2, 1) +  &
              &psi(cells(i, 2))*grad(i, 2, 2) + psi(cells(i, 3))*grad(i, 2, 3))/3.0d0
      enddo
   enddo

   call sparse_cg(mass, Con, xdot(:, 1), B_1, NoOfNodes, MaxNoOfSurroundingCells, 2)
   call sparse_cg(mass, Con, xdot(:, 2), B_2, NoOfNodes, MaxNoOfSurroundingCells, 2)

   return
end subroutine lagrangian
