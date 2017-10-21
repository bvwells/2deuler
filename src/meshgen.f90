program MeshGenerator

! This code generates a uniform mesh on the region [xmin(1),xmax(1)]*[xmin(2),xmax(2)]
   implicit none

   integer :: n(2), i
   double precision :: xmin(2), xmax(2), dx(2), time1, time2
   logical qtriangle

   ! Start timing
   call cpu_time(time1)

   ! Set the parameters for the mesh
   n(1) = 100
   n(2) = 100
   xmin(1) = -1.0d0
   xmin(2) = -1.0d0
   xmax(1) = 1.0d0
   xmax(2) = 1.0d0
   do i = 1, 2
      dx(i) = (xmax(i) - xmin(i))/dble(n(i) - 1)
   enddo

   qtriangle = .true.
   if (qtriangle) then
      ! Generate mesh with triangular elements
      call Triangles(n, xmin, xmax, dx)
   else
      ! Generate mesh with quadrilateral elements
      call Quads(n, xmin, xmax, dx)
   endif

   ! Finish timing
   call cpu_time(time2)

   write (6, *) 'Total time of computation is ', time2 - time1

end program MeshGenerator

subroutine Triangles(n, xmin, xmax, dx)
   implicit none
   integer, intent(in) :: n(2)
   double precision, intent(in) :: xmin(2), xmax(2), dx(2)
   double precision :: x(2)
   integer :: tri(1:4, 1:3), nnodes, ntris, reg, i, j, k, l

   write (6, *) 'Generating mesh with triangular elements'
   write (6, *) 'Mesh dimensions:'
   write (6, *) '  xmin=', xmin(1)
   write (6, *) '  ymin=', xmin(2)
   write (6, *) '  xmax=', xmax(1)
   write (6, *) '  ymax=', xmax(2)
   write (6, *) '  dx  =', dx(1)
   write (6, *) '  dy  =', dx(2)

   nnodes = (n(1) + 1)*(n(2) + 1)
   ntris = 2*n(1)*n(2)
   reg = n(1)/2

   open (unit=10, file='square.node')
   write (10, *) nnodes

   open (unit=20, file='square.cell')
   write (20, *) ntris, 3

   do i = 1, n(2)
      k = (2*nint(i/2.0d0))
      do j = 1, reg
         if (k == i) then !even
            tri(1, 1) = 2*j - 1 + (i - 1)*(n(1) + 1)
            tri(1, 2) = tri(1, 1) + 1
            tri(1, 3) = tri(1, 2) + n(1) + 1
            tri(2, 1) = tri(1, 1)
            tri(2, 2) = tri(1, 3)
            tri(2, 3) = tri(1, 3) - 1
            tri(3, 1) = tri(1, 2)
            tri(3, 2) = tri(1, 2) + 1
            tri(3, 3) = tri(1, 3)
            tri(4, 1) = tri(3, 2)
            tri(4, 2) = tri(1, 3) + 1
            tri(4, 3) = tri(1, 3)
         else
            tri(1, 1) = 2*j - 1 + (i - 1)*(n(1) + 1)
            tri(1, 2) = tri(1, 1) + 1
            tri(1, 3) = tri(1, 1) + n(1) + 1
            tri(2, 1) = tri(1, 2)
            tri(2, 2) = tri(1, 3) + 1
            tri(2, 3) = tri(1, 3)
            tri(3, 1) = tri(1, 2)
            tri(3, 2) = tri(1, 2) + 1
            tri(3, 3) = tri(2, 2) + 1
            tri(4, 1) = tri(1, 2)
            tri(4, 2) = tri(3, 3)
            tri(4, 3) = tri(2, 2)
         endif
         do L = 1, 4
            write (20, *) tri(L, 1), tri(L, 2), tri(L, 3)
         enddo
      enddo
      k = (i - 1)*(n(1) + 1)
      do j = 1, n(1) + 1
         x(1) = xmin(1) + dx(1)*(j - 1)
         x(2) = xmin(2) + dx(2)*(i - 1)
         write (10, *) x(1), x(2)
      enddo
   enddo

   do j = 1, n(1) + 1
      x(1) = xmin(1) + dx(1)*(j - 1)
      x(2) = xmin(2) + dx(2)*n(2)
      write (10, *) x(1), x(2)
   enddo

   close (10)
   close (20)

   return
end subroutine Triangles

subroutine Quads(N, xmin, xmax, dx)
   integer, intent(in) :: n(2)
   double precision, intent(in) :: xmin(2), xmax(2), dx(2)
   integer :: i, j

   write (6, *) 'Generating mesh with quadrilateral elements'
   write (6, *) 'Mesh dimensions:'
   write (6, *) '  xmin=', xmin(1)
   write (6, *) '  ymin=', xmin(2)
   write (6, *) '  xmax=', xmax(1)
   write (6, *) '  ymax=', xmax(2)
   write (6, *) '  dx  =', dx(1)
   write (6, *) '  dy  =', dx(2)

   open (unit=10, file='square.cell')
   write (10, *) (n(1) - 1)*(n(2) - 1), 4
   do j = 1, n(2) - 1
      do i = 1, n(1) - 1
         write (10, *) i + (j - 1)*n(1), i + 1 + (j - 1)*n(1), i + 1 + j*n(1), i + j*n(1)
      enddo
   enddo
   close (10)

   open (unit=10, file='square.node')
   write (10, *) n(1)*n(2)
   do j = 1, n(2)
      do i = 1, n(1)
         write (10, *) xmin(1) + (i - 1)*dx(1), xmin(2) + (j - 1)*dx(2)
      enddo
   enddo
   close (10)

   return
end subroutine Quads
