module routines

  implicit none

  contains

    subroutine bubble(vec,n)

      implicit none
   !==============================================================================
      integer,intent(in) :: n
      integer,intent(inout), dimension(1:n) :: vec
   !==============================================================================
      integer :: i,j,swap,low,lowi
   !==============================================================================
      do i=1,n-1
         low=vec(i)
         lowi=i
         do j=i+1,n
            if( vec(j)<low ) then
               low=vec(j)
               lowi=j
            endif
         enddo
         swap=vec(i)
         vec(i)=vec(lowi)
         vec(lowi)=swap
      enddo
      return
   end subroutine bubble

   subroutine gradient(x,x_1,x_2,x_3,y_1,y_2,y_3,jac)
      implicit none
   !==============================================================================
      double precision, intent(inout), dimension(1:2,1:3) :: x
      double precision, intent(in) :: x_1,x_2,x_3,y_1,y_2,y_3,jac
   !==============================================================================
      double precision :: TOL=1d-10
   !==============================================================================

      if( dabs(jac)<=TOL ) then
         x(1,1) = 0.0d0
         x(2,1) = 0.0d0
         x(1,2) = 0.0d0
         x(2,2) = 0.0d0
         x(1,3) = 0.0d0
         x(2,3) = 0.0d0
      else
         x(1,1) =-(y_3 - y_2)/jac
         x(2,1) = (x_3 - x_2)/jac
         x(1,2) =-(y_1 - y_3)/jac
         x(2,2) = (x_1 - x_3)/jac
         x(1,3) =-(y_2 - y_1)/jac
         x(2,3) = (x_2 - x_1)/jac
      endif

      return
   end subroutine gradient



   subroutine phi_function(phi,eta,nu)

      implicit none
   !==============================================================================
      double precision,dimension(1:3) :: phi
      double precision,intent(in) :: eta,nu
   !==============================================================================

      phi(1)=1-eta-nu
      phi(2)=eta
      phi(3)=nu

      return
   end subroutine phi_function

end module routines



module matrix_solvers

   implicit none

   contains

   subroutine sparse_cg(A,C,X,B,nodes,max_tris,pre)

      implicit none
   !==============================================================================
      integer,intent(in) :: nodes, max_tris, pre
      integer,intent(in),dimension(1:nodes,1:max_tris+1) :: C
      double precision,intent(in),dimension(1:nodes,1:max_tris+1) :: A
      double precision,intent(inout),dimension(1:nodes) :: X
      double precision,intent(in),dimension(1:nodes) :: B
   !==============================================================================
      double precision,dimension(1:nodes,1:max_tris+1) :: PA
      double precision,dimension(1:nodes) :: r, w, z, PB
      double precision :: alpha, beta, ans, denom
      logical :: Converged
      integer :: i, j, count
   !------------------------------------------------------------------------------
      count=0

      !Use value of X at previous time-step as the initial guess for the CG method.
      !precondition

      if( pre==1 ) then
         PA=A; PB=B
      elseif( pre==2 ) then
         PA=0; PB=0
         do i=1,nodes
            !find diagonal elememt in A
            do j=1,max_tris+1
               if( C(i,j)==i ) then
                  denom = A(i,j)
                  exit
               endif
            enddo

            PB(i) = B(i)/denom
            PA(i,:) = A(i,:)/denom
         enddo
      endif

      do i=1,nodes
         ans=0.0d0
         do j=1,max_tris+1
            if( c(i,j)/=0 ) then
               ans = ans + PA(i,j)*X(C(i,j))
            endif
         enddo
         r(i) = PB(i) - ans
      enddo

      Converged = sqrt(dot_product(r,r)).lt.1d-10
      if( Converged ) return
      w = -r

      do i=1,nodes
         ans  = 0.0d0
         z(i) = 0.0d0
         do j=1,max_tris+1
            if( c(i,j)/=0 ) then
               ans = ans + PA(i,j)*w(C(i,j))
            endif
         enddo
         z(i) = z(i) + ans
      enddo
      
      alpha = dot_product(r,w) / dot_product(w,z)
      x = x + alpha*w 

      do
         ! Increment the number of iterations
         count = count + 1
         r = r - alpha*z;
         ! Check the convergence
         Converged = sqrt(dot_product(r,r)).lt.1d-10
         if( Converged.or.count.gt.1000 ) exit
         beta  = dot_product(r,z) / dot_product(w,z)
         w     = -r + beta*w
         do i=1,nodes
            ans=0.0d0
            z(i) =0.0d0
            do j=1,max_tris+1
               if( c(i,j)/=0 ) then
                  ans = ans + PA(i,j)*w(C(i,j))
               endif
            enddo
            z(i) = z(i) + ans
         enddo
         alpha = dot_product(r,w) / dot_product(w,z)
         x     = x + alpha*w
      enddo

!      print*,'converged in ',count,'iterations'

      return

   end subroutine sparse_cg

end module matrix_solvers
