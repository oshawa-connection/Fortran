module solver 
  contains
  !-----------------------------------------------------------------------------------------------------------------------------------------------
  !-----Possion solver part
  recursive function Vcycle_2DPoisson(u_f,rhs,h, C) result (resV)
    implicit none
    integer, parameter :: long=selected_real_kind(15,307)
    real  resV
    real (long),intent(inout):: u_f(:,:)  ! arguments
    real (long),intent(in)   :: rhs(:,:),h
    integer         :: nx,ny,nxc,nyc, i,j  ! local variables
    real (long),allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
    real (long)           :: alpha=0.7, res_rms, C

    nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size

    if (min(nx,ny)>5) then  ! not the coarsest level so will continue

       allocate(res_f(nx,ny),corr_f(nx,ny), &
            corr_c(nxc,nyc),res_c(nxc,nyc))

       !---------- take 2 iterations on the fine grid--------------
       
      res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c) 
      res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)
      !Also updates the u_f 

       !---------- restrict the residue to the coarse grid --------
      !residue_f is residue on the fine grid. 

       call residue_2DPoisson(u_f,rhs,h,res_f,C) !Gets the residue from the fine grid.
       call restrict(res_f,res_c) !Takes it to coarser grid.

       !---------- solve for the coarse grid correction ----------- 
       corr_c = 0.  
       res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2,C) ! *RECURSIVE CALL* REcalls itself until it reaches the coarsest grid (5)

       !Now it has calculated and updated corr_c, res_c

       !---- prolongate (interpolate) the correction calculated on the coarse grid to the fine grid 
       call prolongate(corr_c,corr_f) !Corr_f is the output of this subroutine

       !---------- correct the fine-grid solution -----------------
       u_f = u_f - corr_f   
       !U_ff is the correction. Not sure why it is negative. 
       !---------- two more smoothing iterations on the fine grid---
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)

       deallocate(res_f,corr_f,res_c,corr_c)

    else  

       !----- coarsest level (ny=5): iterate to get 'exact' solution

       do i = 1,100
          res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)
        
       end do

    end if

    resV = res_rms   ! returns the rms. residue

  end function Vcycle_2DPoisson

!-----------------------------------------------------------------------------------------------------------------------------------------------
!-----Iterations function
  function iteration_2DPoisson(u_ff,rhsf,hf,alphaf,Cf)
      implicit none
      integer, parameter :: long=selected_real_kind(15,307)
      integer :: i,j,    imax, jmax
      real (long):: u_ff(:,:), rhsf(:,:), hf, alphaf, coeff, Cf, af, hfsq
      real (long), allocatable :: residue(:,:), residuesq(:,:), corr_cf(:,:)
      real (long):: iteration_2DPoisson
!When calculating the residue, make sure you are using the GS style. 
    imax = size(u_ff(:,1))
    jmax = size(u_ff(1,:))
    allocate(residue(imax,jmax), residuesq(imax,jmax),corr_cf(imax,jmax))
    residue = 0.0
    !coeff =  (alphaf*(hf**2)/((4.0+Cf)*(hf**2)))
    coeff = (alphaf * (hf**2) )/ (4.0+(Cf*(hf**2)))
    hfsq = hf**2

      do i = 2, imax-1
        do j = 2,jmax-1
          !Calulates the resudue from the previous runs U_ff
         !call residue_2DPoisson(u_ff,rhsf,h,res_f) 
        
!residue(i,j)= ((u_ff(i,j+1) + u_ff(i,j-1) + u_ff(i+1,j) + u_ff(i-1,j) - 4*u_ff(i,j))/(hf**2)) - rhsf(i,j)
             
residue(i,j)= ((u_ff(i,j+1) + u_ff(i,j-1) + u_ff(i+1,j) + u_ff(i-1,j) - (4.0+(Cf*hfsq))*u_ff(i,j))/(hfsq)) - rhsf(i,j)


!residue(i,j)= ((u_ff(i,j+1) + u_ff(i,j-1) + u_ff(i+1,j) + u_ff(i-1,j) - (4+(cf*(hf**2)))*u_ff(i,j))/(hf**2)) - rhsf(i,j)
        
!u_ff(i,j) = u_ff(i,j) + (alphaf*residue(i,j) * ((hf**2)/(4+(Cf*(hf**2)))))


u_ff(i,j) = u_ff(i,j) + (alphaf*residue(i,j) * ((hf**2)/(4+(Cf*hfsq))))
        end do
      end do

         residuesq = residue * residue
         iteration_2DPoisson = ((sum(residuesq)) / hf**2)**0.5

        
  end function iteration_2DPoisson
!-----------------------------------------------------------------------------------------------------------------------------------------------
!-----Restrict - Calculate the residue on the coarse grid.  ! Input the fine residue, output the coarse residue.
subroutine restrict(res_f,res_c)

      implicit none
      integer, parameter :: long=selected_real_kind(15,307)
      integer :: i, it
      real (long):: res_f(:,:), res_c(:,:)
      
      integer :: nxc, nyc , nx, ny

    nx = size(res_f(:,1)) !Must be like 2^integer + 1 
    ny = size(res_f(1,:))
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2 

    !Coarses the grid. Takes every other value of the grid.
    res_c = 0
    do i  =2, nxc-1
        do it = 2,nyc-1

            res_c(i,it) = res_f((i*2)-1,(it*2)-1)    

        end do
    end do
      
  end subroutine restrict

!-----------------------------------------------------------------------------------------------------------------------------------------------
!-----Prologation - Apply the correction calculated on the coarse grid to the fine grid.

subroutine prolongate(coarse, fine)
    implicit none
    integer, parameter :: long=selected_real_kind(15,307)
    integer :: nx, nym, nxc,nyc, i,it, ny, j, ad 
    real (long):: fine(:,:), coarse(:,:)
    

     nx=size(fine,1); ny=size(fine,2) 

     nxc=size(coarse,1); nyc = size(coarse,2)

  do  i = 1,nxc
    do it = 1,nyc
            fine(i*2-1,it*2-1) = coarse(i ,it)

        end do 
  end do

    do i = 1, nxc
        do it =1, nyc-1

                  fine(i*2-1, it*2) = (coarse(i, it) + coarse(i, it+1)) / 2.0

          end do
    end do

         
          fine(1,:) = 0
          fine(:,1) = 0
          fine(nx,:) = 0
          fine(:,ny) =0

   do i = 2, nx-1,2
        do it =2, ny-1

                  fine(i, it) = (fine(i+1, it) + fine(i-1, it)) / 2.0

          end do
          end do
         
end subroutine prolongate

!-----------------------------------------------------------------------------------------------------------------------------------------------
!-----Calculates the residue and outputs it into the res_f array. 
 subroutine residue_2DPoisson(u_fff,rhsfff,hfff,res_fff, cff) 
      implicit none
      integer, parameter :: long=selected_real_kind(15,307)
      integer :: imax, jmax, i, j, ad
      real (long) :: u_fff(:,:), hfff, rhsfff(:,:), res_fff(:,:), cff, hsq
      !Need to change reisude here 
      res_fff = 0

    imax = size(u_fff(:,1))
    jmax = size(u_fff(1,:))
    hsq = (hfff**2)

    do i = 2, imax-1
     do j = 2,jmax-1

 !   res_fff(i,j) = ((u_fff(i,j+1) + u_fff(i,j-1) + u_fff(i+1,j) + u_fff(i-1,j) -&
! (4+(cff*(hfff**2))*u_fff(i,j))/(hfff**2)) - rhsfff(i,j)
res_fff(i,j) = ((u_fff(i,j+1) + u_fff(i,j-1) + u_fff(i+1,j) + u_fff(i-1,j) -&
 (4+(Cff*hsq)   )*u_fff(i,j))/(hfff**2)) - rhsfff(i,j)

!res_fff(i,j) = ((u_fff(i,j+1) + u_fff(i,j-1) + u_fff(i+1,j) + u_fff(i-1,j)-&
 ! (4+(cff*(hfff**2))*u_fff(i,j)))/(hfff**2)) -  rhsfff(i,j)

!residue(i,j) = ((u_ff(i,j+1) + u_ff(i,j-1) + u_ff(i+1,j) + u_ff(i-1,j) - (4*(cf*(hf**2))) * u_ff(i,j))/(hf**2))- rhsf(i,j)

    end do
      end do 
   
  end subroutine residue_2DPoisson
!-------------------------------------------------------------------------------------------------------------------------------- 
!-------------------------------------------------------------------------------------------------------------------------------- 

end module solver


