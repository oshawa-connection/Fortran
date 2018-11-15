module advcon
    contains
!--------------------------------------------------------------------------------------------------------------------------------   
!-----Calculates Calculates the X Velocity
    function XVelo(W,h) 
            implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer ::  its, itts, xsize, ysize
        real (long), allocatable:: XVelo(:,:)
        real (long) :: h, W(:,:), hsq
    
        xsize = size(W(:,1))
        ysize = size(W(1,:))
        allocate (XVelo(xsize,ysize))           
        hsq = 2*h
        do its = 2,(xsize-1)
            do itts = 2, (ysize-1)
                    
                    XVelo(its, itts) = (W(its,itts+1) - W(its,itts-1))/(hsq)
                    
            end do 
        end do     
        XVelo(1,:) = 0
        XVelo(:,1) = 0
        XVelo(:,ysize) = 0
        XVelo(xsize,:) = 0

    end function XVelo 
!--------------------------------------------------------------------------------------------------------------------------------   
!-----Calculates Calculates the Y Velocity
     function YVelo(W,h) 
            implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer ::  its, itts, xsize, ysize
        real (long), allocatable:: YVelo(:,:)
        real (long):: h, W(:,:), hsq
    
        xsize = size(W(:,1))
        ysize = size(W(1,:))
        allocate (YVelo(xsize,ysize))           
        hsq =h*2
        do its = 2,(xsize-1)
            do itts = 2, (ysize-1)
                    
                    YVelo(its, itts) = (W(its+1,itts) - W(its-1,itts))/(hsq)   
            end do 
        end do 
        yvelo = (-1)*yvelo    
        YVelo(1,:) = 0
        YVelo(:,1) = 0
        YVelo(:,ysize) = 0
        YVelo(xsize,:) = 0

    end function YVelo 
!--------------------------------------------------------------------------------------------------------------------------------   
 
!--------------------------------------------------------------------------------------------------------------------------------               
    function ADVECTIONdtdx(currenttemps, h, vxs) 
        implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer :: its, itts, xsize, ysize
        real (long) ::  vxs(:,:), currenttemps(:,:)
        real (long), allocatable:: ADVECTIONdtdx(:,:) 
        real (long) :: h
    
        xsize = size(currenttemps(:,1))
        ysize = size(currenttemps(1,:))
        allocate (ADVECTIONdtdx(xsize,ysize))    

        do its = 1,(xsize)
            do itts = 1, (ysize)
                    if (vxs(its,itts)>0) then
                        
        ADVECTIONdtdx(its,itts) = vxs(its,itts) * (((currenttemps(its,itts) - currenttemps(its-1, itts))/ h))
                    
                    end if 

                    if (vxs(its,itts)<0) then
                    
        ADVECTIONdtdx(its,itts) = vxs(its,itts) * (((currenttemps(its+1,(itts)) - currenttemps(its, itts))/ h))
                    
                    end if 

           
            end do 
        end do 

        ADVECTIONdtdx(1,:) = 0
        ADVECTIONdtdx(:,1) = 0
        ADVECTIONdtdx(:,ysize) = 0
        ADVECTIONdtdx(xsize,:) = 0 
        
    end function ADVECTIONdtdx
!--------------------------------------------------------------------------------------------------------------------------------       
    function ADVECTIONdtdy(currenttemps, h, vys) 
        implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer :: its, itts, xsize, ysize
        real (long) :: vys(:,:)
        real (long), allocatable:: ADVECTIONdtdy(:,:) 
        real (long):: h , currenttemps(:,:)
        xsize = size(currenttemps(:,1)) 
        ysize = size(currenttemps(1,:)) 

        allocate (ADVECTIONdtdy(xsize,ysize))            
        
        do its = 1,(xsize)
            do itts = 1, (ysize)
                    if ((vys(its,itts))>0) then
    ADVECTIONdtdy(its,itts) = (vys(its,itts) * (((currenttemps(its,itts) - currenttemps(its, itts-1))/ h)))
                    end if 

                    if (vys(its,itts)<0) then
    ADVECTIONdtdy(its,itts) = (vys(its,itts) * (((currenttemps(its,(itts+1)) - currenttemps(its, itts))/ h)))

                    end if 
   
            end do 
        end do 

        ADVECTIONdtdy(1,:) = 0
        ADVECTIONdtdy(:,1) = 0
        ADVECTIONdtdy(:,ysize) = 0
        ADVECTIONdtdy(xsize,:) = 0 
    end function ADVECTIONdtdy
!--------------------------------------------------------------------------------------------------------------------------------       
    function laplace(original, ht) 
        implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer :: xgrids, ygrids, its, itts, xsize, ysize
        real (long):: original(:,:)
        real (long), allocatable:: laplace(:,:)
        real (long):: ht, htsq
        xsize = size(original(:,1)) 
        ysize = size(original(1,:))
        allocate (laplace(xsize,ysize))   
     
        htsq= ht**2
        do its = 2,(xsize-1)
            do itts = 2, (ysize-1)            
laplace(its, itts) = (((original(its-1,itts)+original(its+1,itts)+original(its,itts-1)&
    +original(its,itts+1)-4.0*original(its,itts)))/(htsq))          
            end do 
        end do     
        laplace(1,:) = 0
        laplace(:,1) = 0
        laplace(:,ysize) = 0
        laplace(xsize,:) = 0
      end function laplace

    function diff(original, ht) 
        implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer :: xgrids, ygrids, its, itts, xsize, ysize
        real (long):: original(:,:)
        real (long), allocatable:: diff(:,:)
        real (long):: ht, htsq
        xsize = size(original(:,1)) 
        ysize = size(original(1,:))
        allocate (diff(xsize,ysize))   
     
        htsq= ht**2
        do its = 2,(xsize-1)
            do itts = 2, (ysize-1)            
diff(its, itts) = (((original(its-1,itts)+original(its+1,itts)+original(its,itts-1)&
    +original(its,itts+1)-4.0*original(its,itts)))/(htsq))          
            end do 
        end do     
        diff(1,:) = 0
        diff(:,1) = 0
        diff(:,ysize) = 0
        diff(xsize,:) = 0
      end function diff

!--------------------------------------------------------------------------------------------------------------------------------   
function CONVECTIONdtdx(currenttemps, h) 
        implicit none
        integer, parameter :: long=selected_real_kind(15,307)
        integer :: its, itts, xsize, ysize
        real (long)::  xmax, ymax
        real (long), allocatable:: CONVECTIONdtdx(:,:), xgridsystem(:,:), ygridsystem(:,:), currenttemps(:,:)
        real (long) :: h, hsq
    
        xsize = size(currenttemps(:,1))
        ysize = size(currenttemps(1,:))
        allocate (CONVECTIONdtdx(xsize,ysize))            
        hsq = 2*h
        do its = 2,(xsize-1)
            do itts = 1, (ysize)
                   
                CONVECTIONdtdx(its,itts) = ((((currenttemps(its+1,itts) - currenttemps(its-1, itts))/ hsq)))
            end do 
        end do  

        convectiondtdx(1,:) = 0
        convectiondtdx(:,1) = 0
        convectiondtdx(:,ysize) = 0
        convectiondtdx(xsize,:) = 0
    end function CONVECTIONdtdx


end module advcon
    