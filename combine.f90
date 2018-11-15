program Final
  use solver
  use advcon
  use solverT
  implicit none
  integer, parameter :: long=selected_real_kind(15,307)
  real (long):: hp, results, resv, f_error, error, a_adv, a_dif,&
   kappa, rayleigh, totaltime, timecount, dt, dtd, dta, Pr,&
   start, finish ,CStream =0.0, Ctemp, beta, comega
  real (long),allocatable, target :: Temperature(:,:)
  real (long), pointer :: base(:), top(:), lhs(:), rhs(:)
  
  real (long), allocatable :: rhsp(:,:), LHSd(:,:), xtemperaturegrad(:,:),xomegagrad(:,:),&
   vx(:,:), vy(:,:),advectionx(:,:), diffusion(:,:), advectiony(:,:),omega(:,:),streamfunc(:,:),&
   advectionxW(:,:), advectionyW(:,:), deltOmega(:,:), rhsomega(:,:), rhsTEMP(:,:), test(:,:)
    
  integer :: ito, itto, ny, nx, n, a,input, currentiteration, iteration
  CHARACTER(LEN=32) :: Formats

  namelist /inputs/  Pr, nx, ny, rayleigh, totaltime, kappa, a_adv, a_dif, beta
  open(3, file='nmlistt.txt', status = 'old')
  read(3, inputs) 
  close(3)
  
  formats = "(10F10.4, 5X)" 
  write(formats, '(a, i0, a)') '(',nx, 'f10.3,5X)' 
 
!--------------------------------------------------------------------------------------------------------------------------------
!-----allocations
allocate(streamfunc(nx,ny), rhsp(nx,ny), Temperature(nx,ny), LHSd(nx,ny),&
  xtemperaturegrad(nx,ny), vx(nx,ny), vy(nx,ny), omega(nx,ny),&
  advectionxW(nx,ny),advectionyW(nx,ny),xomegagrad(nx,ny),rhsomega(nx,ny), rhsTEMP(nx,ny), test(nx,ny))

    !Assign pointers

    call random_number(temperature)
    top => temperature(:,1)
    base => temperature(:,ny)
    lhs => temperature(1,:) 
    rhs => temperature(nx,:)

    top = 1.0
    base = 0.0
    lhs = temperature(2,:)
    rhs = temperature(nx-1,:) 


    call random_number(omega)
 


    hp=1.0/(ny-1)


!--------------------------------------------------------------------------------------------------------------------------------
!-----START THE LOOP NOW

open(2, file='toplot.txt')
open(4, file='toplotomega.txt')
open(7, file='toplotstream.txt')


timecount = 0.0
iteration = 100
   


call cpu_time(start)
Do while (timecount<totaltime)
   !Using the centred difference here
   xtemperaturegrad = CONVECTIONdtdx(temperature, hp)
  
  

!--------------------------------------------------------------------------------------------------------------------------------
!-----Solve for the streamfunction from omega  
  
  streamfunc = 0.0
  f_error = ((sum(omega*omega)) / nx*ny)**0.5

  error = 1 !DO NOT REMOVE 
  do while (error>0.001)
                 !Arguments(Solvefor,  From, spacing, C) 
   resv =  Vcycle_2DPoisson(streamfunc,omega,hp,CStream)
  
    error = resv/ f_error

  end do

!We now have the values of streamfunc. Now need to solve for velocity then streamfunction....

!--------------------------------------------------------------------------------------------------------------------------------
!-----Calculate the velocity components

  vx = XVelo(streamfunc,hp) 
  vy = YVelo(streamfunc,hp)

  diffusion = (laplace(Temperature, hp))
  deltOmega = (laplace(omega,hp))

  advectionx = ADVECTIONdtdx(Temperature, hp, vx)
  advectiony = ADVECTIONdtdy(Temperature, hp, vy) 

  !advection of omega itself? 
  advectionxW = ADVECTIONdtdx(omega, hp, vx)
  advectionyW = ADVECTIONdtdy(omega, hp, vx)
  
  
!--------------------------------------------------------------------------------------------------------------------------------
!-----Timestep Calculation   

  !Advection timestep  
  dtd = a_dif*(((hp**2)/max(Pr,1.0))) !This is constant but is here for the sake of readability
  dta = a_adv*(min(hp/(maxval(abs(vx))),hp/(maxval(abs(vy)))))

  if (beta<0.5) then
     dt = min(dtd,dta)
  else

    dt = dta
  end if 
!--------------------------------------------------------------------------------------------------------------------------------
!-----Solve for C and right hand sides
 
  ! read in from name file instead!

  if (beta==0) then
      CTemp = 0
      Comega = 0
      temperature = temperature + dt* ((diffusion) - (advectionx+advectiony))
      omega = omega + dt*( Pr *(laplace(omega,hp)) - (advectionxW+advectionyW) - rayleigh*Pr*xtemperaturegrad)

      lhs = temperature(2,:)
      rhs = temperature(nx-1,:)
      top = 1.0
      base = 0 


  else

    CTemp = (1.0)/(beta*dt)
    Comega = (1.0)/(Pr*beta*dt)


    rhsOMEGA = (-1.0)/(pr*beta*dt)*(omega+dt*(pr*(1-beta)*laplace(omega,hp)-(advectionxW+advectionyW)-rayleigh*pr*xtemperaturegrad))
    rhsTEMP =  ((-1.0)/(beta*dt)) *(Temperature +dt*((1.0-beta)*diffusion - (advectionx+advectiony)))


!--------------------------------------------------------------------------------------------------------------------------------
!-----Solve for the new temperature
  

  Temperature = 0
  error = 1

  do while (error>0.001)
                 !Arguments(Solvefor,  From, spacing, C) 
   resv =  Vcycle_2DPoissonT(temperature,rhsTEMP,hp,CTemp)
  !  resv = iteration_2DPoisson(LHS,rhsp,hp, 0.7)
    error = resv/ f_error

  end do


  
!--------------------------------------------------------------------------------------------------------------------------------
!-----Solve for the new omega
  omega = 0
  error = 1

  do while (error>0.001)
                 !Arguments(Solvefor,  From, spacing, C) 
   resv =  Vcycle_2DPoisson(omega,rhsomega,hp,comega)
  !  resv = iteration_2DPoisson(LHS,rhsp,hp, 0.7)
    error = resv/ f_error

  end do



end if 


!--------------------------------------------------------------------------------------------------------------------------------
!-----Write to files every 100 iterations 

  if (iteration == 100) then 
      write(2, formats) temperature
      write(4, formats) omega
      write(7, formats) streamfunc
      !write(4, formats), streamfunc
      !write(7, formats), timecount
      iteration = 0
      print*, "Time",timecount
      print*, "dt", dt
      
  end if 


  
  

 timecount = timecount + dt
 iteration = iteration +1
end do
call cpu_time(finish)
print*, "------------------------------------------"
print*, "finished in", finish-start, "seconds"




end program Final



