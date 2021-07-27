program writeinput
implicit none 

! Generate input files for QCHEM for a sets of configurations

real,parameter :: au2ang=0.529177
real,parameter :: ang2rad=0.01745329252
real           :: x1,y1,z1
real           :: x2,y2,z2
real           :: x3,y3,z3
real           :: x_min,y_min,z_min
real           :: x_max,y_max,z_max
real           :: dx,dy,dz
real           :: temp_x,temp_y,temp_z
real           :: pi
real           :: R_IC,theta
integer        :: nptx,npty,nptz
integer        :: i,j,k
integer        :: icount
character(20)  :: filename

!=========================================================
! set some parameters (in Angstrom and degree)
! x: represents the variable R_jacobi 
! y: represents the variable theta_jacobi
! z: represents the variable r_CN 

pi = acos(-1.)

nptx  = 35
npty  = 18 
nptz  = 7 

write(6,*) 'nptx =',nptx
write(6,*) 'npty =',npty
write(6,*) 'nptz =',nptz

x_min = 2.5
x_max = 6.
dx    = (x_max-x_min)/nptx 

write(6,*) 'x_min =',x_min
write(6,*) 'x_max =',x_max
write(6,*) 'dx =',dx

y_min = 0.0
y_max = 180.
dy    = (y_max-y_min)/npty

write(6,*) 'y_min =',y_min
write(6,*) 'y_max =',y_max
write(6,*) 'dy =',dy

z_min = 1.0
z_max = 1.3
dz    = (z_max-z_min)/nptz

write(6,*) 'z_min =',z_min
write(6,*) 'z_max =',z_max
write(6,*) 'dz =',dz

!=========================================================
! define unscanned positions

! define position of I at the coordinates origin
x1 = 0.0
y1 = 0.0
z1 = 0.0

! define position of C along x axis
y2 = 0.0 
z2 = 0.0

! define position of N in the xy plane
z3 = 0.0

!=========================================================
! start loop over coordinates

open(100,file='movie.xyz',action='write')
open(200,file='position.dat',action='write')

icount = 0
do j=1,npty
  do k=1,nptz
    do i=1,nptx

      write(6,*) icount

      temp_x = x_min + (i-1)*dx
      temp_y = y_min + (j-1)*dy
      temp_z = z_min + (k-1)*dz

      write(200,*) icount,temp_x,temp_y,temp_z

      temp_y = temp_y * ang2rad

      CALL jacobi2internal(temp_x,temp_y,temp_z,R_IC,theta)

      x2 = R_IC 
      x3 = x2 + temp_z*cos(theta)  
      y3 = y2 + temp_z*sin(theta)  


      filename = 'qchem_'//char(48+icount/1000) &
                         //char(48+mod(icount,1000)/100) &
                         //char(48+mod(mod(icount,1000),100)/10) &
                         //char(48+mod(mod(mod(icount,1000),100),10)) &
                         //'.coor' 

      open(10,file=filename,action='write')

      write(10,*) '0 1'   !charge and spin  
      write(10,*) 'I',x1,y1,z1
      write(10,*) 'C',x2,y2,z2
      write(10,*) 'N',x3,y3,z3
  
      close(10)

      write(100,*) '3'   !number of particles
      write(100,*) 
      write(100,*) 'I',x1,y1,z1
      write(100,*) 'C',x2,y2,z2
      write(100,*) 'N',x3,y3,z3

      icount = icount + 1
    enddo
  enddo
enddo

close(100)
end program writeinput

subroutine jacobi2internal(temp_x,temp_y,temp_z,R_IC,theta)
implicit none 

! Transform from Jacobi to internal

real :: temp_x,temp_y,temp_z
real :: mass_C,mass_N,CM
real :: q1,q2,pi
real :: R_IC,theta,costheta

! Define some parameters 

mass_N = 14.0067 ! mass of N atom 
mass_C = 12.0107 ! mass of C atom

CM = mass_N*temp_z/(mass_C+mass_N)  ! center of mass (computed from C atom)  

! Periodic boundary conditions are used for pi<|theta|<2*pi

pi=acos(-1.0)
IF(temp_y.gt.pi) THEN
  temp_y=2*pi-temp_y
ELSEIF(temp_y.lt.-pi) THEN
  temp_y=-2*pi-temp_y
END IF

! Transform Jacodi coords to potential coords

q1 = temp_x**2

q2 = 2.0 * CM * temp_x * COS(temp_y)

R_IC = SQRT(q1 - q2 + CM**2)

costheta = (temp_x * COS(temp_y) - CM)/R_IC

if(costheta.gt.1) costheta=1

theta = ACOS(costheta)

return
end subroutine jacobi2internal 
