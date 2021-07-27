program interpolator
implicit none 

!================================================================
!   Read spectra matrix ouput (from QChem)			 
!    and perform an interpolation using 			
!    'Inverse distance weighting' method		
!================================================================

integer :: i,j,k,joffset
integer :: nlines,ndim_gridx,ndim_gridy,ndim_gridz,ndim_spec
integer :: nptx,npty,nptz,nx,ny,nz
integer :: ndim,npoints,icount,vol_grid
real*8  :: pvalue,weight,temp
real*8  :: x,y,z,dx,dy,dz,dx1,dy1,dz1,dx2,dy2,dz2,dist
real*8  :: Rmin(5),Rmax(5)
real*8,allocatable :: data(:,:),big_grid(:,:)
real*8,allocatable :: energy(:,:),intens(:,:)
real*8  :: tol

!================================================================
! Define some parameters

nlines=33     !number of transition lines in input file

nptx  = 18    !number of points in x dimension of the original grid 
npty  = 35    !number of points in y dimension of the original grid
nptz  = 7     !number of points in z dimension of the original grid

nx = 0  !order of interpolation in x dimension
ny = 0  !order of interpolation in y dimension
nz = 0  !order of interpolation in z dimension

ndim_gridx=nptx+nx*(nptx-1)  !number of points in x dimension of interpolated grid
ndim_gridy=npty+ny*(npty-1)  !number of points in y dimension of interpolated grid
ndim_gridz=nptz+nz*(nptz-1)  !number of points in z dimension of interpolated grid

write(6,*) 'ndim_gridx = ',ndim_gridx
write(6,*) 'ndim_gridy = ',ndim_gridy
write(6,*) 'ndim_gridz = ',ndim_gridz
write(6,*) 

pvalue=16     !power parameter for the interpolation scheme
tol=1E-5      !tolerance for distance to grid point for the interpolation scheme

write(6,*) 'pvalue = ',pvalue
write(6,*) 

!================================================================
! Set  some parameters

vol_grid=ndim_gridx*ndim_gridy*ndim_gridz

write(6,*) 'vol_grid =',vol_grid
write(6,*) 

!#================================================================
! Import data
! column #1: jacobi angles (rad)
! column #2: jacobi R (au)
! column #3: r_CN (au)
! column #4: energy (eV)
! column #5: intensity (arb. units)

open(1,file='../response_mat.dat',action='read')

ndim=0
do i=1,1000000
  read(1,*,END=10)
  ndim=ndim+1
enddo
write(6,*) 'increase reading parameter for response_mat.txt' 
stop
10 continue

write(6,*) 'ndim = ',ndim
write(6,*) 

npoints=ndim/nlines
write(6,*) 'npoints = ',npoints
write(6,*) 

allocate(data(ndim,5))

rewind(1)
do i=1,ndim
  read(1,*) (data(i,j),j=1,5)
enddo

close(1)

!================================================================
! Generate interpolated grid
! x: jacobi theta
! y: jacobi R
! z: r_CN

Rmin=minval(data,1)
Rmax=maxval(data,1)

write(6,*) 'Rmin =',Rmin
write(6,*) 'Rmax =',Rmax

dx=(Rmax(1)-Rmin(1))/(ndim_gridx-1)
dy=(Rmax(2)-Rmin(2))/(ndim_gridy-1)
dz=(Rmax(3)-Rmin(3))/(ndim_gridz-1)

write(6,*) 'x range = ',Rmin(1),Rmax(1),dx
write(6,*) 'y range = ',Rmin(2),Rmax(2),dy
write(6,*) 'z range = ',Rmin(3),Rmax(3),dz

allocate(big_grid(vol_grid,3))

icount = 1
do i=1,ndim_gridx
  x = Rmin(1) + (i-1)*dx
  do j=1,ndim_gridy
    y = Rmin(2) + (j-1)*dy
    do k=1,ndim_gridz
      z = Rmin(3) + (k-1)*dz
      big_grid(icount,1)=x
      big_grid(icount,2)=y
      big_grid(icount,3)=z
      icount = icount + 1
    enddo
  enddo
enddo

!================================================================
! Interpolate energy and intensity in big grid

allocate(energy(vol_grid,nlines))
allocate(intens(vol_grid,nlines))

do i=1,vol_grid

  if(mod(i,100).eq.0) write(6,*) 'step ',i,' of ',vol_grid

  weight=0.

loop: do j=1,npoints
        joffset = (j-1)*nlines

        dx1 = big_grid(i,1) - data(joffset+1,1)
        dy1 = big_grid(i,2) - data(joffset+1,2)
        dz1 = big_grid(i,3) - data(joffset+1,3)

!        dx2 = 0.5*(big_grid(i,1) + data(joffset+1,1))
!        dy2 = 0.5*(big_grid(i,2) + data(joffset+1,2))
!        dz2 = 0.5*(big_grid(i,3) + data(joffset+1,3))

        dx2 = dx
        dy2 = dy
        dz2 = dz

!        dx2 = 1.
!        dy2 = 1.
!        dz2 = 1.

        dx1 = dx1/dx2
        dy1 = dy1/dy2
        dz1 = dz1/dz2

        if(dx2.eq.0) dx1=0.0
        if(dy2.eq.0) dy1=0.0
        if(dz2.eq.0) dz1=0.0

        dist = dx1**2 + dy1**2 + dz1**2
        dist = dsqrt(dist)

        if(dist.lt.tol) then
          do k=1,nlines
            energy(i,k)=data(joffset+k,4)
            intens(i,k)=data(joffset+k,5)
          enddo
          exit loop
        else
          temp = 1./(dist**pvalue)
          weight=weight+temp
          do k=1,nlines
            energy(i,k)=energy(i,k)+temp*data(joffset+k,4)
            intens(i,k)=intens(i,k)+temp*data(joffset+k,5)
          enddo
        endif

        if(j.eq.npoints) then
          do k=1,nlines
            energy(i,k)=energy(i,k)/weight
            intens(i,k)=intens(i,k)/weight
          enddo
        endif
      enddo loop
enddo

!================================================================
! Save interpolate spectra to file

open(2,file='response_mat_interp.txt',action='write')

do i =1,vol_grid
  do j = 1,nlines
    write(2,*) (big_grid(i,k),k=1,3),energy(i,j),intens(i,j)
  enddo
enddo
close(2)

end program
