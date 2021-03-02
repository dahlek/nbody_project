module nbody

! define constants:
real, parameter :: Pi = 3.1415927, Grav = 231.0509 ! Gravitational constant G is in terms of M_saturn, R_saturn, and t_saturn


! set parameters:

! set numbers of particles
integer,parameter :: N = 20009, N_moon = 9, N_ring = N-N_moon ! Number of particles, number of moons, number of ring particles
integer :: N_above, N_above_plus1, N_moon_plus1 ! values of particles above and below plane

! define variables relevant to position, time, energy
real*8,dimension(3,N) :: X,V,g ! 2D arrays for position(X), velocity(V), and acceleration(g)
real*8,dimension(1,N) :: phi,radius,dummy_variable ! a radius and angle for each particle, a dummy variable for finding random positions
real*8,dimension(N) :: mass ! Masses of the particles
real*8,parameter :: eps = 0.001, eps2 = eps**2, num_orbits = 500.0, step_orbits = 1000.0 ! eps = softening
real*8,parameter :: r_inner_ring = 1.2361, r_outer_ring = 2.32288 ! From the inner ring bound to where the f-ring is
real*8 :: E_kin, E_pot, total_E, accel, t, dt, t_orb, t_end ! energy, time

! Set numbers,names for output files:
integer, parameter :: outunit1 = 1 ! position`
integer, parameter :: outunit2 = 2 ! velocity
integer, parameter :: outunit3 = 3 ! acceleration
integer, parameter :: outunit4 = 4 ! total energy
integer, parameter :: outunit5 = 5 ! time
integer, parameter :: outunit7 = 7 ! kinetic energy **change; this output is reserved for the screen
integer, parameter :: outunit8 = 8 ! potential energy


!------------------------------------------------------!

! Subroutines:

contains

subroutine acceleration

do i = 1, N         ! loop over all particles
   do k = 1, 3      ! loop over all dimensions
      g(k,i) = 0.0  ! set velocities of all particles over 3 dimensions to 0
   end do 
end do 

do i = 1, N_moon       ! Loop over Saturn and moons
   do j = i+1,N   ! Calculate their contributions on everything else (ignoring ring particles' effect on each other)
      r = (sqrt(( x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2 + (x(3,i)-x(3,j))**2 + eps2 ))**3 ! distance between particles
         do k=1,3   ! for all three dimensions
            p = (x(k,i) - x(k,j))/r ! Just a constant
            g(k,i) = g(k,i) - mass(j)*p*Grav ! Include G gravitational constant here in the future. new accel for particle i in dimension k
            g(k,j) = g(k,j) + mass(i)*p*Grav ! new accel for particle j in dimension k
         end do 
   end do 
end do 

return

end subroutine acceleration

!------------------------------------------------------!

subroutine energy

   E_pot = 0. ! set potential energy = 0 for all particles
   E_kin = 0.5*SUM(mass*SUM(V**2,1)) ! calculate kinetic energy
   Do i=1,N-1
      Do j=i+1,N
         E_pot = E_pot-((mass(i)*mass(j))/(sqrt((SUM((x(:,j)-x(:,i))**2+eps2)))))
      EndDo
   EndDo

total_E = E_kin+E_pot ! for writing to output file

end subroutine energy

end module nbody


!------------------------------------------------------!


program code

use nbody ! call the module


! --- INITIAL CONDITIONS --- !


! DEFINING MASSES:
mass(1) = 1.0 ! "Saturn" = 1 M_saturn
mass(2) = 2.817e-10 ! "Prometheus" (in M_saturn)
mass(3) = 2.41197e-10 ! Pandora
mass(4) = 8.71478e-12 ! Pan
mass(5) = 1.16197e-11 ! Atlas
mass(6) = 1.35563e-13 ! Daphnis
mass(7) = 0.0002368 ! Titan
mass(8) = 3.34e-9 ! Janus
mass(9) = 9.27e-10 ! Epimetheus
mass(10:N) = (5.28e-8)/(N-N_moon) ! Ring particles (total M_sat/(N-prometheusandsaturn))


! DEFINING MOON, PARTICLE XY POSITIONS:
x(:,1) = (/ 0. , 0. , 0. /) ! place "Saturn" (N = 1) at the center
x(:,2) = (/ 0. , 2.31267, 0. /) ! Placing Prometheus at appropriate position
x(:,3) = (/ 2.312593, 0. , 0. /) ! Pandora
x(:,4) = (/ -2.216426, 0., 0. /) ! Pan
x(:,5) = (/ 0., -2.284221, 0. /) ! Atlas
x(:,6) = (/ 1.0, 2.032183, 0. /) ! Daphnis
x(:,7) = (/ 20.2732, 0., 0. /) ! Titan
x(:,8) = (/ -1.0, 2.302524, 0. /) ! Janus
x(:,9) = (/ 1.0, 2.304588, 0. /)! Epimetheus

! set random spherical coordinates within ring
call random_number(radius)
radius = r_inner_ring + radius * (r_outer_ring-r_inner_ring) ! put inner and outer limits on radius
call random_number(phi)
phi = phi*2*Pi

N_moon_plus1 = N_moon + 1

!Main ring positions (convert spherical to cartesian):
x(1,N_moon_plus1:N) = radius(1,N_moon_plus1:N)*cos(phi(1,N_moon_plus1:N))
x(2,N_moon_plus1:N) = radius(1,N_moon_plus1:N)*sin(phi(1,N_moon_plus1:N))



! DEFINE VELOCITIES:
!Velocities:
v(1,N_moon_plus1:N) = -(sqrt(Grav*mass(1)/radius(1,N_moon_plus1:N)))*SIN(phi(1,N_moon_plus1:N)) ! v_x
v(2,N_moon_plus1:N) = (sqrt(Grav*mass(1)/radius(1,N_moon_plus1:N)))*COS(phi(1,N_moon_plus1:N)) ! v_y

! Velocities of the moons:
v(1,4:N_moon) = -( sqrt(Grav*mass(1)/sqrt(x(1,4:N_moon)**2 +x(2,4:N_moon)**2)) ) *SIN(ATAN(x(2,4:N_moon)/x(1,4:N_moon))) ! v_x
v(2,4:N_moon) = ( sqrt(Grav*mass(1)/sqrt(x(1,4:N_moon)**2 +x(2,4:N_moon)**2)) ) *COS(ATAN(x(2,4:N_moon)/x(1,4:N_moon))) ! v_y

! Prometheus:
v(1,2) = -sqrt(Grav*mass(1)/abs(x(2,2)))
v(2,2) = 0.0

! Pandora:
v(1,3) = 0.0
v(2,3) = sqrt(Grav*mass(1)/x(1,3))

! Saturn:
v(:,1) = (/0.,0.,0./) ! set "Saturn" to be stationary

v(3,:) = 0.0 ! initial velocity in the z-dir is 0. Random motion will be induced by random spatial dist.


! DEFINE Z-POSITIONS:
! split the ring particles into two:
N_above = 4+(N/2)
N_above_plus1 = 5+(N/2)

call random_number(x(3,:))
call random_number(dummy_variable)
! put Saturn, moons in the plane:
x(3,1:N_moon) = 0.0

! use dummy variable to randomly distribute particles in z-direction:
do i = N_moon_plus1,N
   if (dummy_variable(1,i) .ge. 0.5) then
      ! put particles above plane
      x(3,i) = x(3,i)*(3.318e-7)
   end if

   if (dummy_variable(1,i) .lt. 0.5) then
      ! put particles below plane
      x(3,i) = x(3,i)*(3.318e-7)*(-1.0)
   end if
end do


! DEFINING TIME:
t = 0.0
t_orb = 2.*Pi*(((1.0**3.)/(Grav*sum(mass)))**0.5) ! from Kepler's law. Orbital period at Saturn's surface (~10 hours)
dt = t_orb/step_orbits
t_end = num_orbits*t_orb
print *, sum(mass)
print *, '1 orbit takes: ',t_orb
print *, 'dt has been defined as: ',dt
print *, 'The simulation will end when t = ',t_end



! Open the files to save arrays to:
open(unit=outunit1,file="/home/users/dahlek/nbody/code/project/output/pos_draft17.txt",action="write",status="replace")
open(unit=outunit2,file="/home/users/dahlek/nbody/code/project/output/vel_draft17.txt",action="write",status="replace")
open(unit=outunit3,file="/home/users/dahlek/nbody/code/project/output/accel_draft17.txt",action="write",status="replace")
open(unit=outunit4,file="/home/users/dahlek/nbody/code/project/output/tot_E_draft17.txt",action="write",status="replace")
open(unit=outunit5,file="/home/users/dahlek/nbody/code/project/output/time_draft17.txt",action="write",status="replace")
open(unit=outunit7,file="/home/users/dahlek/nbody/code/project/output/kin_E_draft17.txt",action="write",status="replace")
open(unit=outunit8,file="/home/users/dahlek/nbody/code/project/output/pot_E_draft17.txt",action="write",status="replace")


do

   ! Save first time step:
   if (t == 0) then
      call energy
      write(outunit1,*), x
      write(outunit2,*), v
      write(outunit3,*), g
      write(outunit4,*), total_E
      write(outunit5,*), t
      write(outunit7,*), E_kin
      write(outunit8,*), E_pot
   end if
   

   t = t + dt 
   call acceleration

   ! set acceleration to 0 for "Saturn" to make it stationary
   do k = 1,3
      g(k,1) = 0.
   end do

   v = v + g*dt
   x = x + v*dt


   ! Save stuff every 10 iterations
   if (mod(int(t/dt),1000) == 0) then 
      call energy
      write(outunit1,*), x
      write(outunit2,*), v
      write(outunit3,*), g
      write(outunit4,*), total_E
      write(outunit5,*), t
      write(outunit7,*), E_kin
      write(outunit8,*), E_pot
   end if
  

   ! Countdowns

   if (real(mod((t/t_end),10.0)) == 0.005) then
      print *, 'The run is 0.5% complete.'
   end if


   if (real(mod((t/t_end),10.0)) == 0.1) then
      print *, 'The run is 10% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.2) then
      print *, 'The run is 20% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.3) then
      print *, 'The run is 30% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.4) then
      print *, 'The run is 40% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.5) then
      print *, 'The run is 50% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.6) then
      print *, 'The run is 60% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.7) then
      print *, 'The run is 70% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.8) then
      print *, 'The run is 80% complete.'
   end if

   if (real(mod((t/t_end),10.0)) == 0.9) then
      print *, 'The run is 90% complete.'
   end if


   if (t .ge. t_end) exit ! If all of the required orbits are completed, exit loop
   
end do

print *, 'Saving output files...'

! close output files
close(outunit1)
close(outunit2)
close(outunit3)
close(outunit4)
close(outunit5)
close(outunit7)
close(outunit8)

print *, ' '
print *, "CALCULATIONS CORRECT!"


end program code
