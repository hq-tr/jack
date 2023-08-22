program jacks_test
  use jacks_m  
  use iso_fortran_env

  implicit none

  integer, parameter :: mytype = selected_int_kind (16)

  integer(kind=16), parameter :: max_storage = 500000000   ! maximum allowed basis set  size
  integer :: rational_alpha(2), jack_list,rational_jack_list
  integer, parameter :: dp = kind(1.0d0), long = 8   !double precision and long (64bit) integers

  real (kind=dp), pointer :: djack_p(:)
  complex  (kind=dp), pointer :: zjack_p(:)
  integer (kind=16), pointer :: ijack_p(:,:)

  real (kind=dp) alpha,x 

 integer, allocatable :: num(:),partition(:),numlist(:,:),con(:)
 integer :: norb, stat,alpha_num, alpha_den,i,ii,nel,nbasis,df,j,k_sum,nf,nbasis0, count
 integer :: nbasis1,stat1,m,norb1   
 character*4 name
 logical exist
 integer (kind=16) coeff
 integer option

  write(6,10)
10 format(' give norb, stat = 1/-1 (bose/fermi)')
  read(5,*) norb, stat

  if(stat /= -1) stat = 1

  allocate (num(norb),con(norb))
  con=0;
  write (6,20) norb
20 format(' give occupation pattern num(1),..,num(',i3,')')
  read(5,*) num
 
  if(norb <= 0) stop

  write(*,*) 'give output name'
  read(*,*) name

  !write(1000,*) 'energy of the state:  0'



  write(6,30) stat, num
30 format(' stat=',i2,3x,/25i3,/25i3)

  call make_squeezed_basis(num,norb,max_storage,stat,nbasis)
  write(6,40) nbasis
40 format('basis size:',i12)

write(*,*) "do you want coefficients? 1 for yes, 0 for no"
read(*,*) option

  call state_table(6,name,option)  ! print a basis set table

  call get_basis_info(nbasis1,norb1,stat1,nel,m)
  write(6,35) nbasis1,norb1,stat1,nel,m
35 format(' get_basis_info: nbasis=',i12,' norb=',i3,' stat=',i3,' nel=',i3,' M=',i5)
  



if (option==0) stop

inquire(file=name,exist=exist)
if(exist) open(1000, file=name, status='old')
if(.not.exist)open(1000, file=name, status='new')
write(1000,*) nbasis

  
  write(6,50)
50 format(' rational Jack parameter alpha =  num/den; give alpha')
  read(5,*) rational_alpha
   call rational_jack(rational_alpha,jack_list,ijack_p,rational_jack_list)
   if(jack_list/=rational_jack_list) then
      write(6,90) nbasis, rational_jack_list
90    format(' rational_jack failed:',i12,' out of ',i12' coefficients were calculated')
   endif

   alpha = real(rational_alpha(1),kind=dp)/rational_alpha(2)
   call jack(alpha,jack_list,djack_p)

   call jack(cmplx(alpha,kind=dp),jack_list,zjack_p)
   write(6,60) nbasis, jack_list, rational_jack_list

60 format('computed jack_list', 3i20)

   allocate (partition(nel))

   
   do i = 1,jack_list
      call get_partition(i, partition,nel)
      x = real(ijack_p(1,i),kind=dp)/ijack_p(2,i)
      write (6,70) i, ijack_p(1,i), ijack_p(2,i),x, djack_p(i),zjack_p(i),partition

70    format(i10,2i20,4d14.6,5x,20i3)

      do ii=1,nel

         con(partition(ii)+1)=1

      enddo
      write(1000,61) con
61    format(30i1) 

      if (option==1) then
         write(1000,*)  djack_p(i)
62       format(1i10)
      endif

      con=0
    enddo

   stop
 end program jacks_test

