! code written by F. Duncan M. Haldane, haldane@princeton.edu
! v0.3,  24  February 2011
!=================================================================
!   A FORTRAN 95 module for computing the monomial expansion of Jack symmetric polynomials, as well their
!   antisymmetric counterparts (following B.A. Bernevig and N. Regnault, Phys. Rev. Lett. (2009) /  arXiv:0902.4320v2)
!   
!   A fermionic Jack  Q^{alpha}_{lambda+lambda_0}(z) =  V(z)J^{alpha}_{lambda}(z), where J is a symmetric Jack
!   polynomial with Jack parameter alpha, and root partition lambda. V(z) is the Vandermonde determinant and lambda_0
!   is its root partition.      Jacks here are normalized so the coefficient of the root monomial in the monomial expansion
!   is 1.
!----------------------------------------------------------------------------------------------------------------------   
! public subroutines
! 
!  MAKE_SQUEEZED_BASIS(NUM,NORB,MAX_NBASIS,STAT,NBASIS)  creates a searchable basis of states dominated by NUM(1:NORB)
!  which is a bosonic (STAT=1) or fermionic (STAT=-1) occupation configuration of NORB orbital.  NUM(I),
!  I = 1,.., NORB is the  occupation.   
!  MAX_NBASIS is a maximum allowed basis-set size used to limit allocated size  of a temporary workspace for constructing  the basis set.    
!  NBASIS is returned as the actual size found.   (MAX_NBASIS must be at least  NBASIS)
!  A call to  MAKE_SQUEEZED_BASIS with NORB=0 deletes the basis-set  data, and any Jack data.
!
!  subroutine GET_CONFIG(ID,NUM,NORB) returns the occupation number configuration  NUM(1:norb) of basis state ID = 1:NBASIS
!
!  subroutine GET_PARTITION(ID,PARTITION,N) returns the equivalent padded partition PARTITION(1:N) of N = sum(NUM) parts.  
! 
!  integer function CONFIG_ID(NUM,NORB) returns the basis-state number of configuration  NUM(1:NORB),  or 0 if 
!  it is not present in the basis set.
!
!  subroutine GET_BASIS_INFO(NBASIS,NORB,STAT,N,M) retrieves information about the stored basis set: in addition to
!  NBASIS, NORB, and STAT, it returns:
!        N  = SUM(NUM) = sum_{i=1:NORB} NUM(i)   (this is the number of particles (independent  variables))  
!        M = sum_{i=1:NORB}NUM(i)*(i-1)          (equal to |partition| )
!  which are common to all members of the basis set.     NBASIS = 0 is returned in no basis set is present. 
!  (GET_CONFIG with ID = 1 gives the root configuration used in the most recent  call to MAKE_SQUEEZED_BASIS.)
!
!  A call to STATE_TABLE(UNIT) prints a basis-set table to output unit=UNIT.
!--------------------------------------------------------------------------------------------------
! Calculation of JACKS (after the "squeezed basis" derived from the Jack root partion has been created):
!
!   subroutine JACK(ALPHA_JACK,NBASIS_JACK,JACK_P)
!   subroutine RATIONAL_JACK(ALPHA_JACK,NBASIS_JACK,JACK_P,NFOUND_JACK)
!
!   These  produce the monomial expansions
!
!                         J(alpha) = sum_(i=1:NBASIS_JACK) C(i)*m_i,     
!
!   where m_i are N-particle monomials in the  basis where  m_i , i> i, have partions dominated by that of m_1.
!   The "monic normalization" where C(1) = 1 is used.  Fermionic jacks are defined as a Vandermonde determinant times a symmetric
!   Jack with jack parameter alpha.
!
!
!    The partition of m_i (in padded form) is given by PARTITION(N) returned by GET_PARTITION(1,PARTITION,N)
!    The partition of the Jack is that of m_1.     
!
!    ALPHA_JACK is the symmetric-Jack "alpha" parameter.
!
!   in JACK, either:
!      ALPHA_JACK is double precision real, and JACK_P is a dimension(:) pointer to a double precision real array which on return is associated to an
!      array containing C(1:NBASIS_JACK.) 
!   or
!      ALPHA_JACK is double precision complex, and JACK_P is a dimension(:) pointer to a double precision complex array which on return is
!       assiocaietd to an array containing  C(1:NBASIS_JACK).
!      
!
!   in  RATIONAL_JACK
!       ALPHA_JACK is integer(2), ("alpha" = ALPHA_JACK(1)/ALPHA_JACK(2)) and JACK_P is a dimension(:,:) pointer to a "long" integer array, which on
!       return is associqted to an array containing  C(2,NBASIS_JACK),  where "C(i)" = C(1,i)/C(2,i).    "long" is integer (kind=8) (64bit).
!      
!   In this rational case, long integers may have insufficient length  to represent all monomial expansion coefficients.    If integer overflow 
!   is encountered, only the first NFOUND_JACK coefficents are calculated.   This condition is signalled by NFOUND_JACK \= NBASIS_JACK on return.
!
!   In all cases, the "monic-normalization" Jacks for some root partitions are singular when "alpha" is a negative rational.  
!   These "bad"  cases cannot be  calculated by this procedure.  (The Jack must have an alternative normalization where C(1) = 0).   
!   If this condition (a "divide-by-zero") is encounted, JACK_P is returned disassociated, and NBASIS_JACK = -1 is returned.        
!
!   subroutine DELETE_JACK can be called to deallocate the stored array originally pointed to by JACK_P in case this pointer has become disassocaited
!   from it.   Any call to JACK or RATIONAL_JACK will deallocate any previously-allocated Jack expansion data, so a copy must be made of it and its basis
!   set if it is is wished to retain it.
!
! Examples:
!  The Jack parameter  alpha = -2 produces nu = 1/2 bosonic Laughlin states or nu = 1/3 fermionic ones.
!  The Jack parameter  alpha = -3 produces nu = 2/2 bosonic Moore-Read states, or nu = 2/4 fermonic ones.  
!----------------------------------------------------------------------------------------------------- 
module jacks_m
   use iso_fortran_env
  private
  integer, parameter :: dp = kind(1.0d0)    !double precision real or complex.
  integer, parameter :: long=8              !long integers


  public make_squeezed_basis, get_config, get_partition, config_id, state_table, get_basis_info 
  public jack, rational_jack, delete_jack

!-----------------------------------------------------------------------------------------------
! nothing else public 

  interface jack
     module procedure real_jack, complex_jack
  end interface

  integer, allocatable :: row_elements(:,:)
  integer, allocatable :: state_list_temp(:,:),state_list(:,:)
  integer, allocatable :: secondary_index_temp(:,:),secondary_index(:,:)
  integer, allocatable :: sector_sum(:), sectors (:,:) 
  integer, allocatable :: num1(:),ibin1(:,:),num2(:)
  real(kind=dp) , allocatable, target :: djack(:)
  complex(kind=dp), allocatable, target  :: zjack(:)
  integer (kind=16), allocatable, target :: ijack(:,:)


  integer :: len, central_orbital, min_sum, max_sum, nsectors,count
  integer :: allocated_list_space
  integer :: nbasis_stored, norb_stored, n_stored, m_stored
  logical :: fermions
  logical :: debug = .false.
  
contains
  subroutine delete_jack
    implicit none
    if(allocated (djack)) deallocate(djack)
    if(allocated (zjack)) deallocate(zjack)
    if(allocated (ijack)) deallocate(ijack)
    return
  end subroutine delete_jack

  subroutine get_basis_info(nbasis,norb,stat,n,m)
    implicit none
    integer, intent(out) :: nbasis,stat,norb,n,m
!-------------------------------------------------
! check if a stored squeezed basis exists, and if so,
! report nbasis, stat, n = sum_(i=1:norb} num(i), 
! m = sum_(i=1:norb} (i-1)*num(i), 
! if no stired basis exists, nbasis =istat = norb = m = n = 0 is returned
!-----------------------------------------------
    nbasis = 0
    norb = 0
    n = 0
    m = 0

    if(.not.allocated(state_list)) return

    nbasis = nbasis_stored
    norb = norb_stored
    stat = 1
    if(fermions) stat = -1
    n = n_stored
    m = m_stored
    return
  end subroutine get_basis_info

  subroutine make_squeezed_basis(num0,norb,max_nbasis,stat,nbasis)
    implicit none
    integer, intent(in) :: norb,stat
    integer(kind=16),intent(in) :: max_nbasis
    integer, intent(in) ::  num0(norb)
    integer, intent(out) :: nbasis
    !-----------------------------------------------------------------------------
    !  Sets up  a searchable list of the NBASIS bosonic (STAT=1) or fermionic (STAT=-1)
    !  configurations dominated by NUM0((NORB) (plus NUM0(NORB) itself).
    !
    !  After the list is created, the integer function CONFIG_ID(NUM,NORB) 
    !  returns the index (or zero if NUM is not in the list).
    !  If config_id(num1) >= config_id(num2) ,then  num1 cannot be dominated by num2.
    
    !  The subroutine GET_CONFIG(ID,NUM,NORB)  returns the occupation numbers 
    !  NUM(1),..,NUM(NORB) of the configuration with index ID, ID = 1,2,...,NBASIS.
    !
    !  A call to MAKE_SQUEEZED_BASIS with NORB = 0 will release all allocated storage.
    !  The calling argument MAX_NBASIS sets an upper limit for NBASIS.  The program
    !  terminates if this is too small.
    !---------------------------------------------------------------------------------
    integer, allocatable :: sum1(:),order(:)
    integer :: i,j,k,l,s,id,i0,i1,i2,i3,pos_new
    integer :: sublist, sector, sector_length, pos1,pos2,pos



    if(allocated (row_elements)) deallocate(row_elements)
    if(allocated (state_list)) deallocate(state_list)
    if(allocated (sectors)) deallocate(sectors)
    if(allocated (sector_sum)) deallocate(sector_sum)
    if(allocated (secondary_index)) deallocate(secondary_index)
    
    call delete_jack    ! make sure anything allocated in a previous call is deallocated.
    

    allocated_list_space = 0
    nbasis_stored = 0
    norb_stored = 0
    fermions = .false.
    n_stored = 0
    m_stored = 0

    len = 0
    central_orbital = 0

    nbasis = 0
    if(norb==0) return

    norb_stored = norb
    if(norb<0.or.max_nbasis<1) then
       write(6,1) norb,max_nbasis
1      format('MAKE_SQUEEZED_BASIS: invalid calling arguments NORB=',i12,' MAX_NBASIS=',i12)
       stop
    endif
    
    if(stat.eq.-1) then
       fermions = .true.
    else if(stat.ne.1) then
       write(6,2) stat
2      format('MAKE_SQUEEZED_BASIS: called with STAT =',i3,' ONLY +1 or -1 are valid')
       stop
    endif

    n_stored = sum(num0)
    if(n_stored == 0) then
       write(6,3) 
3      format('MAKE_SQUEEZED_BASIS: called with empty configuration NUM0 = 0, (invalid)')
       stop
    endif

    if(minval(num0) < 0) then
       write(6,4) num0
4      format('MAKE_SQUEEZED_BASIS: invalid configuration NUM0= ',/20i3,5(/20i3))
    endif
    
    if(fermions.and. maxval(num0) > 1) then
       write(6,5) num0
5      format('MAKE_SQUEEZED_BASIS: invalid fermion configuration NUM0= ',/20i3,5(/20i3))
    endif

    
    ! determine a "central orbital" so the values of the  code used to order the
    ! list will be  as small as possible.

    i0 = n_stored
    i1  = 0
    do i = 1,norb
       i1 = i1 + num0(i)*(i-1)
    enddo
    m_stored = i1

    central_orbital = (i1+i0)/i0
    
    !each configuration is stored in a binary code as LEN integers
    !  the details of the encoding are known only to subroutine CODE
    ! a call to CODE with input LEN = 0 returns the value of LEN
 
    call code(norb,num0,len,ibin1,.false.)  ! ibin1 not referenced when called with .false. 

    allocate (ibin1(len,1),num1(norb))  ! small arrays used in creating lists


    ! storage for up to MAX_NBASIS configurations is allocated
    ! if this is not enough the program terminates.
    ! when the list is complete, it has NBASIS entries, which are  reordered 
    ! and copied to newly-allocated space of the correct size.  The inital allocation
    ! is then deallocated.
    
    allocate (state_list_temp(len,max_nbasis))
    allocated_list_space = max_nbasis

    ! Call a loop that generates all configurations dominated by 
    ! the root configuration NUM0(1),..,NUM0(NORB).
    ! This is a recursive loop with a depth equal to the number of particles, 
    ! with no limit on the number of orbitals NORB.
    allocate (num2(norb)) ! used recursively in loop/center
    call loop(norb,num0)
    deallocate(num2)

    nbasis = nbasis_stored


    !sort codes first by  decreasing value of ORDERING_SUM(NUM,NORB)
    ! this guarantees they are in (increasing) dominance order,
    ! with NUM0 as the top code, with ID = 1.
   
   !****THIS IS PRIMARY ORDERING IN REVERSE ORDER OF THE "FREE-PARTICLE KINETIC ENERGY" ****
    allocate(sum1(nbasis),order(nbasis))
    do i = 1,nbasis
       call code(norb,num1,len,state_list_temp(1,i),.true.)  
       sum1(i) = -ordering_sum(num1,norb)      
    enddo
    ! call a routine to index the array sum1:
    call list_order(nbasis,sum1,order)
    min_sum = -sum1(order(nbasis));  max_sum = -sum1(order(1))  
    nsectors = 1 + ((max_sum-min_sum)/2) ! values of sum1  are  min_sum, min_sum + 2,  ...  max_sum}
    ! efficient hash code is 1 + (max_sum+sum1)/2 = 1+(max_sum-ordering_sum)/2.

    allocate(state_list(len,nbasis))     
    ! copy the codes in order  to STATE_LIST(LEN,NBASIS)
    forall(i=1:nbasis) state_list(:,i) = state_list_temp(:,order(i))
    deallocate(state_list_temp)      ! release the initial temp  array


    ! break up the list into sectors with the same value of the state-ordering code.
    ! states with this  ordering-code = i occupy positions  sector(1,i),...,sector(2,i) in the 
    ! list of configurations.
    ! If there are no configurations with ordering-code = i, sectors(2,i) = 0


    allocate (sectors(4,nsectors))
    sectors = 0
    pos =  1
    sectors(1,1:pos) = 1
    sectors(2,1:pos) = 0
    do i = 1,nbasis
       pos_new = 1 + ((max_sum + sum1(order(i)))/2)
       if(pos_new > pos) then    
          pos = pos_new
          sectors(1,pos) = i
       endif
       sectors(2,pos) = i
    enddo
    deallocate(sum1,order)


! secondary organization within a sector of states with  common value of K = sum_i (i**2)*num(i),
! by D = sum_{j<i} (i-j)*num(i)*num(j)  (bosons) or D' = sum{j<i} (i-j-1)*num(i)*num(j) (fermions)


    allocate (secondary_index_temp(3,nbasis))
    sublist = 0
    sector_loop:    do i = 1,nsectors
       if (sectors(2,i) == 0 ) cycle
       pos1 = sectors(1,i)
       pos2 = sectors(2,i)
       sector_length = pos2 - pos1 + 1
       allocate (sum1(sector_length),order(sector_length))
  
       count = 0
       do j = pos1,pos2
          count = count + 1
          ibin1(:,1) = state_list(:,j)
          call code(norb,num1,len,ibin1,.true.)
          sum1(count) = -d_sum(norb,num1)  ! calculate the secondary index
       enddo 
       call list_order(sector_length,sum1,order)
       sum1 = -sum1   ! sum(order(i+1)  <= sum(order(i)) (decreasing order of d_sum, highest value first)
       ! reorder within the primary sum sector, in decreasing order of  the secondary index d_sum(norb,num)   
       allocate(state_list_temp(len,sector_length))
       state_list_temp = state_list(:,pos1:pos2)
       forall (j=1:sector_length) state_list(:,pos1+j-1) = state_list_temp(:,order(j))
       deallocate (state_list_temp)

       ! add the new secondary index values to  secondary_index_temp; (there is at least one)
       ! sectors(3,i) is the first subsector in sector i, sectors(4,i) is the last  one 
       sublist = sublist + 1
       sectors(3,i) = sublist ; sectors(4,i) = sublist
       secondary_index_temp(1,sublist) = pos1 
       secondary_index_temp(2,sublist) = pos1 
       secondary_index_temp(3,sublist) = sum1(order(1))
       do j = 2,sector_length
          if(sum1(order(j)) == sum1(order(j-1))) then 
             secondary_index_temp(2,sublist) = pos1 + j -1
             cycle
          endif
          sublist = sublist+1
          sectors(4,i) = sublist
          secondary_index_temp(1,sublist) = pos1 + j - 1
          secondary_index_temp(2,sublist) = pos1 + j - 1
          secondary_index_temp(3,sublist) = sum1(order(j))
       enddo
       deallocate (sum1,order)
    enddo sector_loop
    
    allocate (secondary_index(3,sublist))
    secondary_index = secondary_index_temp(:,1:sublist)
    deallocate (secondary_index_temp)
    

    ! now reorder the codes within each  sub-sector sector so they are ordered by their
    ! binary code.  This makes the list quickly searchable by bisection
    ! a "heapsort" (sort in place) is used.
    do i = 1,sublist
       pos1 = secondary_index(1,i)
       pos2 = secondary_index(2,i)
       call binsort(state_list(1:len,pos1:pos2),pos2-pos1+1)
    enddo
    
    deallocate (ibin1,num1)   ! the process is complete
    
    ! on exit,  this subroutine has allocated  STATE_LIST(LEN,NBASIS) and
    ! SECTORS(2,MAXSUM), where MAXSUM is the root state value of the state-ordering code.
    ! The values of LEN, NBASIS_STORED,MAXSUM,NORB_STORED and  CENTRAL_ORBITAL  are also saved.
    ! a call to MAKE_SQUEEZED_BASIS with NORB = 0 will release all allocated storage
    
    return
  end subroutine make_squeezed_basis

  subroutine state_table(unit,name,option)
    implicit none
    integer, intent(in) :: unit
    
    integer i,j,k,L,i0,i1,i2,i3, id, sector,s,norb,nbasis,option
    character (len=60) :: form
    character*4 name
    logical exist


    norb = norb_stored
    nbasis = nbasis_stored
    
    ! print state table produced by make_squeezed_basis

    allocate( num1(norb),ibin1(len,1))
    call get_config(1,num1,norb)

    write(unit,10)
10  format (40('='))

    write(unit,20) num1
20  format(' squeezed basis with root: ',/25i3,20(/27x,25i3))
    if(.not.fermions) then
       write(unit,30)
30     format(/' SPINLESS BOSONS')
    else
       write(unit,40)
40     format(/' SPINLESS FERMIONS')
    endif
    write(unit,50) nbasis
50 format(/' nbasis =',i8)
    
    write(unit,60)
60 format(/' N     = sum_{i=1,norb} num(i) ',& 
         &//' L     = sum_{i=1:norb} (i-1)*num(i) ',&
         &//' K_sum = sum_{i=1:norb} ((i-1)**2)*num(i) ')
    if(.not.fermions) then
       write(unit,70)
70     format(/' D_sum = sum_{j<i=1:norb} (i-j*)num(i)*num(j)   (BOSONS) ')
    else
       write(unit,75)
75     format(/' D_sum = sum_{j<i=1:norb} (i-j-1)*num(i)*num(j)   (FERMIONS) ')
    endif
    write(unit, 80)
80  format(//6x,'index ',10x,' sector ',10x,'N',5x,'L',2x,'K_sum',2x,'D_sum',2x,' binary code',20x,'occupation configuration')

    i = 3+ ceiling(log10(1.0*huge(state_list(1,1))))
    write(form,90) len,i,norb
90 format('(i16," =",i16,i5," =",i5,5x,4i6,1x,',i2,'i',i2,',5x,',i2,'i3)')
    do i = 1,nsectors
       if(sectors(2,i) == 0) cycle 
       write(6,100)
100    format(40('-'))
       do j = sectors(1,i), sectors(2,i)
          call  get_config(j,num1,norb)
          sector = 1 + ((max_sum-ordering_sum(num1,norb))/2)
          i0 = 0
          i1 = 0
          i2 = 0
          do k = 1,norb
             i0 = i0 + num1(k)
             i1 = i1 + num1(k)*(k-1)
             i2 = i2 + num1(k)*(k-1)**2
          enddo
          i3 = d_sum(norb,num1)
          id = config_id(num1,norb)
          write(unit,fmt=form) j,id,i,sector,i0,i1,i2,i3,(state_list(s,j),s=1,len),num1
       enddo
    enddo

    if (option==1) then 
      write(unit,10)
      deallocate  (num1,ibin1)
      return
    endif

    inquire(file=name,exist=exist)
    if(exist) open(1000, file=name, status='old')
    if(.not.exist)open(1000, file=name, status='new')
    write(1000,*) nbasis

    do i = 1,nsectors
      if(sectors(2,i) == 0) cycle 
      do j = sectors(1,i), sectors(2,i)
         call  get_config(j,num1,norb)
         sector = 1 + ((max_sum-ordering_sum(num1,norb))/2)
         i0 = 0
         i1 = 0
         i2 = 0
         do k = 1,norb
            i0 = i0 + num1(k)
            i1 = i1 + num1(k)*(k-1)
            i2 = i2 + num1(k)*(k-1)**2
         enddo
         i3 = d_sum(norb,num1)
         id = config_id(num1,norb)        
         write(1000,88) num1
88       format(30i1) 
         write(1000,*)  0


      enddo
   enddo
 

    return
  end subroutine state_table


  integer pure function d_sum(norb,num)
    implicit none
    integer , intent(in) :: norb
    integer , intent(in) :: num(norb)
    integer ::i,j

    d_sum = 0
    do i = 2,norb
       do j = 1,i-1
          if(.not.fermions) then
             d_sum = d_sum + num(i)*num(j)*(i-j) !symmetric case
          else
             d_sum = d_sum + num(i)*num(j)*(i-j-1) !antisymmetric case
          endif
       enddo
    enddo
    return
  end function d_sum

  subroutine get_partition(id,partition,n)
    implicit none
    integer, intent(in) :: id,n
    integer , intent(out) :: partition(n)
    !-----------------------------------------------
    !  get the padded partition of n non-negative parts  corresponding to basis state id
    !------------------------------------------------------
    integer , allocatable :: occ_num(:)
    integer count, i

    if(n/=n_stored) then
       write(6,10) n, n_stored
10     format(' get_partition (jack_list_m): parts number mismatch, n=',2i5)
       stop
    endif


    partition = 0
    allocate (occ_num(norb_stored))
    count = 0
    call get_config(id,occ_num,norb_stored)
    do i = norb_stored, 1, -1
       do while (occ_num(i) > 0) 
          count = count + 1
          partition(count) = i-1
          occ_num(i) = occ_num(i) - 1
       enddo
    enddo
    deallocate (occ_num)
    return
  end subroutine get_partition
 
  subroutine get_config(id,num,norb)
    implicit none
    integer, intent (in) :: id,norb
    integer , intent (out) :: num(norb)
    !================================
    ! given the state id (position in basis set), return the
    ! configuration NUM(1),...,NUM(NORB)
    ! (uses the list constructed by MAKE_SQUEEZED_BASIS)
    !=================================
    
    if(norb/=norb_stored) then
       write(6,10) norb,norb_stored
10     format('GET_CONFIG: mismatch between NORB and NORB_STORED',2i6)
       stop
    endif
    
    if(id<1.or.id>nbasis_stored) then
       num(1:norb) = 0
       return
    endif
    call code(norb,num,len,state_list(1,id),.true.)   
    
    return
  end subroutine get_config
  

  


  integer function config_id(num,norb)
    implicit none
    integer, intent(in) :: norb
    integer, intent(in)  :: num(norb)
    !========================================
    ! given a configuration NUM(1),...,NUM(NORB),
    ! returns its position 1,2,..,NBASIS  in the list, 
    ! or 0 if it is absent.
    !=======================================
    integer, allocatable :: ib(:,:),ib1(:,:)
    integer :: sum, pos,imin,imax,imid,id,test,i0,i1,i,j,i3
    integer :: subsector, val(1)

    config_id = 0
    
    if(norb/=norb_stored) then
       write(6,10) norb,norb_stored
10     format('CONFIG_ID: mismatch between NORB, NORB_STORED:',2i6)
       stop
    endif
    
    val = minval(num)
    if(val(1) < 0 ) then
       write(6,1) num
1      format('config_id: invalid occupation numbers :',20(25i3))
       stop
    endif

    if(fermions) then
       val  = maxval(num)
       if(val(1) > 1) then
          write(6,1) num
          stop
       endif
    endif
    
    i0 = 0
    i1 = 0
    do i = 1,norb
       i0 = i0 + num(i)
       i1 = i1 + num(i)*(i-1)
    enddo
    if(i0/=n_stored.or.i1.ne.m_stored) return
    
    sum = ordering_sum(num,norb)
    ! "hash code" for sector, based on sum    pos =  1+ (max_sum-sum)/2
    pos = 1 + ((max_sum-sum)/2)
    if(sectors(4,pos) == 0 ) return ! there no states in the basis  with this ordering code
    imin = sectors(3,pos)  ! d_sum values are in decreasing order within a sector
    imax = sectors(4,pos)
       
    ! now identity subsector in range sectors(3,pos)... sectors(4,pos)
    i3 = d_sum(norb,num)
    subsector = 0
    do  while(imin <= imax)   ! binary search for subsector
       imid = imin + (imax-imin)/2
       if(secondary_index(3,imid)==i3) then
          subsector = imid
          exit
       else if (i3 <  secondary_index(3,imid)) then
           imin = imid + 1
        else
           imax = imid - 1
        endif
      enddo
    if(subsector == 0) return

    ! binary search  of ordered binary codes in range subsector(1,i):subsector(2:i)
    imin = secondary_index(1,subsector)
    imax = secondary_index(2,subsector)
    allocate(ib(len,1),ib1(len,1))

    call code(norb,num,len,ib,.false.)
    ib1(1:len,1) = state_list(1:len,imin)
    
    call compare(ib,ib1,len,test)
    if(test==-1) then 
       id = 0
       goto 1000
    else if(test==0)  then
       id = imin
       goto 1000
    endif
    
    ib1(1:len,1) = state_list(1:len,imax)
    
    call compare(ib,ib1,len,test)
    if(test==1) then 
       id = 0
       goto 1000
    else if(test==0)  then
       id = imax
       goto 1000
    endif
    
100 continue
    if(imax-imin<=1) then
       id = 0
       goto 1000
    endif
    imid = (imax+imin)/2
    
    ib1(1:len,1) = state_list(1:len,imid)
    
    call compare(ib,ib1,len,test)
    if(test>0) then
       imin = imid
    else if (test<0) then
       imax = imid
    else if(test==0) then
       id = imid
       goto 1000
    endif
    goto 100
    
1000 deallocate(ib,ib1)
    config_id = id
    
    return
  end function config_id
  
  
  subroutine compare(ibin_1,ibin_2,size,test)
    implicit none
    integer size,test
    integer, intent(in)::  ibin_1(size,1),ibin_2(size,1)
    integer i
    
    ! compare two binary codes (size = number of integers used for coding)
    ! test = -1,0, 1
    do i = 1,len
       if(ibin_1(i,1)>ibin_2(i,1)) then
          test = 1
          return
       else if (ibin_1(i,1)<ibin_2(i,1))then
          test = -1 
          return
       endif
    enddo
    test = 0
    return
  end subroutine compare
  
  
  integer function ordering_sum(num,norb)
    implicit none
    integer norb
    integer num(norb)
    
    integer i
    ordering_sum=1
    do i = 1,norb
       ordering_sum = ordering_sum + num(i)*(i-central_orbital)**2
    enddo
    
    return
  end function ordering_sum
  
  
  
  
  subroutine code(norb,num,len,ibin,decode)
    implicit none
    integer norb,len
    integer ibin(len,1),num(norb)
    
    integer pos,i,j,k,n,bits,id
    logical decode
    
    integer,  save:: bitsize
    
    ! codes a boson occupation into n + norb - 1 bits,
    ! with up to bitsize bits in each of integer ibin(1),..,ibin(len)
    
    ! codes a fermion occupation into norb bits.
    
    
    if(decode) then
       if(fermions) then
          j = 0
          do i = 1,len
             do k = 0,bitsize-1
                j = j + 1
                if(j.gt.norb) exit
                if(btest(ibin(i,1),k)) then
                   num(j) = 1
                else
                   num(j) = 0
                endif
             enddo
          enddo
       else       ! bosons
          j = 1
          num = 0
          do i = 1,len
             do k = 0,bitsize-1
                if(btest(ibin(i,1),k)) then
                   num(j) = num(j) + 1
                else
                   j = j+1
                   if(j>norb) return
                endif
             enddo
          enddo
       endif
    endif
    
    
    ! initialize using  a call with len = 0, decode = .false.
    if(len==0) then
       n = 0
       do i = 1,norb
          n = n + num(i)
       enddo
       bitsize = bit_size(ibin(1,1))
       if(fermions) then
          bits = norb
       else
          bits = norb + n - 1
       endif
       len = 1 + ((bits-1)/bitsize )
       return
    endif
    
    if(fermions) then
       i = 0
       do id = 1,len
          ibin(id,1) = 0
          do pos = 0,bitsize-1
             i = i + 1
             if(i.gt.norb) exit
             if(num(i)==1) ibin(id,1) = ibset(ibin(id,1),pos)
          enddo
       enddo
    else
       ! bosons
       pos = 0
       id = 1
       ibin(id,1) = 0
       do i = 1,norb
          do j = 1,num(i)
             if(pos>=bitsize) then
                id  = id  + 1
                if(id>len) then
                   write(6,200) len
200                format('CODE called with too-small len=',i5)
                   stop
                endif
                ibin(id,1) = 0
                pos = pos - bitsize
             endif
             ibin(id,1) = ibset(ibin(id,1),pos)
             pos = pos + 1
          enddo
          pos = pos + 1
       enddo
       do i = id+1,len
          ibin(i,1) = 0
       enddo
    endif
    
    return
  end subroutine code
  
  
  subroutine binsort(list,nlist) 
    implicit none
    integer nlist
    integer , dimension(len,nlist) :: list
    !---------------------------------------------------
    !  heap sort to reorder a list of binary code  arrays, representing
    !  configurations, to allow searching in numerical order of codes.
    !
    ! code adapted from hpsort of numerical recipes (W. H. Press et al.)
    ! 2nd. ed. p329
    !
    !  sorted list is searchable.
    !----------------------------------------------------
    
    integer i,j,l,r,itest,s
    
    
    if(nlist.lt.2) return

    l = (nlist/2) + 1
    r = nlist
    
10  continue
    if(l.gt.1) then
       l = l-1
       ibin1(1:len,1) = list(1:len,l)
    else
       ibin1(1:len,1) = list(1:len,r)
       list(1:len,r) = list(1:len,1)
       r = r -1
       if(r.eq.1) then
          list(1:len,1) = ibin1(1:len,1)
          return
       endif
    endif
    i = l
    j = l+l
20  if(j.le.r) then
       if(j.lt.r) then
          do s = 1,len
             if(list(s,j+1)>list(s,j)) then
                j = j + 1
                exit
             else if (list(s,j+1)<list(s,j)) then
                exit
             endif
          enddo
       endif
       
       itest = 0
       do s = 1,len
          if(list(s,j)>ibin1(s,1)) then
             itest = 1
             exit
          else if(list(s,j)<ibin1(s,1)) then
             exit
          endif
       enddo
       if(itest.eq.1) then
          list(1:len,i) = list(1:len,j)
          i = j
          j = j + j
       else
          j = r + 1
       endif
       goto 20
    endif
    list(1:len,i) = ibin1(1:len,1)
    goto 10
  end subroutine binsort
  
  
  subroutine loop(norb,num)
    implicit none
    integer norb
    integer num(norb)
    integer, allocatable :: sum1(:),sum0(:),lambda(:),mu(:)
    integer n,i,j,sum
    integer imin,imax
    integer :: s=0
    
    if(fermions) s=1
    n = 0
    do i = 1,norb
       n = n + num(i)
    enddo

    allocate (sum1(n),sum0(n),lambda(n),mu(n))
    
    n = 0
    do i = norb,1,-1
       do j = 1,num(i)
          n = n + 1
          lambda(n) = i
       enddo
    enddo
    
    sum = 0
    do i = 1,n
       sum = sum + lambda(i)
       sum0(i) = sum
    enddo
    
    imax = sum0(1)
    if(n.eq.1) then
       imin = sum0(1)
    else
       imin = 1
    endif
    do i = imin,imax
       mu(1) = i
       sum1(1) = i
       if(n.eq.1) then
          call center(n,norb,mu)
          cycle
       endif
       call loop1(2)
    enddo
    deallocate (sum1,sum0,lambda,mu)

    return
    
  contains
    
    recursive subroutine loop1(level)
      implicit none
      integer, intent(in) :: level
      integer :: jj,jmin,jmax
      

      jmax = min(mu(level-1)-s,sum0(level) - sum1(level-1))
      if(n.eq.level) then
         jmin = sum0(level) -sum1(level-1)
      else
         jmin = 1
      endif
      
      do jj = jmin,jmax
         mu(level) = jj
         sum1(level) = sum1(level-1) + jj
         if(n.eq.level) then
            call center(n,norb,mu)
            cycle
         endif
         call loop1(level+1)
      enddo
      return
    end subroutine loop1
    
  end subroutine loop
  
  subroutine center(n,norb,mu)
    implicit none
    integer n, mu(n),i,norb
    
    num1(1:norb) = 0
    do i = 1,n
       num1(mu(i)) = num1(mu(i)) + 1
    enddo
    
    nbasis_stored = nbasis_stored  + 1
    
    if (nbasis_stored> allocated_list_space)  then
!       write(6,10) allocated_list_space
!       write(*,*) 
!10     format(' nbasis is larger than storage =',i12)
!       stop
    endif
    
    call code(norb,num1,len,ibin1,.false.)
    
    state_list_temp(:,nbasis_stored) = ibin1(:,1)
    
    if(debug) then
       ! check that code/decode from config to binary and back is working
       call code(norb,num2,len,ibin1,.true.)
       do i = 1,norb
          if(num2(i)/=num1(i)) then
             write(6,100) nbasis_stored
100          format(' CODE: error when nbasis =',i12)
             write(6,101) ibin1(1:len,1)
101          format(20i16)
             write(6,102) num1(1:norb)
             write(6,102) num2(1:norb)
102          format(30i3)
             stop
          endif
       enddo
    endif

    return
  end subroutine center
  
  subroutine list_order(n,list,order)
    implicit none
    integer n
    integer list(n),order(n)
    !-------------------------------------
    ! sort an integer list, list(1),...,list(n) so
    ! list(order(i)) <=list(order(i+1)) for i = 1,...,n-1
    ! based on  indexx  from "Numerical Recipes" by Press et al.
    !------------------------------------------
    integer low,high,id
    integer item,i,j
    
    
    !  do nothing
    if(n<=0) return
    
    
    ! initial ordering  
    forall (j=1:n) order(j) = j
    if(n==1) return
    
    
    low=n/2+1
    high=n
    
10  continue
    if(low>1)then
       low=low-1
       id=order(low)
       item=list(id)
    else
       id=order(high)
       item=list(id)
       order(high)=order(1)
       high=high-1
       if(high==1)then
          order(1)=id
          return
       endif
    endif
    
    i=low
    j=low+low
    do while (j<=high)
       if(j/=high)then
          if(list(order(j))<list(order(j+1))) j=j+1
       endif
       if(item<list(order(j)))then
          order(i)=order(j)
          i=j
          j=j+j
       else
          j=high+1
       endif
    enddo
    order(i)=id
    
    goto 10
    
  end subroutine list_order
  
  subroutine d_matrix_row(id,k_sum,n_found)
  implicit none
  integer, intent(in) :: id
  integer, intent(out) :: n_found,k_sum
  !----------------------------------------------------
  !    D_b =  \sum_{i<j}  ( (z_i +  z_j)/(z_i - z_j) (z_i d_i - z_j d_j),      d_i = d/d/z_i
  !    D_f = D_b -  \sum_{i<j}  (z_i + z_j)^2/(z_i-z_j)^2 
  !
  !    action of D = D_b  on symmetric monomials, or D = D_f on antisymmetric monomials:
  !
  !    D m(num) =  sum_{num1} D(num1,num) m(num1),    id(num1) >= id(num)  (squeezing increases id(num)) 
  ! 
  !    monomials are first-quantized forms of fermion or boson occupation number states 
  !    (with a modified normalization in the boson case):
  !   
  !    |m(num)>  =  (N!)^{1/2}  X(norb)X(norb-1)....X(2)X(1)|vac>
  !    
  !    X(k)  =   (\frac{1}{num(k)!}) (c^+(k))^{num(k)}, N = sum_k num(k).
  !
  !    and "normalized one-particle wavefunctions"  w_k(z) = z^{k-1} k = 1,2,..., norb.
  !
  !    nonzero entries D(i,j), with row index i = id <= j, of the lower-triangular matrix D(i,j), i,j = 1:nbasis
  !    are given by D(i,j(k)) = row_elements(1,k), j(k) = row_elements(2,k), k = 1,n_found.
  !  
  !    row_elements(1:2,1:n_found) is made public by the module.
  !-----------------------------------------------------
  integer ::  root_sum, root_d_sum, norb
  integer :: i,i1,j1,i2,j2,k,nn
  integer :: unsqueezed_id, mxel, exchange
  
  integer, allocatable :: num(:),num0(:),num1(:),list(:,:),order(:),data(:)
  logical, parameter :: reorder = .true.

  n_found = 0

  if(id < 1.or.id > nbasis_stored) return
  
  norb = norb_stored
  allocate (num(norb),num0(norb),num1(norb))
  
  !  pair_unsqueezing cannot find more than norb**3 unsqueezed configurations

  if(allocated(row_elements)) deallocate (row_elements)
  allocate(row_elements(2,norb**3))

! get ordering_sum of root configuration
! this provides a limit to how much unsqueezing is needed.  
  call get_config(1,num0,norb)
  call get_config(id,num,norb)



  root_d_sum = d_sum(norb,num0)
  root_sum = ordering_sum(num0,norb)
  k_sum = ordering_sum(num,norb) - root_sum

  ! diagonal matrix element for this row.
  n_found = 1
  row_elements(1,n_found) = id
  row_elements(2,n_found) = d_sum(norb,num) - root_d_sum

  ! loop over unsqueezes  (j2,i2) -> (j2-k,i2+k) = (j1,i1),   
  !   (i1,j1) with i1 - j1 > 1 squeezes to (i2,j2) with i1 > i2 >= j2 > j1, i1 + j1 = i2 + j2.
   
  
  do i2 = 1,norb-1
     if (num(i2) == 0) cycle
     do j2 = 2,i2
        if(num(j2) == 0) cycle
        if(i2 == j2) then
           if(fermions) cycle
           if(num(i2) == 1) cycle
        endif

        ! remove a particle from orbitals i2 and j2
        num0 = num
        num0(i2) = num0(i2) - 1
        num0(j2) = num0(j2) - 1
        
        do  k = 1,norb 
           i1 = i2 + k
           if(i1 > norb) exit
           j1 = j2 - k
           if(j1 < 1) exit
           
           ! update fermionic exchange factor for this unsqueeze
           if(fermions) then
              if (k == 1) then
                 exchange = 1
              else
                 if(num0(i1-1) == 1) exchange = -exchange
                 if(num0(j1+1) == 1) exchange = -exchange
              endif
           endif
           
           if(fermions .and. ((num(i1)==1).or.(num(j1)==1))) cycle
           
           ! add a particle to orbitals i1 and j1
           num1 = num0
           num1(i1) = num0(i1) + 1
           num1(j1) = num0(j1) + 1
           if(ordering_sum(num1,norb) > root_sum) exit
           unsqueezed_id = config_id(num1,norb)
           if(unsqueezed_id==0) cycle
           ! add to list 
           n_found = n_found + 1
           row_elements(1,n_found) = unsqueezed_id
           if(fermions) then      ! fermionic case, need exchange factor
              mxel = 2*(i2-j2)
              row_elements(2,n_found)  = mxel*exchange
           else                   ! bosonic case, need multiple occupancy factor
              mxel = 2*(i1-j1)
              if(i2 == j2) then
                 row_elements(2,n_found) = mxel*((num(i2)*(num(i2)-1))/2)
              else
                 row_elements(2,n_found) = mxel*num(i2)*num(j2)
              endif
           endif
        enddo
     enddo
  enddo
  deallocate (num,num0,num1)
      
    if(reorder) then     ! order by unsqueezed_id!
       allocate (order(n_found),data(n_found),list(2,n_found))
       list = row_elements(:,1:n_found)
       deallocate (row_elements)
       data = list(1,:)
       call list_order(n_found,data,order)
       allocate (row_elements(2,n_found))
       forall(i=1:n_found) row_elements(:,i) = list(:,order(i))
       deallocate(list,order,data)
    endif
    return
  end subroutine d_matrix_row

  subroutine get_gcd(aa, bb, gcd)
   use iso_fortran_env
    implicit none
    integer (kind=16) :: aa, bb, gcd
    !-----------------------------------------
    ! greatest common divisor of two long integers aa,bb
    !---------------------------------------
    integer (kind=16)::  p, q, r
    p = abs(aa)
    q = abs(bb)
    if(p == 0 .and. q == 0) then
       write(6,10)
10     format(' get_gcd(0,0,gcd) called  (INVALID) ')
       stop
    endif
    if(q == 0) then
       gcd = p
       return
    endif
    do while (p /= 0)
       r = p
       p = mod(q,p)
       q = r
    enddo
    gcd = q
    return
  end subroutine get_gcd
  
  subroutine rational_axpy(a,x,y,fail)
   use iso_fortran_env
    implicit none
    integer,  intent (in) :: a
    integer (kind=16),  intent(in) :: x(2)
    integer (kind=16),  intent(inout) :: y(2)
    logical, intent(out) :: fail 
    !-----------------------------------
    !  rational y <- a*x + y, 
    !   Y(1)/Y(2) = Y(1)/Y(2) + a*(X(1)/X(2)) 
    !   returns with fail = .true. if long integer overflow 
    !   prevents successful calculation   
    !______________________________________
    integer (kind=16) :: a0, ax, ay, x1, x2, gcd

    fail = .false.
    
    if(y(2) == 0) then
       write(6,10) 
10     format('rational_long_axpy: y had denominator = 0')
       stop
    endif
    
    if(x(2) == 0) then
       write(6,20) 
20     format('rational_long_axpy: x had denominator = 0')
       stop
    endif
    
    if(a == 0 .or. x(1) == 0) return
    
    x1 = x(1)
    x2 = x(2)
    
    a0 = int(a,kind=16)
    call get_gcd(a0,x2,gcd)
    a0 = a0 / gcd
    x2 = x2 / gcd
    
    call get_gcd(x1,x2,gcd)
    x1 = x1 / gcd
    x2 = x2 / gcd
    
    call check_mult(a0,x1,fail);  if(fail) return   
    x1 = x1*a0
    
    if(y(1) == 0) then
       y(1) = x1
       y(2) = x2
       if(y(2) < 0)  y = -y
       return
    endif
    
    call get_gcd(x2,y(2),gcd)
    ax = x2 / gcd
    ay = y(2) / gcd
    
    call check_mult(ax,y(2),fail); if(fail) return   
    y(2) = ax * y(2)
    
    call check_mult(ax,y(1),fail); if(fail) return
    y(1) = ax * y(1)
    
    call check_mult(ay,x1,fail);  if(fail) return
    x1 = ay * x1
    
    call check_add(x1,y(1),fail);  if(fail) return
    y(1) = y(1) + x1
    
    call get_gcd(y(1),y(2),gcd)
    y(1) = y(1)/gcd
    y(2) = y(2)/gcd

    if(y(2) < 0)  y = -y
    
    return
  end subroutine rational_axpy    

  subroutine check_mult(aa, bb,fail)
    implicit none
    integer (kind=16), intent (in):: aa, bb
    logical, intent(out) :: fail
    !---------------------------------------------
    ! test whether two long integers can be multiplied
    !----------------------------------------------
    fail = .false.
    if(abs(aa) > huge(aa) / abs(bb)) fail = .true.
    if (fail) write(*,*) "check_mult",aa,bb
    return
  end subroutine check_mult
  
  subroutine check_add(aa, bb,fail)
   use iso_fortran_env
    implicit none
    integer (kind=16), intent(in):: aa,bb
    logical, intent(out) :: fail
    !-------------------------------------------
    ! test whether two long integers can be added
    !------------------------------------------
    integer :: sign = 1
    fail = .false.
    if (aa < 0) sign = - sign
    if (bb < 0) sign = - sign
    if(sign == -1) return
    if(abs(aa) > huge(aa) - abs(bb)) fail = .true.
    if (fail) write(*,*) "check_add",aa,bb

    return
  end subroutine check_add

  subroutine rational_jack(alpha_jack,nbasis_jack,jack_p,nfound_jack)
   use iso_fortran_env
    implicit none
    integer, intent(in) :: alpha_jack(2)
    integer , intent(out):: nbasis_jack, nfound_jack
    integer (kind=16), pointer, dimension(:,:):: jack_p
!-----------------------------------------
!  if (.not.fermions) build symmetric jacks with alpha = alpha_num/alpha_den
!  or if (fermions) build antisymmetric Jacks given by Vandermonde x symmetric Jack(alpha)
!----------------------------------------------------------------
!   bose case, alpha_bose =  alpha_num/alpha_den
!   fermion case alpha_fermi = alpha_bose/(1-alpha_bose) = alpha_num/(alpha_den-alpha_num)
!----------------------------------------------------------------
!   If the calculation fails because of long-integer overflow, fail = .true. is returned.
!   monomial expansion coefficients (ijack(1,i)/ijack*2,i)) were correctly calculated
!   for i = 1:rational_jack_list
!
!  The monomial expansion coefficents are computed in rational form
!  
!  jack = = sum_{i=1,jack_list) (ijack(1,i)/ijack(2,i))*monomial(i)
!
!  if the jack does not have a monic form, divide-by-zero occurs, and
!  public integer rational_jack_list is set to 0
!
!  The rational form is represented in terms of long integers ijack(2,jack_list)
!-------------------------------------------------------------------
    integer :: i,j,k_sum,n_found,count,num,den,d_sum,factor
    integer (kind=16), parameter :: zero = int(0,kind=16), one = int(1,kind=16)
    integer (kind=16) :: gcd, div
    logical :: fail

    nfound_jack = 0

    if(associated(jack_p)) then
       deallocate(jack_p)
       nullify(jack_p)
    endif

    nbasis_jack = nbasis_stored
    if(nbasis_jack == 0)  return

    if(allocated(ijack)) deallocate(ijack)
    allocate (ijack(2,nbasis_jack))
    jack_p => ijack
    ijack(1,:) = zero
    ijack(2,:) = one
    num = alpha_jack(1)
    den = alpha_jack(2)
    if(fermions) den = den - num ! fermionic case: vandermonde x symmetric jack

    ijack(2,1) = one
    ijack(1,1) = one  ! monic normalization
    nfound_jack = 1

    do i = 2,nbasis_jack
       ijack(1,i) = zero
       ijack(2,i) = one
       call  d_matrix_row(i,k_sum,n_found)
       k_sum = k_sum*num
       do  count = 1,n_found
          j = row_elements(1,count)
          if(j==i) then
             d_sum = row_elements(2,count)*den
             cycle
          endif
          factor = den*row_elements(2,count)
          call rational_axpy(-factor,ijack(1,j),ijack(1,i),fail)
          if(fail) then
             ijack(1,i) = zero
             ijack(2,i) = one
             write(*,*) "fail1"
             return
          endif
       enddo
       if(ijack(1,i)==zero) then 
          nfound_jack= i
          cycle
       endif

       div = int(k_sum + d_sum,kind=16)   ! alpha-dependent divisor must not be 0
       if(div==zero) then
          if(debug) write(6,100) 
100       format(' rational_jack*jack_list_m): divide-by-zero')
          deallocate(ijack)
          nullify(jack_p)
          nbasis_jack = -1
          nfound_jack = nbasis_jack
          write(*,*) "fail2"
          return
       endif

       call get_gcd(ijack(1,i),div,gcd)
       ijack(1,i) = ijack(1,i)/gcd
       div = div/gcd
       call check_mult(ijack(2,i),div,fail) ;
       if(fail) then
          ijack(1,i) = zero
          ijack(2,i) = one
          write(*,*) "fail3"
          return
       endif
       ijack(2,i) = ijack(2,i)*div
       if(ijack(2,i) < zero) ijack(:,i) = -ijack(:,i)
       nfound_jack = i
    enddo

    return
    
  end subroutine rational_Jack
  
  subroutine real_Jack(alpha_jack,nbasis_jack,jack_p)
    implicit none
    real(kind=dp), intent(in):: alpha_jack
    integer, intent(out):: nbasis_jack
    real (kind=dp) , pointer, dimension(:) :: jack_p   ! pointer to the monomial expansion coefficients.
!-----------------------------------------------------------------
!  if (.not.fermions) build symmetric jacks with real parameter alpha 
!  or if (fermions)  antisymmetric Jacks given by Vandermonde x symmetric Jack(alpha)
!----------------------------------------------------------------
!   bose case, alpha  =  alpha_jack
!   fermion case alpha_fermi = alpha_jack/(1-alpha_jack) 
!----------------------------------------------------------------
    integer ::  i, j, k_sum, n_found, count, d_sum,temp
    real(kind=dp), parameter :: one = 1_dp, zero = 0_dp
    real(kind=dp) :: alpha, div
    real(kind=dp) :: tiny
    if(associated(jack_p)) then
       !deallocate(jack_p)
       jack_p => null()
    endif
   write(*,*) "here"

    nbasis_jack = nbasis_stored
    if(nbasis_jack==0) return

    tiny = epsilon(real(1,kind=dp))

    nbasis_jack = nbasis_stored
    if(allocated(djack)) deallocate(djack)
    allocate (djack(nbasis_jack))
    jack_p => djack
 
    djack(1) = one  ! monic normalization

    alpha = alpha_jack 
    if(fermions) alpha = alpha/(one-alpha) ! fermionic case: vandermonde x symmetric jack

    do i = 2,nbasis_jack
       djack(i) = zero
       call  d_matrix_row(i,k_sum,n_found)
       do  count = 1,n_found
          j = row_elements(1,count)
          if(j==i) then
             d_sum = row_elements(2,count)
             cycle
          endif
          djack(i) = djack(i) - djack(j)*row_elements(2,count)
       enddo
       if(abs(djack(i)) < tiny) then
          djack(i) = zero
          cycle
       endif
       div = alpha*k_sum + one*d_sum
       if(abs(div)< tiny) then
          if(debug) write(6,100) 
100       format(' real_jack (nbasis_jack_m): divide-by-zero')
          deallocate(djack)
          jack_p => null()
          nbasis_jack = -1
          return
       endif
       djack(i) = djack(i)/div
    enddo
    return

  end subroutine real_Jack
  
  subroutine complex_Jack(alpha_jack,nbasis_jack,jack_p)
    implicit none
    complex(kind=dp), intent(in):: alpha_jack
    integer, intent(out):: nbasis_jack
    complex (kind=dp) , pointer, dimension(:) :: jack_p   ! pointer to the monomial expansion coefficients.
!------------------------------------------------------------------
!  if (.not.fermions) build symmetric jacks with complex parameter alpha 
!  or if (fermions) build antisymmetric Jacks given by Vandermonde x symmetric Jack(alpha)
!----------------------------------------------------------------
!   bose case, alpha  =  alpha_jack
!   fermion case alpha_fermi = alpha_jack/(1-alpha_jack) 
!----------------------------------------------------------------
    integer ::  i, j, k_sum, n_found, count, d_sum
    complex(kind=dp), parameter :: one = (1_dp,0_dp), zero = (0_dp,0_dp)
    complex(kind=dp) :: alpha, div
    real(kind=dp) :: tiny

    if(associated(jack_p)) then
       deallocate(jack_p)
       jack_p => null()
    endif

    nbasis_jack = nbasis_stored
     
    if(nbasis_jack==0) return

    tiny = epsilon(real(1,kind=dp))
 
    nbasis_jack = nbasis_stored
    if(allocated(zjack)) deallocate(zjack)
    allocate (zjack(nbasis_jack))
    jack_p => zjack
 
    zjack(1) = one  ! monic normalization

    alpha = alpha_jack 
    if(fermions) alpha = alpha/(one-alpha) ! fermionic case: vandermonde x symmetric jack

    do i = 2,nbasis_jack
       zjack(i) = zero
       call  d_matrix_row(i,k_sum,n_found)
       do  count = 1,n_found
          j = row_elements(1,count)
          if(j==i) then
             d_sum = row_elements(2,count)
             cycle
          endif
          zjack(i) = zjack(i) - zjack(j)*row_elements(2,count)
       enddo
       if(abs(zjack(i)) < tiny) then
          zjack(i) = zero
          cycle
       endif
       div = alpha*k_sum + one*d_sum
       if(abs(div)< tiny) then
          if(debug) write(6,100) 
100       format(' complex_jack (nbasis_jack_m): divide-by-zero')
          deallocate(zjack)
          jack_p => null()
          nbasis_jack = -1
          return
       endif
       zjack(i) = zjack(i)/div
    enddo
    return

  end subroutine complex_jack

end module jacks_m

