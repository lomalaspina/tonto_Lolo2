! A fragment representative must carry the symmetries of its own position.
program test_fragment_stabilizers
  use TYPES_MODULE, only: MOLECULE_TYPE, ATOM_TYPE
  use SYSTEM_MODULE, only: create_system => create_, tonto
  use TEXTFILE_MODULE, only: create_, open_for_, stdout, stdin, std_err
  use MOLECULE_MAIN_MODULE, only: run_
  use CRYSTAL_MODULE, only: update_fragment_info_, transform_position_
  use ATOM_MODULE, only: set_position_
  implicit none
  type(MOLECULE_TYPE), allocatable :: m
  type(ATOM_TYPE), allocatable :: frag(:)
  real(8) :: p(3),q(3),delta(3)
  integer :: n,f,a,expected_count,bad_old_operations,variant
  logical :: expected,reported,found
  character(1024) :: input,basis
  call get_command_argument(1,input)
  call get_command_argument(2,basis)
  call create_system(tonto)
  call create_(stdout,'fragment-stabilizers.stdout')
  call open_for_(stdout,'write')
  call create_(std_err,'fragment-stabilizers.stderr')
  call open_for_(std_err,'write')
  call create_(stdin,trim(input))
  call open_for_(stdin,'read')
  call run_(m,trim(basis))
  frag=m%crystal%fragment_atom
  p=m%crystal%asymmetric_unit_geometry(:,1)
  found=.false.
  do n=1,m%crystal%spacegroup%n_seitz
    q=p
    call transform_position_(m%crystal,q,n,to_unit_cell=.false.)
    if(abs(dot_product(q,p))<0.5d0*sum(p*p)) then
      found=.true.
      exit
    endif
  enddo
  if(.not.found) error stop 'test did not find a differently oriented site'
  call set_position_(frag(1),matmul(m%crystal%unit_cell%direct_mx,q))
  ! Also permute the atom order and translate a whole cell in later variants.
  do variant=1,3
    if(variant==2) frag=frag([2,1])
    if(variant==3) then
      do f=1,size(frag)
        call set_position_(frag(f),frag(f)%position + real(f,8)*m%crystal%unit_cell%direct_mx(:,1))
      enddo
    endif
    call update_fragment_info_(m%crystal,frag)
    bad_old_operations=0
    do f=1,size(frag)
      a=m%crystal%asym_atom_for_frag_atom(f)
      p=m%crystal%fragment_geometry(:,f)
      expected_count=0
      do n=1,m%crystal%spacegroup%n_seitz
        q=p
        call transform_position_(m%crystal,q,n,to_unit_cell=.false.)
        delta=q-p
        delta=delta-anint(delta)
        expected=maxval(abs(delta))<1d-10
        reported=any(m%crystal%frag_atom_stabilizer(f)%element==n)
        if(expected.neqv.reported) error stop 'fragment stabilizer belongs to the wrong representative'
        if(expected) expected_count=expected_count+1
        if(any(m%crystal%asym_atom_stabilizer(a)%element==n).and..not.expected) &
            bad_old_operations=bad_old_operations+1
      enddo
      if(expected_count/=m%crystal%frag_atom_sym_factor(f)) error stop 'wrong site multiplicity'
      if(expected_count/=m%crystal%asym_atom_sym_factor(a)) error stop 'changed orbit multiplicity'
    enddo
    if(bad_old_operations==0) error stop 'test did not distinguish the old copied stabilizers'
  enddo
  print '(A)', 'FRAGMENT STABILIZERS PASS: rotated, reordered and translated representatives'
end program
