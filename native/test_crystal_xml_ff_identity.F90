! XML density stays in source AO order; form-factor columns follow the fit.
! The fixture has two H atoms with different primitive exponents, so it also
! checks that sharing quadrature by element cannot depend on atom ordering.
program test_crystal_xml_ff_identity
#include "macros"
  use TYPES_MODULE
  use SYSTEM_MODULE
  use MOLECULE_MAIN_MODULE, only: run_
  use MOLECULE_RHO_MODULE, only: make_C23_Hirshfeld_atom_FFs_
  use MOLECULE_CE_MODULE, only: phar_defragment_
  use CRYSTAL_MODULE, only: update_fragment_info_
  use DIFFRACTION_DATA_SET_MODULE, only: create_data => create_
  use ATOM_MODULE, only: set_position_
  use TEXTFILE_MODULE, only: create_std_err, create_stdin, create_stdout, &
      open_for_, close_, std_err, stdin, stdout
  implicit none
  type(MOLECULE_TYPE), allocatable :: m
  type(ATOM_TYPE), allocatable :: original(:),frag(:)
  complex(8), allocatable :: reference(:,:),actual(:,:)
  real(8) :: source_xyz(3,3),difference
  integer :: i
  integer, parameter :: permutation(3)=[3,1,2]
  integer :: original_asym_map(3)
  character(1024) :: basis
  call get_command_argument(1,basis)
  TONTO_CREATE
  call create_std_err
  call open_for_(std_err,'write')
  call create_stdout
  call open_for_(stdout,'write')
  call set_stdout_file_(tonto,stdout)
  call close_(std_err)
  call set_std_err_file_(tonto,stdout)
  call create_stdin
  call open_for_(stdin,'read')
  call run_(m,trim(basis))
  ! CIF order is H/He/H; imported AO order remains He/H/H.
  if(any(m%crystal%asymmetric_unit_atom%atomic_number/=[1,2,1])) error stop 'CIF atom order lost on import'
  if(any(m%atom%atomic_number/=[2,1,1])) error stop 'internal AO order changed'
  ! Defragmentation is checked after the independent import/FF checks below.
  if(trim(m%crystal%asymmetric_unit_atom(1)%label)/='H1') error stop 'first CIF label changed'
  if(trim(m%crystal%asymmetric_unit_atom(2)%label)/='He1') error stop 'second CIF label changed'
  if(trim(m%crystal%asymmetric_unit_atom(3)%label)/='H2') error stop 'third CIF label changed'
  if(.not.allocated(m%crystal%xray_data)) call create_data(m%crystal%xray_data)
  allocate(m%crystal%xray_data%reflections(3))
  do i=1,3
    m%crystal%xray_data%reflections(i)%h1=0
    m%crystal%xray_data%reflections(i)%h2=0
    m%crystal%xray_data%reflections(i)%h3=0
  enddo
  m%crystal%xray_data%reflections(1)%h1=1
  m%crystal%xray_data%reflections(2)%h2=1
  m%crystal%xray_data%reflections(3)%h1=1
  m%crystal%xray_data%reflections(3)%h2=1
  m%crystal%xray_data%reflections(3)%h3=1
  m%crystal%xray_data%stockholder_model='cluster'
  original=m%crystal%fragment_atom
  original_asym_map=m%crystal%asym_atom_for_frag_atom
  source_xyz=reshape([(m%atom(i)%position,i=1,3)],[3,3])
  call make_C23_Hirshfeld_atom_FFs_(m,reference)
  if(maxval(abs(reference))<1d-6) error stop 'empty form-factor test'
  call make_C23_Hirshfeld_atom_FFs_(m,actual)
  if(maxval(abs(actual-reference))>1d-10) error stop 'nonrepeatable XML form factors'
  frag=original(permutation)
  call update_fragment_info_(m%crystal,frag)
  if(any(m%crystal%asym_atom_for_frag_atom/=original_asym_map(permutation))) error stop 'test atom map was not constructed'
  call make_C23_Hirshfeld_atom_FFs_(m,actual)
  difference=maxval(abs(actual-reference(:,permutation)))
  print *, 'XML FF permutation max error:',difference
  if(difference>1d-10) error stop 'XML FF columns do not follow the fragment order'
  do i=1,3
    call set_position_(frag(i),original(permutation(i))%position + &
        real(i-2,8)*m%crystal%unit_cell%direct_mx(:,2))
  enddo
  call update_fragment_info_(m%crystal,frag)
  call make_C23_Hirshfeld_atom_FFs_(m,actual)
  difference=maxval(abs(actual-reference(:,permutation)))
  print *, 'XML FF translated-representative max error:',difference
  if(difference>1d-10) error stop 'XML FF changed under integer translations'
  do i=1,3
    if(maxval(abs(m%atom(i)%position-source_xyz(:,i)))>1d-14) error stop 'source basis centre was changed'
  enddo
  call update_fragment_info_(m%crystal,original)
  call phar_defragment_(m)
  do i=1,size(m%crystal%fragment_atom)
    if(.not.allocated(m%crystal%fragment_atom(i)%basis)) error stop 'defragmented atom lost imported basis'
  enddo
  call make_C23_Hirshfeld_atom_FFs_(m,actual)
  if(maxval(abs(actual))<1d-6) error stop 'empty defragmented form factors'
  if(any(m%crystal%asym_atom_for_frag_atom/=original_asym_map)) error stop 'fixture fragment order changed'
  difference=maxval(abs(actual-reference))
  print *, 'XML FF defragmentation max error:',difference
  if(difference>1d-10) error stop 'XML FF changed after defragmentation'
  if(any(m%crystal%asymmetric_unit_atom%atomic_number/=[1,2,1])) error stop 'ASU order changed after defragmentation'
  print *, 'XML FF IDENTITY PASS'
end program
