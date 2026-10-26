! Preserve atom identities through periodic partition setup and fit output.
program test_refinement_identity
  use TYPES_MODULE, only: MOLECULE_TYPE, ATOM_TYPE
  use SYSTEM_MODULE, only: create_system => create_, tonto
  use TEXTFILE_MODULE, only: create_, open_for_, stdout, stdin, std_err
  use MOLECULE_MAIN_MODULE, only: run_
  use MOLECULE_SCF_MODULE, only: make_HA_info_
  use VEC_ATOM_MODULE, only: set_pADP_errors_to_
  use ATOM_MODULE, only: set_isotropic_ADP_
  use CRYSTAL_MODULE, only: update_fragment_info_, initialize_fit_data_, &
      refresh_refined_fragment_info_, set_fragment_esds_, set_asym_from_ufrag_cov_mx_
  implicit none
  type(MOLECULE_TYPE), allocatable :: m
  type(ATOM_TYPE), allocatable :: frag(:),asym(:)
  type(ATOM_TYPE) :: source_atom
  character(1024) :: input,basis
  real(8), allocatable :: source_cov(:,:),covariance(:,:),x(:)
  integer, allocatable :: mapping(:),pmap(:)
  real(8) :: expected
  integer :: a,i,j,n,u
  call get_command_argument(1,input)
  call get_command_argument(2,basis)
  call create_system(tonto)
  call create_(stdout,'refinement-identity.stdout')
  call open_for_(stdout,'write')
  call create_(std_err,'refinement-identity.stderr')
  call open_for_(std_err,'write')
  call create_(stdin,trim(input))
  call open_for_(stdin,'read')
  call run_(m,trim(basis))

  ! Imported source atoms may not yet have refinement-error storage.
  source_atom=m%atom(1)
  if(allocated(source_atom%pADP_errors)) deallocate(source_atom%pADP_errors)
  source_atom%position_error=[1d-4,2d-4,3d-4]
  source_atom%U_iso_error=4d-4
  expected=sum([(source_atom%ADP_tensor(i,i),i=1,3)])/3d0
  call set_isotropic_ADP_(source_atom)
  if(.not.allocated(source_atom%pADP_errors)) error stop 'missing source-atom isotropic errors'
  if(maxval(abs(source_atom%pADP_errors(1:3)-[1d-4,2d-4,3d-4]))>1d-15) &
      error stop 'source-atom coordinate errors were lost'
  if(maxval(abs(source_atom%pADP_errors(4:6)-4d-4))>1d-15) error stop 'wrong source Uiso error'
  do i=1,3
    if(abs(source_atom%ADP_tensor(i,i)-expected)>1d-15) error stop 'wrong source isotropic ADP'
  end do

  ! Source order O1,O2,O3; separate HAR fragment order O3,O1,O2.
  frag=m%atom([3,1,2])
  call update_fragment_info_(m%crystal,frag)
  mapping=m%crystal%asym_atom_for_frag_atom
  if(any(mapping/=[3,1,2])) error stop 'test did not construct the required atom map'
  call make_HA_info_(m,skip_ANOs=.true.,skip_moments=.true.)
  if(any(m%crystal%asym_atom_for_frag_atom/=mapping)) &
      error stop 'periodic partition setup replaced the HAR fragment'
  do a=1,3
    if(maxval(abs(m%crystal%fragment_atom(a)%position-frag(a)%position))>1d-12) &
        error stop 'periodic partition setup changed an atom representative'
  end do

  ! No fit is needed: use a positive covariance with distinctive entries.
  if(.not.allocated(m%crystal%xray_data%reflections)) &
      allocate(m%crystal%xray_data%reflections(1))
  call initialize_fit_data_(m%crystal)
  n=size(m%crystal%xray_data%X_fit)
  if(n/=27) error stop 'unexpected synthetic parameter count'
  allocate(source_cov(n,n),covariance(n,n),x(n),pmap(n))
  source_cov=0d0
  do i=1,n
    source_cov(i,i)=1d-8*i
  end do
  do i=1,n-1
    source_cov(i,i+1)=0.1d0*sqrt(source_cov(i,i)*source_cov(i+1,i+1))
    source_cov(i+1,i)=source_cov(i,i+1)
  end do
  do u=1,3
    a=mapping(u)
    pmap(9*a-8:9*a)=[(9*u-9+i,i=1,9)]
  end do
  covariance=source_cov(pmap,pmap)
  m%crystal%xray_data%covariance_mx=source_cov
  call set_asym_from_ufrag_cov_mx_(m%crystal)
  if(maxval(abs(m%crystal%xray_data%covariance_mx-covariance))>1d-18) &
      error stop 'covariance was assigned to the wrong asymmetric atoms'
  do j=1,n
  do i=1,n
    expected=covariance(i,j)/sqrt(covariance(i,i)*covariance(j,j))
    if(abs(m%crystal%xray_data%correlation_mx(i,j)-expected)>1d-12) &
        error stop 'correlations do not follow the asymmetric parameter labels'
  end do
  end do

  ! The transform was checked against the independent permutation above.
  ! Now require bitwise preservation of its actual result across refresh.
  covariance=m%crystal%xray_data%covariance_mx
  call set_pADP_errors_to_(m%crystal%asymmetric_unit_atom,covariance,.false.)
  call set_fragment_esds_(m%crystal)
  x=[(real(i,8)*0.01d0,i=1,n)]
  m%crystal%xray_data%X_fit=x
  asym=m%crystal%asymmetric_unit_atom
  frag=m%crystal%fragment_atom
  call refresh_refined_fragment_info_(m%crystal)
  if(any(m%crystal%xray_data%covariance_mx/=covariance)) &
      error stop 'refresh changed covariance identity'
  if(any(m%crystal%xray_data%X_fit/=x)) error stop 'refresh changed fitted-vector identity'
  if(any(m%crystal%asym_atom_for_frag_atom/=mapping)) error stop 'refresh changed atom mapping'
  do a=1,3
    if(m%crystal%asymmetric_unit_atom(a)%tag/=asym(a)%tag) &
        error stop 'refresh relabeled asymmetric atoms'
    if(maxval(abs(m%crystal%asymmetric_unit_atom(a)%position-asym(a)%position))>1d-12) &
        error stop 'refresh changed asymmetric representatives'
    if(maxval(abs(m%crystal%fragment_atom(a)%position-frag(a)%position))>1d-12) &
        error stop 'refresh changed fragment representatives'
    if(maxval(abs(m%crystal%asymmetric_unit_atom(a)%pADP_vector-asym(a)%pADP_vector))>1d-12) &
        error stop 'refresh changed asymmetric parameters'
    if(maxval(abs(m%crystal%asymmetric_unit_atom(a)%pADP_errors-asym(a)%pADP_errors))>1d-12) &
        error stop 'refresh assigned asymmetric uncertainties to the wrong atom'
    if(maxval(abs(m%crystal%fragment_atom(a)%pADP_errors-frag(a)%pADP_errors))>1d-12) &
        error stop 'refresh assigned fragment uncertainties to the wrong atom'
    if(maxval(abs(m%crystal%asymmetric_unit_atom(a)%position_error-asym(a)%position_error))>1d-12) &
        error stop 'refresh changed asymmetric position uncertainties'
  end do
  print '(A)', 'REFINEMENT IDENTITY PASS: periodic source/HAR separation, covariance, correlations and uncertainties'
end program
