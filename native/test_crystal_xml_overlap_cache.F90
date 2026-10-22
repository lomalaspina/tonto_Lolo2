! Integration contract for the production XML overlap-validation cache.
! Inputs are copied to temporary directories by the Python test driver.
program test_crystal_xml_overlap_cache
# include "macros"
  use TYPES_MODULE
  use SYSTEM_MODULE
  use MOLECULE_MAIN_MODULE, only: run_, cleanup_
  use TEXTFILE_MODULE, only: create_std_err, create_stdin, create_stdout, &
       open_for_, close_, destroy_, flush_, text_, std_err, stdin, stdout
  use MOLECULE_RHO_MODULE, only: validate_CRYSTAL23_basis_overlap_
  use MOLECULE_READ_MODULE, only: process_CIF_AND_C23_XML_
  use MOLECULE_SET_MODULE, only: update_shell_info_
  use GAUSSIAN_DATA_MODULE, only: CRYSTAL23_ordering
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  implicit none
  type(MOLECULE_TYPE), allocatable :: m
  real(8), allocatable :: original_overlap(:,:), original_density(:,:,:)
  character(80) :: scenario
  integer :: i

  call get_command_argument(1,scenario)
  TONTO_CREATE
  call create_std_err
  call open_for_(std_err,"write")
  call create_stdout
  call open_for_(stdout,"write")
  call set_stdout_file_(tonto,stdout)
  call close_(std_err)
  call set_std_err_file_(tonto,stdout)
  call create_stdin
  call open_for_(stdin,"read")
  call run_(m)
  call check_counts(1,0,"initial XML import")
  original_overlap=m%crystal23%central_overlap_mx
  original_density=m%crystal23%density_mx

  call validate_CRYSTAL23_basis_overlap_(m)
  call check_counts(1,1,"unchanged validation")
  if (any(original_overlap/=m%crystal23%central_overlap_mx)) &
    error stop "Cache hit modified source overlap"
  if (any(original_density/=m%crystal23%density_mx)) &
    error stop "Cache hit modified source density"
  call text_(stdout,"CACHE MUTATION: "//trim(scenario))
  call flush_(stdout)

  select case(trim(scenario))
  case("unchanged")
    call validate_CRYSTAL23_basis_overlap_(m)
    call check_counts(1,2,"second unchanged validation")
  case("reread")
    call process_CIF_AND_C23_XML_(m)
    call check_counts(1,0,"same XML reread must reset validation")
    call validate_CRYSTAL23_basis_overlap_(m)
    call check_counts(1,1,"reread subsequent cache hit")
  case default
    select case(trim(scenario))
    case("geometry")
      m%atom(1)%position(1)=m%atom(1)%position(1)+1.d-7
    case("element")
      m%atom(1)%atomic_number=m%atom(1)%atomic_number+1
    case("source_overlap")
      m%crystal23%central_overlap_mx(1,1)=m%crystal23%central_overlap_mx(1,1)+1.d-7
    case("order_flag")
      m%use_CRYSTAL23_order=.not.m%use_CRYSTAL23_order
    case("global_order_flag")
      CRYSTAL23_ordering=.not.CRYSTAL23_ordering
    case("spherical_flag")
      m%use_spherical_basis=.not.m%use_spherical_basis
    case("shell_update")
      call update_shell_info_(m)
      if(m%crystal23%xml_overlap_cache_valid) error stop "Shell update retained valid cache"
    case("exponent")
      m%atom(1)%basis%shell(1)%exponent(1)=m%atom(1)%basis%shell(1)%exponent(1)+1.d-7
      do i=1,m%n_basis
        if(m%basis(i)%label==m%atom(1)%basis%label) &
          m%basis(i)%shell(1)%exponent(1)=m%atom(1)%basis%shell(1)%exponent(1)
      enddo
      call update_shell_info_(m)
    case("contraction")
      m%atom(1)%basis%shell(1)%contraction(1)=m%atom(1)%basis%shell(1)%contraction(1)+1.d-8
      do i=1,m%n_basis
        if(m%basis(i)%label==m%atom(1)%basis%label) &
          m%basis(i)%shell(1)%contraction(1)=m%atom(1)%basis%shell(1)%contraction(1)
      enddo
      call update_shell_info_(m)
    case("nan_overlap")
      m%crystal23%central_overlap_mx(1,1)=ieee_value(0.d0,ieee_quiet_nan)
    case("inf_overlap")
      m%crystal23%central_overlap_mx(1,1)=ieee_value(0.d0,ieee_positive_inf)
    case("bad_dimension")
      deallocate(m%crystal23%central_overlap_mx)
      allocate(m%crystal23%central_overlap_mx(m%n_bf+1,m%n_bf))
      m%crystal23%central_overlap_mx=0.d0
    case("nan_position")
      m%atom(1)%position(1)=ieee_value(0.d0,ieee_quiet_nan)
    case("nan_exponent")
      m%atom(1)%basis%shell(1)%exponent(1)=ieee_value(0.d0,ieee_quiet_nan)
    case("inf_contraction")
      m%atom(1)%basis%shell(1)%contraction(1)=ieee_value(0.d0,ieee_positive_inf)
    case("bad_atom_map")
      m%atom_for_shell(1)=0
    case("bad_ao_map")
      m%first_basis_fn_for_shell(1)=2
    case("bad_basis_map")
      m%basis_shell_for_shell(1)=m%n_unique_shells+1
    case("bad_overlap")
      m%crystal23%central_overlap_mx(1,1)=m%crystal23%central_overlap_mx(1,1)+0.01d0
    case default
      error stop "Unknown XML overlap cache scenario"
    end select
    call validate_CRYSTAL23_basis_overlap_(m)
    call check_counts(2,1,"changed import must recompute")
    call validate_CRYSTAL23_basis_overlap_(m)
    call check_counts(2,2,"unchanged changed-state validation")
  end select

  call text_(stdout,"CACHE SCENARIO ACCEPTED: "//trim(scenario))
  call flush_(stdout)
  call cleanup_
  call destroy_(stdout)
  call destroy_(stdin)
  call destroy_(std_err)
  TONTO_DESTROY
contains
  subroutine check_counts(evaluations,hits,where)
    integer, intent(in) :: evaluations,hits
    character(*), intent(in) :: where
    if (.not.m%crystal23%xml_overlap_cache_valid) then
      print *, "Invalid overlap cache at ",where
      error stop "Validation did not cache success"
    endif
    if(m%crystal23%xml_overlap_cache_evaluations/=evaluations.or. &
       m%crystal23%xml_overlap_cache_hits/=hits) then
      print *, "Unexpected counters at ",where,": ", &
        m%crystal23%xml_overlap_cache_evaluations,m%crystal23%xml_overlap_cache_hits, &
        " expected ",evaluations,hits
      error stop "Wrong overlap evaluation/cache-hit counts"
    endif
  end subroutine
end program
