program probe_cif_adps
  use TYPES_MODULE, only: MOLECULE_TYPE
  use SYSTEM_MODULE, only: create_system => create_, tonto
  use TEXTFILE_MODULE, only: create_, open_for_, stdout, stdin, std_err
  use MOLECULE_MAIN_MODULE, only: run_
  implicit none
  type(MOLECULE_TYPE), allocatable :: m
  character(1024) :: input,basis
  integer :: a
  call get_command_argument(1,input)
  call get_command_argument(2,basis)
  call create_system(tonto)
  call create_(stdout,'probe.stdout')
  call open_for_(stdout,'write')
  call create_(std_err,'probe.stderr')
  call open_for_(std_err,'write')
  call create_(stdin,trim(input))
  call open_for_(stdin,'read')
  call run_(m,trim(basis))
  do a=1,size(m%crystal%asymmetric_unit_atom)
    print '(A,1X,12ES24.15)', trim(m%crystal%asymmetric_unit_atom(a)%label), &
      m%crystal%asymmetric_unit_atom(a)%ADP_tensor(1,1), &
      m%crystal%asymmetric_unit_atom(a)%ADP_tensor(2,2), &
      m%crystal%asymmetric_unit_atom(a)%ADP_tensor(3,3), &
      m%crystal%asymmetric_unit_atom(a)%ADP_tensor(1,2), &
      m%crystal%asymmetric_unit_atom(a)%ADP_tensor(1,3), &
      m%crystal%asymmetric_unit_atom(a)%ADP_tensor(2,3), &
      m%crystal%asymmetric_unit_atom(a)%pADP_errors(4:9)
  enddo
end program
