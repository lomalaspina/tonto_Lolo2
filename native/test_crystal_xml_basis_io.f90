program test_crystal_xml_basis_io
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_quiet_nan
  use crystal_xml_basis_io
  implicit none
  character(len=4096) :: directory, filename
  character(len=512) :: error
  character(len=256) :: base(10), rows(10)
  integer :: ns, np, checks, unit, ios
  real(real64) :: nonfinite, empty(0)
  integer, allocatable :: atom(:), angular(:), count(:), first(:)
  real(real64), allocatable :: exponent(:), coefficient(:)

  ! Only a uniquely named scratch file in the build/work directory is written.
  ! Optional first argument selects an already existing work directory.
  directory = '.'
  if (command_argument_count() > 0) call get_command_argument(1, directory)
  filename = trim(directory)//'/xml_companion_basis_parser_test.tmp'
  checks = 0
  call require(xml_basis_finite_vector([0.0_real64,-1.0_real64,huge(1.0_real64)]), 'finite vector')
  nonfinite = ieee_value(0.0_real64,ieee_positive_inf)
  call require(.not.xml_basis_finite_vector([1.0_real64,nonfinite]), 'infinite vector member')
  nonfinite = ieee_value(0.0_real64,ieee_quiet_nan)
  call require(.not.xml_basis_finite_vector([nonfinite,1.0_real64]), 'NaN vector member')
  call require(xml_basis_same_finite_vector([1.0_real64,-2.0_real64], &
                                            [1.0_real64,-2.0_real64]), 'identical finite cache inputs')
  call require(.not.xml_basis_same_finite_vector([1.0_real64], [1.0_real64,2.0_real64]), &
               'different cache input lengths')
  call require(.not.xml_basis_same_finite_vector([1.0_real64], [nearest(1.0_real64,1.0_real64)]), &
               'one ULP cache input change invalidates equality')
  call require(xml_basis_same_finite_vector([0.0_real64], [-0.0_real64]), &
               'signed zeros have equal cache input values')
  call require(xml_basis_same_finite_vector(empty,empty), 'equal empty cache inputs')
  call require(.not.xml_basis_same_finite_vector([nonfinite], [1.0_real64]), &
               'NaN in left cache input')
  call require(.not.xml_basis_same_finite_vector([1.0_real64], [nonfinite]), &
               'NaN in right cache input')
  call require(.not.xml_basis_same_finite_vector([nonfinite], [nonfinite]), &
               'NaN on both sides must not give cache equality')
  nonfinite = ieee_value(0.0_real64,ieee_positive_inf)
  call require(.not.xml_basis_same_finite_vector([nonfinite], [nonfinite]), &
               'infinity on both sides must not give cache equality')
  call require(.not.xml_basis_same_finite_vector([1.0_real64], [nonfinite]), &
               'infinity in right cache input')
  call require(.not.xml_basis_same_finite_vector([-nonfinite], [1.0_real64]), &
               'negative infinity in left cache input')
  base = [character(len=256) :: &
    '6 2 # carbon', '0 1 2 4.0 1.0', '2.0 0.25 -0.2', '0.5 0.75 0.8', &
    '0 3 1 0.0 1.0', '0.2 1.0', '1 1', '0 0 1 1.0 1.0', '1.d0 1.0', '99 0']
  call load_case(base, [6,1], [9,1])
  call require(error == '', 'S/SP/D and D-exponent valid input: '//trim(error))
  call xml_basis_dimensions(ns,np)
  call require(ns == 4 .and. np == 6, 'expanded dimensions')
  allocate(atom(ns), angular(ns), count(ns), first(ns), exponent(np), coefficient(np))
  call xml_basis_copy(atom,angular,count,first,exponent,coefficient,error)
  call require(error == '', 'copy valid basis')
  call require(all(atom == [1,1,1,2]), 'atom ordering')
  call require(all(angular == [0,1,2,0]), 'SP splitting and spherical angular labels')
  call require(all(count == [2,2,1,1]) .and. all(first == [1,3,5,6]), 'primitive indexing')
  call require(maxval(abs(exponent-[2.0_real64,0.5_real64,2.0_real64,0.5_real64, &
                                  0.2_real64,1.0_real64])) < 1.e-14_real64, 'source exponents preserved')
  call require(maxval(abs(coefficient-[0.25_real64,0.75_real64,-0.2_real64,0.8_real64, &
                                     1.0_real64,1.0_real64])) < 1.e-14_real64, 'source coefficients preserved')
  call xml_basis_copy(atom(:1),angular,count,first,exponent,coefficient,error)
  call require(error /= '' .and. all(coefficient == 0.0_real64), 'copy shape failure is safe')

  call load_case(base, [1,6], [1,9])
  call require(error == '', 'unique element mapping allows different XML atom order')
  call xml_basis_copy(atom,angular,count,first,exponent,coefficient,error)
  call require(all(atom == [1,2,2,2]) .and. all(angular == [0,0,1,2]), 'XML atom order determines output')
  rows = base
  rows(1) = '106 2'
  call load_case(rows, [6,1], [9,1])
  call require(error == '', 'NAT=106 means carbon')
  rows(1) = '1006 2'
  call load_case(rows, [6,1], [9,1])
  call require(error == '', 'NAT=1006 means carbon')
  call load_case([character(len=256) :: '6 1','0 0 1 0 1','1 1','106 1','0 0 1 0 1', &
    '2 1','99 0','ATOM_MAP 2','2 6','1 106'], [6,6], [1,1])
  call require(error == '', 'explicit same-element/different-basis map')
  call xml_basis_dimensions(ns,np)
  deallocate(atom,angular,count,first,exponent,coefficient)
  allocate(atom(ns),angular(ns),count(ns),first(ns),exponent(np),coefficient(np))
  call xml_basis_copy(atom,angular,count,first,exponent,coefficient,error)
  call require(all(exponent == [2.0_real64,1.0_real64]), 'explicit map resolves original XML atom indices')
  call load_case([character(len=256) :: '6 1','0 0 1 0 1','1 1','106 1','0 0 1 0 1', &
    '2 1','99 0'], [6,6], [1,1])
  call rejected('ambiguous same-element map', 'ATOM_MAP')
  call load_case([character(len=256) :: '6 1','0 0 1 0 1','1 1','99 0', &
    'ATOM_MAP 2','1 6','1 6'], [6,6], [1,1])
  call rejected('duplicate mapping', 'repeats')
  call load_case([character(len=256) :: '6 1','0 0 1 0 1','1 1','99 0', &
    'ATOM_MAP 1','1 6'], [1], [1])
  call rejected('wrong mapped element', 'atomic number')
  call load_case(base, [6,1], [10,1])
  call rejected('wrong spherical AO count', 'AO count')
  call load_case(base, [8], [1])
  call rejected('missing element', 'missing basis')
  call load_case(base, [6], [1,1])
  call rejected('array length mismatch', 'equal')

  rows = base
  rows(1) = '206 2'
  call load_case(rows,[6,1],[9,1])
  call rejected('ECP atom', 'ECP')
  rows = base
  rows(1) = '6 2 7'
  call load_case(rows,[6,1],[9,1])
  call rejected('extra block token', 'exactly NAT')
  rows = base
  rows(2) = '1 1 2 4 1'
  call load_case(rows,[6,1],[9,1])
  call rejected('internal shell', 'ITYB=0')
  rows(2) = '0 6 2 4 1'
  call load_case(rows,[6,1],[9,1])
  call rejected('unsupported high angular momentum', 'LAT')
  rows(2) = '0 1 2 4 2'
  call load_case(rows,[6,1],[9,1])
  call rejected('implicit exponent scale', 'SCAL')
  rows(2) = '0 1 0 4 1'
  call load_case(rows,[6,1],[9,1])
  call rejected('zero primitives', 'NG')
  rows(2) = '0 1 2 NaN 1'
  call load_case(rows,[6,1],[9,1])
  call rejected('nonfinite shell charge', 'finite')
  rows = base
  rows(3) = '2.0 0.25 -0.2 0.5'
  call load_case(rows,[6,1],[9,1])
  call rejected('extra primitive token', 'primitive record')
  rows(3) = '0 0.25 -0.2'
  call load_case(rows,[6,1],[9,1])
  call rejected('zero exponent', 'positive')
  rows(3) = '-2 0.25 -0.2'
  call load_case(rows,[6,1],[9,1])
  call rejected('negative exponent', 'positive')
  rows(3) = 'Infinity 0.25 -0.2'
  call load_case(rows,[6,1],[9,1])
  call rejected('infinite exponent', 'finite')
  rows(3) = '1e999 0.25 -0.2'
  call load_case(rows,[6,1],[9,1])
  call rejected('overflowed exponent', 'finite')
  rows(3) = '2 1/ -0.2'
  call load_case(rows,[6,1],[9,1])
  call rejected('Fortran list-directed slash', 'finite')
  rows(3) = '2 0 -0.2'
  rows(4) = '0.5 0 0.8'
  call load_case(rows,[6,1],[9,1])
  call rejected('zero S contraction', 'zero')
  rows(3) = '2 0.25 0'
  rows(4) = '0.5 0.75 0'
  call load_case(rows,[6,1],[9,1])
  call rejected('zero SP p contraction', 'zero')
  rows(3) = '2 1 -0.2'
  rows(4) = '2 -1 0.8'
  call load_case(rows,[6,1],[9,1])
  call rejected('identical-exponent cancelling S contraction', 'norm')
  rows(3) = '2 0.25 1'
  rows(4) = '2 0.75 -1'
  call load_case(rows,[6,1],[9,1])
  call rejected('identical-exponent cancelling SP p contraction', 'norm')
  rows = base
  rows(9) = '1 1e200'
  call load_case(rows,[6,1],[9,1])
  call rejected('overflowed contraction norm', 'norm')
  rows = base
  rows(9) = '1e300 1'
  call load_case(rows,[6,1],[9,1])
  call require(error == '', 'finite large exponent norm is evaluated without exponent product overflow')
  call load_case(base(:9),[6,1],[9,1])
  call rejected('missing terminator', '99 0')
  rows = base
  rows(10) = '99 1'
  call load_case(rows,[6,1],[9,1])
  call rejected('malformed terminator', 'reserved')
  call load_case([character(len=256) :: base, '99 0'],[6,1],[9,1])
  call rejected('extra terminator', 'ATOM_MAP')
  call load_case([character(len=256) :: ' ! comment','', '8 2','0 4 1 0 1','1 1', &
    '0 5 1 0 1','2 1','99 0'],[8],[16])
  call require(error == '', 'F/G shells and blank/comment records')
  call xml_basis_clear()
  call xml_basis_dimensions(ns,np)
  call require(ns == 0 .and. np == 0, 'explicit clear')
  open(newunit=unit,file=trim(filename),status='old',iostat=ios)
  if (ios == 0) close(unit,status='delete')
  write(*,'(a,i0,a)') 'CRYSTAL XML companion parser: ',checks,' checks passed.'
contains
  subroutine load_case(lines,z,ao)
    character(len=*), intent(in) :: lines(:)
    integer, intent(in) :: z(:),ao(:)
    integer :: out,i
    open(newunit=out,file=trim(filename),status='replace',action='write')
    do i=1,size(lines)
      write(out,'(a)') trim(lines(i))
    end do
    close(out)
    call xml_basis_load(trim(filename),z,ao,error)
  end subroutine
  subroutine require(condition,message)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: message
    if (.not.condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
    checks=checks+1
  end subroutine
  subroutine rejected(name,part)
    character(len=*),intent(in) :: name,part
    call require(error /= '' .and. index(error,part)>0,trim(name)//': '//trim(error))
    call xml_basis_dimensions(ns,np)
    call require(ns==0 .and. np==0,trim(name)//' clears retained data')
  end subroutine
end program
