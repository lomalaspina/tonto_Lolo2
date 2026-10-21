! Explicit all-electron CRYSTAL basis companion for the legacy XML reader.
! This module has no dependency on any native CRYSTAL wavefunction interface.
module crystal_xml_basis_io
  use, intrinsic :: iso_fortran_env, only: real64, iostat_end, iostat_eor
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: xml_basis_load, xml_basis_dimensions, xml_basis_copy, xml_basis_clear
  public :: xml_basis_finite_vector

  integer, parameter :: max_line = 8192, max_token = 256
  integer, parameter :: max_blocks = 4096, max_shells = 10000
  integer, parameter :: max_primitives = 1000000

  type :: input_shell
    integer :: angular = -1
    integer :: count = 0
    real(real64), allocatable :: exponent(:), coefficient(:)
  end type
  type :: input_block
    integer :: nat = 0, z = 0, count = 0
    type(input_shell), allocatable :: shell(:)
  end type

  integer, allocatable, save :: saved_atom(:), saved_angular(:)
  integer, allocatable, save :: saved_count(:), saved_first(:)
  real(real64), allocatable, save :: saved_exponent(:), saved_coefficient(:)
contains

  ! Keep IEEE validation in this native module, compiled without fast-math.
  ! In particular, callers may use it after converting source contractions
  ! into another code's normalization convention, whose factors can overflow.
  logical function xml_basis_finite_vector(values) result(ok)
    real(real64), intent(in) :: values(:)
    ok = all(ieee_is_finite(values))
  end function

  ! Release all retained data. Failed loads also leave this empty state.
  subroutine xml_basis_clear()
    if (allocated(saved_atom)) deallocate(saved_atom)
    if (allocated(saved_angular)) deallocate(saved_angular)
    if (allocated(saved_count)) deallocate(saved_count)
    if (allocated(saved_first)) deallocate(saved_first)
    if (allocated(saved_exponent)) deallocate(saved_exponent)
    if (allocated(saved_coefficient)) deallocate(saved_coefficient)
  end subroutine

  ! Return expanded, atom-resolved dimensions; both are zero before a valid load.
  subroutine xml_basis_dimensions(n_shell, n_primitive)
    integer, intent(out) :: n_shell, n_primitive
    n_shell = 0
    n_primitive = 0
    if (allocated(saved_atom)) n_shell = size(saved_atom)
    if (allocated(saved_exponent)) n_primitive = size(saved_exponent)
  end subroutine

  ! Copy into caller-owned exact-size arrays. Atom and first-primitive indices
  ! are one-based. Shells follow ORIGINAL XML atom order, then source shell
  ! order; an SP shell becomes S followed by P. angular is l (0 through 4).
  ! Coefficients are the input CRYSTAL normalized-primitive coefficients, NOT
  ! transformed to another code's Gaussian normalization convention.
  ! On an invalid array shape all arrays are zeroed and optional error is set.
  subroutine xml_basis_copy(atom_for_shell, angular, n_prim, first_prim, &
                            exponents, contractions, error)
    integer, intent(out) :: atom_for_shell(:), angular(:), n_prim(:), first_prim(:)
    real(real64), intent(out) :: exponents(:), contractions(:)
    character(len=*), intent(out), optional :: error
    integer :: ns, np
    atom_for_shell = 0
    angular = 0
    n_prim = 0
    first_prim = 0
    exponents = 0.0_real64
    contractions = 0.0_real64
    if (present(error)) error = ''
    call xml_basis_dimensions(ns, np)
    if (ns == 0 .or. size(atom_for_shell) /= ns .or. size(angular) /= ns .or. &
        size(n_prim) /= ns .or. size(first_prim) /= ns .or. &
        size(exponents) /= np .or. size(contractions) /= np) then
      if (present(error)) error = 'XML basis copy: no loaded basis or incorrect output-array dimensions'
      return
    end if
    atom_for_shell = saved_atom
    angular = saved_angular
    n_prim = saved_count
    first_prim = saved_first
    exponents = saved_exponent
    contractions = saved_coefficient
  end subroutine

  ! Read a standalone, explicit CRYSTAL general-basis block terminated by
  ! exactly "99 0". ITYB=0, LAT=0(S),1(SP),2(P),3(D),4(F),5(G), SCAL=1.
  ! ECP and standard/internal-basis shell codes are deliberately unsupported.
  ! CRYSTAL atomic identifiers encode Z=mod(NAT,100). NAT=1..199 and >=1000
  ! are allowed; ECP identifiers 200..999 are rejected. E.g. 6,106,1006 are
  ! carbon. "99 0" is reserved, so element 99 needs an alias, e.g.1099.
  ! Blank lines and # or ! comments are allowed. All records have exact field
  ! counts; extra numbers, non-finite values and zero contractions are errors.
  !
  ! Optional trailing mapping:
  !   ATOM_MAP n
  !   xml_source_atom_index NAT
  !   ... n rows ...
  ! Maps the ORIGINAL XML source atoms, before any Tonto atom permutation.
  ! Required if a present element has multiple basis blocks; otherwise unique
  ! element blocks map automatically. Each mapped atom's spherical AO count is
  ! checked against atom_ao_count (sum of 2*l+1 after splitting SP).
  !
  ! No STOP is used. error='' means success. A failure clears previous data.
  subroutine xml_basis_load(path, atom_z, atom_ao_count, error)
    character(len=*), intent(in) :: path
    integer, intent(in) :: atom_z(:), atom_ao_count(:)
    character(len=*), intent(out) :: error
    type(input_block), allocatable :: blocks(:)
    integer, allocatable :: mapping(:)
    character(len=max_token) :: tok(8)
    character(len=512) :: detail
    integer :: unit, ios, line_number, nt, nb, nat, ns, z, b, i, j, k
    integer :: nmap, atom, mapped_nat, matches, nao, total_shell, total_prim
    integer :: sh, first, last, raw_primitives
    logical :: eof, opened

    call xml_basis_clear()
    error = ''
    detail = ''
    opened = .false.
    line_number = 0
    if (size(atom_z) == 0 .or. size(atom_z) /= size(atom_ao_count)) then
      error = 'XML basis: atom atomic-number and AO-count arrays must have equal, positive size'
      return
    end if
    if (any(atom_z < 1) .or. any(atom_z > 118) .or. any(atom_ao_count <= 0)) then
      error = 'XML basis: invalid source atom atomic number or non-positive atom AO count'
      return
    end if
    allocate(blocks(max_blocks), mapping(size(atom_z)))
    mapping = 0
    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      error = 'XML basis: cannot open companion file '//trim(path)
      return
    end if
    opened = .true.
    nb = 0
    raw_primitives = 0
    do
      call next_record(unit, line_number, tok, nt, eof, detail)
      if (detail /= '') exit
      if (eof) then
        detail = 'missing final 99 0 basis terminator'
        exit
      end if
      if (nt /= 2) then
        detail = 'expected exactly NAT NSHL (or final 99 0), not a full CRYSTAL input or listing'
        exit
      end if
      if (.not. parse_integer(tok(1), nat) .or. .not. parse_integer(tok(2), ns)) then
        detail = 'NAT and NSHL must be integers'
        exit
      end if
      if (nat == 99) then
        if (ns /= 0) detail = 'NAT=99 is reserved for 99 0; use alias 1099 for element 99'
        exit
      end if
      z = atomic_number(nat)
      if (z == 0) then
        detail = 'invalid/all-electron-unsupported NAT; ECP identifiers are not accepted'
        exit
      end if
      if (ns < 1 .or. ns > max_shells) then
        detail = 'NSHL must be between 1 and 10000'
        exit
      end if
      if (nb == max_blocks) then
        detail = 'too many basis blocks (limit 4096)'
        exit
      end if
      do b = 1, nb
        if (blocks(b)%nat == nat) exit
      end do
      if (b <= nb) then
        detail = 'duplicate NAT identifier; use distinct 1000*m+Z aliases for different bases'
        exit
      end if
      nb = nb + 1
      blocks(nb)%nat = nat
      blocks(nb)%z = z
      call read_block(unit, line_number, ns, raw_primitives, blocks(nb), detail)
      if (detail /= '') exit
    end do

    if (detail == '' .and. nb == 0) detail = 'basis contains no element blocks'
    if (detail == '') then
      call next_record(unit, line_number, tok, nt, eof, detail)
      if (detail == '' .and. .not. eof) then
        if (nt /= 2 .or. trim(tok(1)) /= 'ATOM_MAP') then
          detail = 'only an optional ATOM_MAP section is allowed after 99 0'
        else if (.not. parse_integer(tok(2), nmap)) then
          detail = 'ATOM_MAP count must be an integer'
        else if (nmap /= size(atom_z)) then
          detail = 'ATOM_MAP count must equal the number of ORIGINAL XML source atoms'
        else
          do i = 1, nmap
            call next_record(unit, line_number, tok, nt, eof, detail)
            if (detail /= '') exit
            if (eof .or. nt /= 2) then
              detail = 'ATOM_MAP needs exactly atom_index NAT on each row'
              exit
            end if
            if (.not. parse_integer(tok(1), atom) .or. .not. parse_integer(tok(2), mapped_nat)) then
              detail = 'ATOM_MAP atom index and NAT must be integers'
              exit
            end if
            if (atom < 1 .or. atom > size(atom_z)) then
              detail = 'ATOM_MAP atom index is outside the ORIGINAL XML atom list'
              exit
            end if
            if (mapping(atom) /= 0) then
              detail = 'ATOM_MAP repeats an atom index'
              exit
            end if
            do b = 1, nb
              if (blocks(b)%nat == mapped_nat) exit
            end do
            if (b > nb) then
              detail = 'ATOM_MAP references a NAT with no basis block'
              exit
            end if
            if (blocks(b)%z /= atom_z(atom)) then
              detail = 'ATOM_MAP NAT atomic number differs from the XML atom'
              exit
            end if
            mapping(atom) = b
          end do
          if (detail == '') then
            call next_record(unit, line_number, tok, nt, eof, detail)
            if (detail == '' .and. .not. eof) detail = 'unexpected records after complete ATOM_MAP'
          end if
        end if
      else if (detail == '') then
        do i = 1, size(atom_z)
          matches = 0
          do b = 1, nb
            if (blocks(b)%z == atom_z(i)) then
              matches = matches + 1
              mapping(i) = b
            end if
          end do
          if (matches == 0) then
            write(detail, '(a,i0,a,i0)') 'missing basis for XML atom ', i, ' with Z=', atom_z(i)
            exit
          else if (matches > 1) then
            detail = 'multiple basis blocks for one element: an explicit ORIGINAL XML ATOM_MAP is required'
            exit
          end if
        end do
      end if
    end if

    total_shell = 0
    total_prim = 0
    if (detail == '') then
      do i = 1, size(atom_z)
        b = mapping(i)
        if (b == 0) then
          detail = 'ATOM_MAP does not cover every XML atom'
          exit
        end if
        nao = 0
        do j = 1, blocks(b)%count
          nao = nao + 2*blocks(b)%shell(j)%angular + 1
          if (total_prim > max_primitives - blocks(b)%shell(j)%count) then
            detail = 'expanded XML basis exceeds limit of 1000000 primitive entries'
            exit
          end if
          total_prim = total_prim + blocks(b)%shell(j)%count
        end do
        if (detail /= '') exit
        if (nao /= atom_ao_count(i)) then
          write(detail, '(a,i0,a,i0,a,i0)') 'XML atom ', i, ': companion spherical AO count ', nao, &
            ' differs from XML AO count ', atom_ao_count(i)
          exit
        end if
        total_shell = total_shell + blocks(b)%count
      end do
    end if
    if (opened) close(unit)
    if (detail /= '') then
      write(error, '(a,i0,a,a)') 'XML basis companion near line ', line_number, ': ', trim(detail)
      call xml_basis_clear()
      return
    end if

    allocate(saved_atom(total_shell), saved_angular(total_shell), saved_count(total_shell), &
             saved_first(total_shell), saved_exponent(total_prim), saved_coefficient(total_prim))
    sh = 0
    first = 1
    do i = 1, size(atom_z)
      b = mapping(i)
      do j = 1, blocks(b)%count
        sh = sh + 1
        k = blocks(b)%shell(j)%count
        last = first + k - 1
        saved_atom(sh) = i
        saved_angular(sh) = blocks(b)%shell(j)%angular
        saved_count(sh) = k
        saved_first(sh) = first
        saved_exponent(first:last) = blocks(b)%shell(j)%exponent
        saved_coefficient(first:last) = blocks(b)%shell(j)%coefficient
        first = last + 1
      end do
    end do
  end subroutine

  subroutine read_block(unit, line_number, ns, raw_primitives, block, error)
    integer, intent(in) :: unit, ns
    integer, intent(inout) :: line_number
    integer, intent(inout) :: raw_primitives
    type(input_block), intent(inout) :: block
    character(len=*), intent(out) :: error
    character(len=max_token) :: tok(8)
    integer :: nt, s, j, ityb, lat, ng, sh, fields
    real(real64) :: che, scal, exponent, cs, cp
    logical :: eof
    error = ''
    allocate(block%shell(2*ns))
    sh = 0
    do s = 1, ns
      call next_record(unit, line_number, tok, nt, eof, error)
      if (error /= '') return
      if (eof .or. nt /= 5) then
        error = 'expected exactly ITYB LAT NG CHE SCAL for each shell'
        return
      end if
      if (.not. parse_integer(tok(1), ityb) .or. .not. parse_integer(tok(2), lat) .or. &
          .not. parse_integer(tok(3), ng)) then
        error = 'ITYB, LAT and NG must be integers'
        return
      end if
      if (ityb /= 0) then
        error = 'only explicit general shells ITYB=0 are supported; expand standard shells first'
        return
      end if
      if (lat < 0 .or. lat > 5) then
        error = 'unsupported LAT; supported shells are S, SP, P, D, F and G (LAT=0..5)'
        return
      end if
      if (ng < 1 .or. ng > max_primitives) then
        error = 'NG must be positive and no greater than 1000000'
        return
      end if
      if (raw_primitives > max_primitives - ng) then
        error = 'companion file exceeds limit of 1000000 primitive records'
        return
      end if
      raw_primitives = raw_primitives + ng
      if (.not. parse_real(tok(4), che) .or. .not. parse_real(tok(5), scal)) then
        error = 'CHE and SCAL must be finite real numbers'
        return
      end if
      if (che < 0.0_real64) then
        error = 'CHE shell electron charge must be non-negative'
        return
      end if
      if (scal /= 1.0_real64) then
        error = 'SCAL must be 1: supply already-scaled explicit exponents, not an implicit scale factor'
        return
      end if
      sh = sh + 1
      block%shell(sh)%angular = max(0, lat - 1)
      block%shell(sh)%count = ng
      allocate(block%shell(sh)%exponent(ng), block%shell(sh)%coefficient(ng))
      fields = 2
      if (lat == 1) then
        fields = 3
        block%shell(sh+1)%angular = 1
        block%shell(sh+1)%count = ng
        allocate(block%shell(sh+1)%exponent(ng), block%shell(sh+1)%coefficient(ng))
      end if
      do j = 1, ng
        call next_record(unit, line_number, tok, nt, eof, error)
        if (error /= '') return
        if (eof .or. nt /= fields) then
          error = 'wrong primitive record: expected exponent plus one coefficient (two for SP)'
          return
        end if
        if (.not. parse_real(tok(1), exponent) .or. .not. parse_real(tok(2), cs)) then
          error = 'primitive exponent and coefficient must be finite real numbers'
          return
        end if
        if (exponent <= 0.0_real64) then
          error = 'primitive exponent must be positive'
          return
        end if
        block%shell(sh)%exponent(j) = exponent
        block%shell(sh)%coefficient(j) = cs
        if (lat == 1) then
          if (.not. parse_real(tok(3), cp)) then
            error = 'SP p coefficient must be a finite real number'
            return
          end if
          block%shell(sh+1)%exponent(j) = exponent
          block%shell(sh+1)%coefficient(j) = cp
        end if
      end do
      if (all(block%shell(sh)%coefficient == 0.0_real64)) then
        error = 'shell contraction is identically zero'
        return
      end if
      if (.not. valid_contraction(block%shell(sh))) then
        error = 'shell contraction has non-positive or unrepresentable normalized-primitive norm'
        return
      end if
      if (lat == 1) then
        if (all(block%shell(sh+1)%coefficient == 0.0_real64)) then
          error = 'SP p contraction is identically zero'
          return
        end if
        if (.not. valid_contraction(block%shell(sh+1))) then
          error = 'SP p contraction has non-positive or unrepresentable normalized-primitive norm'
          return
        end if
        sh = sh + 1
      end if
    end do
    block%count = sh
  end subroutine

  ! Validate the radial contraction, including cancellation between repeated
  ! exponents. Coefficient scaling avoids overflow in the intermediate sum;
  ! the exponent ratio avoids computing ai*aj or ai+aj at large exponents.
  ! No empirical conditioning threshold is imposed: a representable positive
  ! norm is required, while periodic linear dependence remains a later check.
  logical function valid_contraction(shell) result(ok)
    type(input_shell), intent(in) :: shell
    integer :: i, j
    real(real64) :: scale, norm, ratio, overlap, ci, cj, term
    ok = .false.
    scale = maxval(abs(shell%coefficient))
    if (scale <= 0.0_real64) return
    norm = 0.0_real64
    do i = 1, shell%count
      ci = shell%coefficient(i)/scale
      norm = norm + ci*ci
      do j = 1, i-1
        cj = shell%coefficient(j)/scale
        ratio = min(shell%exponent(i), shell%exponent(j))/max(shell%exponent(i), shell%exponent(j))
        overlap = (2.0_real64*sqrt(ratio)/(1.0_real64+ratio))**(real(shell%angular, real64)+1.5_real64)
        term = 2.0_real64*ci*cj*overlap
        norm = norm + term
      end do
    end do
    if (.not. ieee_is_finite(norm)) return
    if (norm <= 0.0_real64) return
    term = sqrt(norm)
    if (term > 1.0_real64) then
      if (scale > huge(norm)/term) return
    end if
    term = scale*term
    if (term > sqrt(huge(norm))) return
    norm = term*term
    if (.not. ieee_is_finite(norm)) return
    ok = norm > 0.0_real64
  end function

  integer function atomic_number(nat) result(z)
    integer, intent(in) :: nat
    z = 0
    if (nat >= 1 .and. (nat < 200 .or. nat >= 1000)) z = mod(nat, 100)
  end function

  ! List-directed numeric reads alone accept punctuation such as '/' or ','.
  ! First validate each complete token, then use the compiler's numeric parser.
  logical function parse_integer(token, value) result(ok)
    character(len=*), intent(in) :: token
    integer, intent(out) :: value
    integer :: j, start, ios, n
    ok = .false.
    value = 0
    n = len_trim(token)
    if (n == 0) return
    start = 1
    if (token(1:1) == '+' .or. token(1:1) == '-') start = 2
    if (n < start) return
    do j = start, n
      if (token(j:j) < '0' .or. token(j:j) > '9') return
    end do
    read(token, *, iostat=ios) value
    ok = ios == 0
  end function

  logical function parse_real(token, value) result(ok)
    character(len=*), intent(in) :: token
    real(real64), intent(out) :: value
    integer :: j, ios
    ok = .false.
    value = 0.0_real64
    if (len_trim(token) == 0) return
    do j = 1, len_trim(token)
      if (index('0123456789+-.eEdD', token(j:j)) == 0) return
    end do
    read(token, *, iostat=ios) value
    if (ios /= 0) return
    ok = ieee_is_finite(value)
  end function

  subroutine next_record(unit, line_number, tokens, count, eof, error)
    integer, intent(in) :: unit
    integer, intent(inout) :: line_number
    character(len=*), intent(out) :: tokens(:)
    integer, intent(out) :: count
    logical, intent(out) :: eof
    character(len=*), intent(out) :: error
    character(len=max_line) :: line
    integer :: ios, comment, j, start, n, pos, nread
    tokens = ''
    count = 0
    eof = .false.
    error = ''
    do
      read(unit, '(a)', advance='no', size=nread, iostat=ios) line
      if (ios == iostat_end) then
        eof = .true.
        return
      else if (ios /= 0 .and. ios /= iostat_eor) then
        error = 'I/O error while reading companion basis'
        return
      end if
      line_number = line_number + 1
      if (ios == 0) then
        error = 'companion basis record too long (limit below 8192 characters)'
        return
      end if
      comment = len_trim(line) + 1
      pos = index(line, '#')
      if (pos > 0) comment = min(comment, pos)
      pos = index(line, '!')
      if (pos > 0) comment = min(comment, pos)
      line(comment:) = ''
      n = len_trim(line)
      do j = 1, n
        if (line(j:j) == achar(9) .or. line(j:j) == achar(13)) line(j:j) = ' '
      end do
      if (len_trim(line) /= 0) exit
    end do
    j = 1
    do while (j <= n)
      if (line(j:j) == ' ') then
        j = j + 1
        cycle
      end if
      start = j
      do while (j <= n)
        if (line(j:j) == ' ') exit
        j = j + 1
      end do
      count = count + 1
      if (count > size(tokens)) then
        error = 'too many fields in companion basis record'
        return
      end if
      if (j-start > len(tokens(1))) then
        error = 'numeric token too long'
        return
      end if
      tokens(count) = line(start:j-1)
    end do
  end subroutine
end module crystal_xml_basis_io
