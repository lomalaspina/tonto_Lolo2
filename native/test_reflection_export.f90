! Production export-observation helpers; no refreshable golden output.
program test_reflection_export
  use TYPES_MODULE, only: CRYSTAL_TYPE
  use CRYSTAL_MODULE, only: F_exp_scaled_corrected_, F2_exp_scaled_corrected_, F2_sigma_scaled_corrected_
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(CRYSTAL_TYPE) :: c
  real(8) :: f(3),obs(3),sig(3),factor(3),raw(3),unc(3),expected(3),s
  integer :: i,mode
  allocate(c%xray_data)
  allocate(c%xray_data%reflections(3))
  c%xray_data%scale_factor=2d0
  c%xray_data%extinction_factor=0d0
  c%xray_data%wavelength=1d0/0.5291772108d0
  raw=[196d0,0d0,-4d0]; unc=[20d0,0.02d0,2d0]
  do i=1,3
    c%xray_data%reflections(i)%F_calc=cmplx(3d0,4d0,8)
    c%xray_data%reflections(i)%F_disp=cmplx(1d0,2d0,8)
    c%xray_data%reflections(i)%F_exp=sqrt(max(0d0,raw(i)))
    c%xray_data%reflections(i)%F_corr=c%xray_data%reflections(i)%F_exp
    c%xray_data%reflections(i)%F_sigma=1d0
    c%xray_data%reflections(i)%F_exp_set=.true.
    c%xray_data%reflections(i)%I_exp=raw(i)
    c%xray_data%reflections(i)%I_exp_set=.true.
    c%xray_data%reflections(i)%I_sigma=unc(i)
    c%xray_data%reflections(i)%theta=0.1d0*i
  end do
  ! Test no dispersion and added complex dispersion, with/without extinction.
  ! Extinction is the SAME fitted fixed multiplier for observed value/sigma.
  do mode=0,3
    c%xray_data%add_dispersion_to_F_calc=mod(mode,2)==1
    c%xray_data%remove_dispersion_from_F_exp=.false.
    c%xray_data%extinction_factor=0d0
    if(mode>=2) c%xray_data%extinction_factor=0.4d0
    do i=1,3
      s=sin(2*c%xray_data%reflections(i)%theta)
      factor(i)=2d0/(1+0.001d0*c%xray_data%extinction_factor*25d0/s)**0.25d0
    end do
    f=F_exp_scaled_corrected_(c)
    obs=F2_exp_scaled_corrected_(c)
    sig=F2_sigma_scaled_corrected_(c)
    expected=sqrt(max(0d0,raw))/factor
    if(maxval(abs(f-expected))>1d-11) error stop 'dispersion removed twice from exported F'
    if(maxval(abs(obs-raw/factor**2))>1d-11) error stop 'signed/zero observation export changed'
    if(maxval(abs(sig-unc/factor**2))>1d-11) error stop 'intensity uncertainty export changed'
    if(.not.all(ieee_is_finite(sig)).or.any(sig<=0)) error stop 'nonfinite or zero sigma export'
  end do
  ! F-only legacy export propagates sigma without dividing by a zero F.
  c%xray_data%extinction_factor=0d0
  c%xray_data%reflections%I_exp_set=.false.
  sig=F2_sigma_scaled_corrected_(c)
  if(any(abs(sig-[7d0,.25d0,.25d0])>1d-12)) error stop 'F-only sigma propagation failed'
  ! Older corrected-amplitude target must use the selected F_corr rather
  ! than removing the magnitude of a complex dispersion factor a second time.
  c%xray_data%extinction_factor=0d0
  c%xray_data%add_dispersion_to_F_calc=.false.
  c%xray_data%remove_dispersion_from_F_exp=.true.
  c%xray_data%reflections%F_corr=6d0
  f=F_exp_scaled_corrected_(c)
  obs=F2_exp_scaled_corrected_(c)
  sig=F2_sigma_scaled_corrected_(c)
  if(any(abs(f-3d0)>1d-12).or.any(abs(obs-9d0)>1d-12).or.any(abs(sig-3d0)>1d-12)) &
    error stop 'corrected-amplitude export convention changed'
  print *, 'REFLECTION EXPORT PASS: dispersion, signed/zero I, uncertainty, scale, extinction, corrected-amplitude mode'
end program

