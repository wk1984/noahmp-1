module noahmp_utils
  use noahmp_const, only: r4
  implicit none
  private

  public :: error
  public :: warning
  public :: message
  public :: assert
  public :: noahmp_ptable_read_var
  public :: noahmp_ptable_read_tablestr

  interface noahmp_ptable_read_var
     module procedure noahmp_ptable_read_int0d
     module procedure noahmp_ptable_read_int1d
     module procedure noahmp_ptable_read_real0d
     module procedure noahmp_ptable_read_real1d
  end interface noahmp_ptable_read_var

contains
  subroutine error(msg)
    implicit none
    character(*), intent(in) :: msg
    write(*, *) msg
    stop
  end subroutine error


  subroutine warning(msg)
    implicit none
    character(*), intent(in) :: msg
    write(*, *) msg
  end subroutine warning


  subroutine message(msg)
    implicit none
    character(*), intent(in) :: msg
    write(*, *) msg
  end subroutine message


  subroutine assert(good, msg)
    ! assert if GOOD is good.
    ! if GOOD is not good, print MSG and stop program
    implicit none
    logical, intent(in) :: good
    character(*), intent(in) :: msg
    if (.not. good) then
       write(*, *) msg
       stop
    end if
  end subroutine assert


  subroutine noahmp_ptable_read_tablestr(filename, name, tag, nrow, tblstr)
    ! read TAG-tagged table NAME from FILENAME
    ! return the number of rows (NROW) and an array of string represent the table
    implicit none
    integer, parameter :: STRLEN = 512
    integer, parameter :: iounit = 99
    character(*), intent(in) :: filename
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    integer, intent(out) :: nrow
    character(*), intent(out) :: tblstr(:)
    character(STRLEN) :: nametag
    integer :: ii
    logical :: found
    integer :: ierr

    ! nullify outputs
    nrow = 0
    do ii = 1, size(tblstr)
       tblstr(ii) = ' '
    end do

    ! open file and read the table
    open(iounit, file=filename, form='formatted', status='old', action='read', iostat=ierr)
    call assert(ierr == 0, 'Unable to open file parameter file '//trim(filename))

    call noahmp_ptable_nametag(name, tag, nametag)
    call noahmp_ptable_find_name_tag(iounit, name, tag, found)
    call assert(found, 'Unable to find '//trim(nametag)//' in file '//trim(filename))
    read(iounit, *, iostat=ierr) nrow
    call assert(ierr == 0, &
         &      'Unable to read the number of rows for table '//trim(nametag)//' in file '//trim(filename))
    call assert(nrow <= size(tblstr), &
         &      'Table '//trim(nametag)//' in file '//trim(filename)//' is too big to be read')
    do ii = 1, nrow
       read(iounit, '(A)', iostat=ierr) tblstr(ii)
       call assert(ierr == 0, &
            &      'Unable to read table '//trim(nametag)//' from file '//trim(filename))
    end do

    close(iounit)
  end subroutine noahmp_ptable_read_tablestr


  subroutine noahmp_ptable_read_int0d(filename, name, tag, var)
    implicit none
    integer, parameter :: STRLEN = 512
    integer, parameter :: iounit = 99
    character(*), intent(in) :: filename
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    integer, intent(out) :: var
    character(STRLEN) :: nametag = ' '
    logical :: found = .false.
    integer :: ierr

    open(iounit, file=filename, form='formatted', status='old', action='read', iostat=ierr)
    call assert(ierr == 0, 'Unable to open file parameter file '//trim(filename))

    call noahmp_ptable_nametag(name, tag, nametag)
    call noahmp_ptable_find_name_tag(iounit, name, tag, found)
    call assert(found, 'Unable to find '//trim(nametag)//' in file '//trim(filename))
    read(iounit, *, iostat=ierr) var
    call assert(ierr == 0, 'Unable to read table '//trim(nametag)//' from file '//trim(filename))

    close(iounit)
  end subroutine noahmp_ptable_read_int0d


  subroutine noahmp_ptable_read_int1d(filename, name, tag, var)
    implicit none
    integer, parameter :: STRLEN = 512
    integer, parameter :: iounit = 99
    character(*), intent(in) :: filename
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    integer, intent(out) :: var(:)
    character(STRLEN) :: nametag = ' '
    logical :: found = .false.
    integer :: ierr

    open(iounit, file=filename, form='formatted', status='old', action='read', iostat=ierr)
    call assert(ierr == 0, 'Unable to open file parameter file '//trim(filename))

    call noahmp_ptable_nametag(name, tag, nametag)
    call noahmp_ptable_find_name_tag(iounit, name, tag, found)
    call assert(found, 'Unable to find '//trim(nametag)//' in file '//trim(filename))
    read(iounit, *, iostat=ierr) var
    call assert(ierr == 0, 'Unable to read '//trim(nametag)//' from file '//trim(filename))

    close(iounit)
  end subroutine noahmp_ptable_read_int1d


  subroutine noahmp_ptable_read_real0d(filename, name, tag, var)
    implicit none
    integer, parameter :: STRLEN = 512
    integer, parameter :: iounit = 99
    character(*), intent(in) :: filename
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    real(r4), intent(out) :: var
    character(STRLEN) :: nametag = ' '
    logical :: found = .false.
    integer :: ierr

    open(iounit, file=filename, form='formatted', status='old', action='read', iostat=ierr)
    call assert(ierr == 0, 'Unable to open file parameter file '//trim(filename))

    call noahmp_ptable_nametag(name, tag, nametag)
    call noahmp_ptable_find_name_tag(iounit, name, tag, found)
    call assert(found, 'Unable to find '//trim(nametag)//' in file '//trim(filename))
    read(iounit, *, iostat=ierr) var
    call assert(ierr == 0, 'Unable to read '//trim(nametag)//' from file '//trim(filename))

    close(iounit)
  end subroutine noahmp_ptable_read_real0d


  subroutine noahmp_ptable_read_real1d(filename, name, tag, var)
    implicit none
    integer, parameter :: STRLEN = 512
    integer, parameter :: iounit = 99
    character(*), intent(in) :: filename
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    real(r4), intent(out) :: var(:)
    character(STRLEN) :: nametag = ' '
    logical :: found = .false.
    integer :: ierr

    open(iounit, file=filename, form='formatted', status='old', action='read', iostat=ierr)
    call assert(ierr == 0, 'Unable to open file parameter file '//trim(filename))

    call noahmp_ptable_nametag(name, tag, nametag)
    call noahmp_ptable_find_name_tag(iounit, name, tag, found)
    call assert(found, 'Unable to find '//trim(nametag)//' in file '//trim(filename))
    read(iounit, *, iostat=ierr) var
    call assert(ierr == 0, 'Unable to read '//trim(nametag)//' from file '//trim(filename))

    close(iounit)
  end subroutine noahmp_ptable_read_real1d


  subroutine noahmp_ptable_find_name_tag(iounit, name, tag, found)
    implicit none
    integer, parameter :: STRLEN = 512
    integer, intent(in) :: iounit
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    logical, intent(out) :: found

    character(STRLEN) :: strbuf = ' '
    integer :: tagloc = 0
    integer :: ierr = 0

    found = .false.

    if (len_trim(adjustl(name)) == 0) return

    rewind(iounit)
    do while (.true.)
       strbuf = ' '
       read(iounit, *, iostat=ierr) strbuf
       if (ierr /= 0) exit

       if (strbuf(1:1) == '&') then
          tagloc = scan(strbuf, '#')
          if (tagloc == 0 .and. len_trim(adjustl(tag)) == 0) then ! has tag
             if (trim(adjustl(strbuf(2:))) == trim(adjustl(name))) then
                found = .true.
             end if
          elseif (tagloc /= 0 .and. len_trim(adjustl(tag)) /= 0) then ! no tag
             if (trim(adjustl(strbuf(2:tagloc-1))) == trim(adjustl(name)) .and. &
                  & trim(adjustl(strbuf(tagloc+1:))) == trim(adjustl(tag))) then
                found = .true.
             end if
          end if
          if (found) exit
       end if
    end do
  end subroutine noahmp_ptable_find_name_tag


  subroutine noahmp_ptable_nametag(name, tag, nametag)
    implicit none
    character(*), intent(in) :: name
    character(*), intent(in) :: tag
    character(*), intent(out) :: nametag
    if (len_trim(adjustl(name)) == 0) then
       nametag = ' '
       return
    end if
    if (len_trim(adjustl(tag)) == 0) then
       nametag = trim(adjustl(name))
    else
       nametag = trim(adjustl(name)) // '#' // trim(adjustl(tag))
    end if
  end subroutine noahmp_ptable_nametag
end module noahmp_utils