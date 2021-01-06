    Module Policy
    implicit none
    type sparseCOOType
        integer (kind=1) :: col
        integer (kind=1) :: row
        real (kind=4):: val
    end type sparseCOOType

    type policyType
        type(sparseCOOType), allocatable :: COO(:)
    !contains
    !procedure :: write_sample => write_container_sample_impl
    !procedure :: read_sample  => read_container_sample_impl
    !
    !generic   :: write(unformatted) => write_sample
    !generic   :: read(unformatted) => read_sample
    end type policyType
    contains

    subroutine write_container_sample_impl(this, unit, iostat, iomsg)
    class(policyType), intent(in)    :: this
    integer, intent(in)         :: unit
    integer, intent(out)        :: iostat
    character(*), intent(inout) :: iomsg
    integer :: i

    write( unit) size(this%COO)
    write(unit) this%COO
    !write(unit, iostat=iostat, iomsg=iomsg) size(this%COO)
    !do i=1,size(this%COO)
    !    write(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%col
    !    write(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%row
    !    write(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%val
    !end do
    end subroutine write_container_sample_impl

    subroutine read_container_sample_impl(this, unit, iostat, iomsg)
    class(policyType), intent(inout) :: this
    integer, intent(in)         :: unit
    integer, intent(out)        :: iostat
    character(*), intent(inout) :: iomsg
    integer :: i, sizeCOO

    read(unit) sizeCOO
    !write(*,*) sizeCOO
    allocate(this%COO(sizeCOO))
    read(unit) this%COO
    !read(unit, iostat=iostat, iomsg=iomsg) this%COO
    !if (iostat /= 0 ) then
        !write(*,*) "Error"
        !do i=1,sizeCOO
        !    read(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%col
        !    read(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%row
        !    read(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%val
        !end do
        !deallocate(this%COO)
        !allocate(this%COO(1000))
        !read(unit, iostat=iostat, iomsg=iomsg) this%COO
    !end if

    !read(unit, iostat=iostat, iomsg=iomsg) sizeCOO
    !allocate(this%COO(sizeCOO))
    !
    !do i=1,sizeCOO
    !    read(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%col
    !    read(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%row
    !    read(unit, iostat=iostat, iomsg=iomsg) this%COO(i)%val
    !end do

    end subroutine read_container_sample_impl
    end module