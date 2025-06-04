! xpm.f90
!
! Interface to XPM for Fortran 2003.
!
! Author:  Philipp Engel
! Licence: ISC
module xpm
    use, intrinsic :: iso_c_binding
    use :: xlib, only: x_image
    implicit none
    private

    ! C bindings.
    public :: xpm_read_file_to_image_
    public :: xpm_read_file_to_pixmap

    ! Wrapper functions and routines.
    public :: xpm_read_file_to_image

    interface
        ! int XpmReadFileToImage(Display *display, char *filename, XImage *image_return, XImage *shapeimage_return, XpmAttributes *attributes)
        function xpm_read_file_to_image_(display, file_name, image_return, shape_image_return, attributes) &
                bind(c, name='XpmReadFileToImage')
            import :: c_char, c_int, c_ptr
            implicit none
            type(c_ptr),            intent(in),  value :: display
            character(kind=c_char), intent(in)         :: file_name
            type(c_ptr),            intent(out)        :: image_return
            type(c_ptr),            intent(out)        :: shape_image_return
            type(c_ptr),            intent(in),  value :: attributes
            integer(kind=c_int)                        :: xpm_read_file_to_image_
        end function xpm_read_file_to_image_

        ! int XpmReadFileToPixmap(Display *display, Drawable d, char *filename, Pixmap *pixmap_return, Pixmap *shapemask_return, XpmAttributes *attributes)
        function xpm_read_file_to_pixmap(display, d, file_name, pixmap_return, shapemask_return, attributes) &
                bind(c, name='XpmReadFileToPixmap')
            import :: c_char, c_int, c_long, c_ptr
            implicit none
            type(c_ptr),            intent(in), value :: display
            integer(kind=c_long),   intent(in), value :: d
            character(kind=c_char), intent(in)        :: file_name
            integer(kind=c_long),   intent(out)       :: pixmap_return
            integer(kind=c_long),   intent(out)       :: shapemask_return
            type(c_ptr),            intent(in), value :: attributes
            integer(kind=c_int)                       :: xpm_read_file_to_pixmap
        end function xpm_read_file_to_pixmap
    end interface
contains
    function xpm_read_file_to_image(display, file_name, image_return, shape_image_return, attributes)
        type(c_ptr),            intent(in)  :: display
        character(len=*),       intent(in)  :: file_name
        type(x_image), pointer, intent(out) :: image_return
        type(x_image), pointer, intent(out) :: shape_image_return
        type(c_ptr),            intent(in)  :: attributes

        type(c_ptr) :: ptr1, ptr2
        integer     :: xpm_read_file_to_image

        xpm_read_file_to_image = xpm_read_file_to_image_(display, file_name, ptr1, ptr2, attributes)

        call c_f_pointer(ptr1, image_return)
        call c_f_pointer(ptr2, shape_image_return)
    end function xpm_read_file_to_image
end module xpm
