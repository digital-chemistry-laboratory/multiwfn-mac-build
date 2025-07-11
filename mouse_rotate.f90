!!------------------- Mouse rotate module. The subroutine mouse_rotate defined in this module is used to rotate system in GUI
!Link user32.lib
!Written by Yujie Liu, 2025-04-08 (http://bbs.keinsci.com/thread-52873-1-1.html)
!Modified and greatly extended by Tian Lu, 2025-4-28 (http://bbs.keinsci.com/thread-53146-1-1.html)
!Ref.: https://learn.microsoft.com/en-us/windows/win32/inputdev/wm-lbuttondown
#ifndef disable_mouse_rotate
#ifdef _WIN32
	!Use Intel fortran ifwin
#	ifdef __INTEL_COMPILER
		module mouse_rotate_mod
		implicit none
    
		contains
			subroutine mouse_rotate(id)
				use defvar
				use plot
				use IFWIN
				integer(WORD) :: startX = 0, startY = 0 !Mouse position
				logical :: isDragging = .false.   !Mouse dragging status
				logical :: isCtrlPressed = .false.   !Ctrl key status
				logical :: isShiftPressed = .false.   !Shift key status
				integer, intent(in) :: id
				type(T_MSG) :: current_msg 
                INTEGER :: wParam
				integer(WORD) :: currentX, currentY, delx, dely
				real*8 :: yvutemp,ZVUtmp
				integer(HANDLE) :: active_window
				integer(BOOL) :: dummy1
				integer(HANDLE) :: dummy2
				character tmpstr*20
        
				active_window = GetActiveWindow()
				if (active_window == NULL) return
				dummy2 = SetCapture(active_window)
        
				do while (.true.)
					if (PeekMessage(current_msg, active_window, 0, 0, PM_REMOVE) == TRUE) then
						dummy1 = TranslateMessage(current_msg)
						dummy2 = DispatchMessage(current_msg)
                
						select case (current_msg%message)
						case (WM_LBUTTONDOWN)
							startX = LOWORD(int(current_msg%lParam,kind=4))  !Sobereva modification
							startY = HIWORD(int(current_msg%lParam,kind=4))  !Sobereva modification
							isDragging = .true.
							isCtrlPressed = (IAND(current_msg%wParam, MK_CONTROL) /= 0)
							isShiftPressed = (IAND(current_msg%wParam, MK_SHIFT) /= 0)
						case (WM_MBUTTONDOWN)
							
						case (WM_MOUSEMOVE)
							if (isDragging) then
								currentX = LOWORD(int(current_msg%lParam,kind=4))  !Sobereva modification
								currentY = HIWORD(int(current_msg%lParam,kind=4))  !Sobereva modification
								delx = currentX - startX
								dely = currentY - startY
                        
								if (abs(delx) > 2 .or. abs(dely) > 2) then
									startX = currentX
									startY = currentY
                                    
									!Plot actions
                                    IF (isCtrlPressed) THEN
										!Zoom in/out
										if (iorthoview==0) then
											ZVUtmp = ZVU + dely * 0.01D0
											if (ZVUtmp>2) ZVU=ZVUtmp
                                        else
											XFAC = XFAC - dely * 0.002D0
                                        end if
                                        !Rotate along screen
                                        camrotang = camrotang + delx * 0.2D0
                                    else if (isShiftPressed) THEN
										ORIGIN_3D_X= ORIGIN_3D_X + delx * 2
										ORIGIN_3D_Y= ORIGIN_3D_Y + dely * 2
                                    else !Rotate
										XVU = XVU + delx * 0.3
										yvutemp = YVU + dely * 0.3
										if (yvutemp >= -90D0 .and. yvutemp < 90D0) YVU = yvutemp
                                    end if
									if (GUI_mode/=2) then
										call drawmol
									else if (GUI_mode==2) then
										call drawplane(dp_init1,dp_end1,dp_init2,dp_end2,dp_init3,dp_end3)
										write(tmpstr,"(f8.2)") XVU
										call SWGTXT(idissetplaneXVU,tmpstr)
										write(tmpstr,"(f8.2)") YVU
										call SWGTXT(idissetplaneYVU,tmpstr)
									end if
								end if
							end if
                
						case (WM_LBUTTONUP)
							if (ReleaseCapture() == TRUE) isDragging = .false.
							return

						end select
					else
						call Sleep(1)
					end if
				end do
			end subroutine mouse_rotate
		end module mouse_rotate_mod

#	else !Use general C binding
		module mouse_rotate_mod
			use, intrinsic :: iso_c_binding
			implicit none
    
			type, bind(c) :: POINT
				integer(c_long) :: x, y
			end type
    
			type, bind(c) :: MSG
				integer(c_intptr_t) :: hwnd
				integer(c_int) :: message
				integer(c_intptr_t) :: wParam ! should 8 bytes on x64
				integer(c_intptr_t) :: lParam ! should 8 bytes on x64
				integer(c_int) :: time
				type(POINT) :: pt
			end type
    
			interface
				function SetCapture(hWnd) bind(c, name='SetCapture')
					use iso_c_binding
					integer(c_intptr_t) :: SetCapture
					integer(c_intptr_t), value :: hWnd
				end function
        
				function ReleaseCapture() bind(c, name='ReleaseCapture')
					use iso_c_binding
					logical(c_bool) :: ReleaseCapture
				end function
        
				function  GetActiveWindow() bind(c, name='GetActiveWindow')
					use iso_c_binding
					integer(c_intptr_t) :: GetActiveWindow
				end function
        
				function PeekMessage(lpMsg, hWnd, wMsgFilterMin, wMsgFilterMax, wRemoveMsg) bind(c, name='PeekMessageA')
					use iso_c_binding
					import :: MSG
					logical(c_bool) :: PeekMessage
					type(MSG) :: lpMsg
					integer(c_intptr_t), value :: hWnd
					integer(c_int), value :: wMsgFilterMin, wMsgFilterMax, wRemoveMsg
				end function
        
				function TranslateMessage(lpMsg) bind(c, name='TranslateMessage')
					use iso_c_binding
					import :: MSG
					logical(c_bool) :: TranslateMessage
					type(MSG) :: lpMsg
				end function
        
				function DispatchMessage(lpMsg) bind(c, name='DispatchMessageA')
					use iso_c_binding
					import :: MSG
					integer(c_intptr_t) :: DispatchMessage
					type(MSG) :: lpMsg
				end function
        
				subroutine Sleep(dwMilliseconds) bind(c, name='Sleep')
					use iso_c_binding
					integer(c_int), value :: dwMilliseconds
				end subroutine
			end interface
    
			! message codes
			integer, parameter :: WM_LBUTTONDOWN = 513  !0x0201
			integer, parameter :: WM_LBUTTONUP = 514    !0x0202
			integer, parameter :: WM_MOUSEMOVE = 512    !0x0200
			integer, parameter :: WM_KILLFOCUS = 8      !0x0008
			integer, parameter :: PM_REMOVE = 1         !0x0001
    
		contains
    
			subroutine mouse_rotate(id)
				use defvar
				use plot
				integer :: startX = 0, startY = 0 ! mouse position
				logical :: isDragging = .false.   ! mouse dragging status
				integer, intent(in) :: id
				type(MSG) :: current_msg 
				integer :: currentX, currentY, delx, dely
				real*8 :: yvutemp
				integer(c_intptr_t) :: active_window
				logical(c_bool) :: dummy1
				integer(c_intptr_t) :: dummy2
				character tmpstr*20
        
				active_window = GetActiveWindow()
				if (active_window == 0) return
				dummy2 = SetCapture(active_window)
        
				do while (.true.)
					if (PeekMessage(current_msg, active_window, 0, 0, PM_REMOVE)) then
						dummy1 = TranslateMessage(current_msg)
						dummy2 = DispatchMessage(current_msg)
                
						select case (current_msg%message)
						case (WM_LBUTTONDOWN)
							startX = iand(current_msg%lParam, 65535_8)
							startY = iand(ishft(current_msg%lParam, -16), 65535_8)
							isDragging = .true.
						case (WM_MOUSEMOVE)
							if (isDragging) then
								currentX = iand(current_msg%lParam, 65535_8)
								currentY = iand(ishft(current_msg%lParam, -16), 65535_8)
								delx = currentX - startX
								dely = currentY - startY
                        
								if (abs(delx) > 2 .or. abs(dely) > 2) then
									startX = currentX
									startY = currentY

									!Plot actions
									XVU = XVU + delx * 0.3
									yvutemp = YVU + dely * 0.3
									if (yvutemp >= -90D0 .and. yvutemp < 90D0) YVU = yvutemp
									if (GUI_mode/=2) then
										call drawmol
									else if (GUI_mode==2) then
										call drawplane(dp_init1,dp_end1,dp_init2,dp_end2,dp_init3,dp_end3)
										write(tmpstr,"(f8.2)") XVU
										call SWGTXT(idissetplaneXVU,tmpstr)
										write(tmpstr,"(f8.2)") YVU
										call SWGTXT(idissetplaneYVU,tmpstr)
									end if
								end if
							end if
                
						case (WM_LBUTTONUP)
							if (ReleaseCapture()) isDragging = .false.
							return

						end select
					else
						call Sleep(1)
					end if
				end do
			end subroutine mouse_rotate
		end module mouse_rotate_mod
#	endif
#else
	! Use X11 & fortran-xlib
	module mouse_rotate_mod
		use, intrinsic :: iso_c_binding
		use xlib
		implicit none

		! key code
		integer(c_long), parameter :: XK_Control_L  = 65507     !0xffe3
		integer(c_long), parameter :: XK_Control_R  = 65508     !0xffe4
		integer(c_long), parameter :: XK_Shift_L    = 65505     !0xffe1
		integer(c_long), parameter :: XK_Shift_R    = 65506     !0xffe2
		integer(c_long), parameter :: None          = 0     
		integer(c_int), parameter :: XC_exchange    = 50        ! exchange shape
		integer(c_int), parameter :: XC_hand2       = 60        ! hand shape
		integer(c_int), parameter :: GrabModeSync   = 0
		integer(c_int), parameter :: GrabModeAsync  = 1
		integer(c_int), parameter :: AnyModifier    = ishft(1, 15)
		integer(c_int), parameter :: AnyButton      = 0
		integer(c_int), parameter :: AnyKey         = 0

		interface
			function XLookupKeysym(xkey_event, index) bind(C, name="XLookupKeysym")
				import :: c_int, c_long, x_key_event
				type(x_key_event), intent(in) :: xkey_event  
				integer(c_int), value :: index               
				integer(c_long) :: XLookupKeysym             
			end function

			function XGrabButton(display, button, modifiers, grab_window, owner_events, &
				event_mask, pointer_mode, keyboard_mode, confine_to, cursor) &
				bind(c, name="XGrabButton")
				import :: c_ptr, c_int, c_long, c_bool
				type(c_ptr), value :: display
				integer(c_int), value :: button
				integer(c_int), value :: modifiers
				integer(c_long), value :: grab_window
				logical(c_bool), value :: owner_events
				integer(c_int), value :: event_mask
				integer(c_int), value :: pointer_mode
				integer(c_int), value :: keyboard_mode
				integer(c_long), value :: confine_to
				integer(c_long), value :: cursor
				integer(c_int) :: XGrabButton
			end function XGrabButton

			function XGrabKey(display, keycode, modifiers, grab_window, owner_events, &
				pointer_mode, keyboard_mode) bind(C, name='XGrabKey')
				import :: c_int, c_ptr, c_long, c_bool
				type(c_ptr), value :: display
				integer(c_int), value :: keycode, modifiers
				integer(c_long), value :: grab_window
				logical(c_bool), value :: owner_events
				integer(c_int), value :: pointer_mode, keyboard_mode
				integer(c_int) :: XGrabKey
			end function

			function XCreateFontCursor(display, shape) bind(c, name="XCreateFontCursor")
				import :: c_int, c_ptr, c_long
				type(c_ptr), value :: display
				integer(c_int), value :: shape
				integer(c_long) :: XCreateFontCursor
			end function

			subroutine XDefineCursor(display, window, cursor) bind(c, name="XDefineCursor")
				import :: c_int, c_ptr, c_long
				type(c_ptr), value :: display
				integer(c_long), value :: window, cursor
			end subroutine

			subroutine XUndefineCursor(display, window) bind(c, name="XUndefineCursor")
				import :: c_int, c_ptr, c_long
				type(c_ptr), value :: display
				integer(c_long), value :: window
			end subroutine XUndefineCursor

			subroutine XFreeCursor(display, cursor) bind(c, name="XFreeCursor")
				import :: c_int, c_ptr, c_long
				type(c_ptr), value :: display
				integer(c_long), value :: cursor
			end subroutine XFreeCursor
		end interface

	contains
		subroutine mypot()
			use plot
			character tmpstr*20

			if (GUI_mode/=2) then
				call drawmol
			else if (GUI_mode==2) then
				call drawplane(dp_init1,dp_end1,dp_init2,dp_end2,dp_init3,dp_end3)
				write(tmpstr,"(f8.2)") XVU
				call SWGTXT(idissetplaneXVU,tmpstr)
				write(tmpstr,"(f8.2)") YVU
				call SWGTXT(idissetplaneYVU,tmpstr)
			end if
		end subroutine mypot
    
		subroutine mouse_rotate(id)
			use defvar
			integer, intent(in) :: id
			type(c_ptr) :: display
			integer(c_long) :: hwnd, cursor
			integer(c_long)  :: mask = NO_EVENT_MASK
			type(x_event) :: event
			integer :: ival
			integer(c_int) :: currx, curry, delx, dely, startx, starty
			logical :: isDragging = .false.
			logical :: isCtrlPressed = .false.
			logical :: isShiftPressed = .false.
			real*8 :: yvutemp, ZVUtmp
			integer(c_long) :: keycode = 0
			integer(c_int) :: status
			logical(c_bool) :: owner_events = .true.

			! check gui type must be X11(1=OpenMotif, 2=GTK, 3=Win32)
			call gwggui(ival)
			if (ival/=1) return
        
			! connect X server
			display = x_open_display(c_null_char)
			if (.not. c_associated(display)) return
			! get window id (X11) for widget idisgraph
			call gwgxid(id, ival)
			hwnd = ival

			! set mouse event listening
			mask = ior(BUTTON_PRESS_MASK, BUTTON_RELEASE_MASK)
			mask = ior(mask, BUTTON1_MOTION_MASK) ! left button

			status = XGrabButton(display, AnyButton, AnyModifier, hwnd, owner_events, &
								int(mask, 4), GrabModeAsync, GrabModeAsync, None, None)
			!print *, 'GrabButton status: ', status

			status = XGrabKey(display, AnyKey, AnyModifier, hwnd, owner_events, &
							GrabModeAsync, GrabModeAsync)
			!print *, 'XGrabKey status: ', status
        
			! change cursor use hand
			cursor = XCreateFontCursor(display, XC_hand2)
			call XDefineCursor(display, hwnd, cursor)

			! event loop
			do while (.true.)
				!print *, 'waiting for event...'
				call x_next_event(display, event)
				!print *, 'event type= ', event%type
				select case (event%type)
					case (BUTTON_PRESS)
						!print *, 'button press'
						isDragging = .true.
						startx = event%x_button%x
						starty = event%x_button%y
						!print *, 'startx=', startx, 'starty=', starty

						! change cursor to exchange shape
						call XUndefineCursor(display, hwnd)
						call XFreeCursor(display, cursor)
						cursor = XCreateFontCursor(display, XC_exchange)
						call XDefineCursor(display, hwnd, cursor)

					case (KEY_PRESS)
						keycode = XLookupKeysym(event%x_key, 0)
						!print  *, "keycode= ", keycode

						! ctrl and shift key
						if (keycode ==  XK_Shift_L .or. keycode == XK_Shift_R) isShiftPressed = .true.
						if (keycode ==  XK_Control_L .or. keycode == XK_Control_R) isCtrlPressed = .true.

					case (KEY_RELEASE, BUTTON_RELEASE)
						! print *, 'key/button release'
						isDragging = .false.
						isShiftPressed = .false.
						isCtrlPressed = .false.
						! delete all and close server
						call XUndefineCursor(display, hwnd)
						call XFreeCursor(display, cursor)
						call x_close_display(display) ! close display
						return

					case (MOTION_NOTIFY)
						if (isDragging .eqv. .false.) cycle
						currx = event%x_motion%x
						curry = event%x_motion%y
						delx = currx - startx
						dely = curry - starty
						if (abs(delx) > 2 .or. abs(dely) > 2) then
							startX = currx
							startY = curry
                        
							if (isCtrlPressed) then
								!Zoom in/out
								if (iorthoview==0) then
									ZVUtmp = ZVU + dely * 0.01D0
									if (ZVUtmp>2) ZVU=ZVUtmp
								else
									XFAC = XFAC - dely * 0.002D0
								end if
								!Rotate along screen
								camrotang = camrotang + delx * 0.2D0
							else if (isShiftPressed) then
								ORIGIN_3D_X= ORIGIN_3D_X + delx * 2D0
								ORIGIN_3D_Y= ORIGIN_3D_Y + dely * 2D0
							else !Rotate
								XVU = XVU + delx * 0.3D0
								yvutemp = YVU + dely * 0.3D0
								if (yvutemp >= -90D0 .and. yvutemp < 90D0) YVU = yvutemp
							end if

							! plot 
							call mypot()
						endif
				end select
			end do
		end subroutine mouse_rotate
	end module mouse_rotate_mod
#endif
#endif
