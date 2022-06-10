        module libs 
          use iso_c_binding 
          
          abstract interface
            function getfunval_proc(r1) result(funval) bind(c)
              import :: c_double 
              real(c_double), intent(in), value :: r1 
              real(c_double)                    :: funval 
            end function 
          end interface 

          contains 

          subroutine integrate(getfunvalfromc, x0, x1, n, outputreal) bind(c, name='Integrate')
            type(c_funptr), intent(in), value  :: getfunvalfromc 
            real(c_double), intent(in)         :: x0, x1 
            integer(c_int), intent(in)         :: n 
            real(c_double), intent(out)        :: outputreal
            procedure(getfunval_proc), pointer :: getfunval
            real*8                             :: res
            call c_f_procpointer(cptr=getfunvalfromc, fptr=getfunval)
            call integr_method(getfunval, n, x0, x1, res)
            outputreal = res
            nullify(getfunval)
          end subroutine  
  
          subroutine get_determinant(n, a, det) bind(c, name='GetDeterminant')
            integer(c_int), intent(in)  :: n 
            real(c_double), intent(in)  :: a(n, n)
            real(c_double), intent(out) :: det
            det = deter(n, a) 
          end subroutine 

          subroutine get_fft(n, x) bind(c, name='FastFourierTransform')
            integer(c_int), intent(in)                             :: n 
            complex(c_double_complex), intent(inout), dimension(n) :: x
            call fft(n, x)
          end subroutine 

          subroutine least_square_fit(n, x, y, m, q, r2, sem, seq) bind(c, name='LeastSquareFit')
            integer(c_int), intent(in)                :: n 
            real(c_double), intent(in), dimension(n)  :: x, y 
            real(c_double), intent(out)               :: m, q, r2, sem, seq
            real*8                                    :: ssxx, ssyy, ssxy, s, xm, ym 
            xm = sum(x(:)) / dble(n)
            ym = sum(y(:)) / dble(n)
            ssxx = sum((x(:) - xm)**2)
            ssyy = sum((y(:) - ym)**2)
            ssxy = sum((x(:) - xm)*(y(:) - ym))
            m = ssxy / ssxx 
            q = ym - m*xm 
            r2 = ssxy**2 / (ssxx*ssyy)
            s = dsqrt((ssyy - m*ssxy)/dble(n-2))
            sem = s / dsqrt(ssxx)
            seq = s * dsqrt(1.d0/dble(n) + xm**2/ssxx)
          end subroutine 

          subroutine jacobi_diag(n, a, w) bind(c, name='JacobiDiagonalization')
            integer(c_int), intent(in)                     :: n 
            real(c_double), intent(inout), dimension(n, n) :: a
            real(c_double), intent(out), dimension(n)      :: w
            real*8, dimension(n)                           :: d, e
            call tred2(n, a, d, e)
            call tqli(n, d, e, a)
            w = d
          end subroutine 

          subroutine gaussjordan(n, m, a, b) bind(c, name='GaussJordan')
          !!! Pag. 1014 Numerical Recipes in Fortran 90
            integer(c_int), intent(in)                    :: n, m 
            real(c_double), intent(inout), dimension(n,n) :: a 
            real(c_double), intent(inout), dimension(n,m) :: b 
            call gaussj(n, m, a, b)
          end subroutine 

          subroutine gaussj(n, m, a, b)
          !!! Pag. 1014 Numerical Recipes in Fortran 90
            integer, intent(in)                   :: n, m 
            real*8, intent(inout), dimension(n,n) :: a 
            real*8, intent(inout), dimension(n,m) :: b 
            integer, dimension(n)                 :: ipiv, indxr, indxc 
            logical, dimension(n)                 :: lpiv 
            real*8                                :: pivinv 
            real*8, dimension(n)                  :: dumc
            integer, target, dimension(2)         :: irc
            integer, pointer                      :: irow, icol 
            irow => irc(1)
            icol => irc(2)
            ipiv = 0 
            do i=1,n
              lpiv = (ipiv == 0)
              irc = maxloc(array=abs(a), mask=outerand(lpiv, lpiv))
              ipiv(icol) = ipiv(icol) + 1
              if(ipiv(icol) .gt. 1) then
                print '(a)', 'nrerror: gaussj: singular matrix (1)' 
              endif  
              if(irow .ne. icol) then
                call swap(a(irow,:), a(icol,:))
                call swap(b(irow,:), b(icol,:))
              endif 
              indxr(i) = irow 
              indxc(i) = icol 
              if(a(icol,icol) .eq. .0d0) then 
                print '(a)', 'nrerror: gaussj: singular matrix (2)'
              endif 
              pivinv = 1.d0 / a(icol, icol)
              a(icol, icol) = 1.d0 
              a(icol, :) = a(icol, :) * pivinv 
              b(icol, :) = b(icol, :) * pivinv 
              dumc = a(:, icol)
              a(:, icol) = .0d0 
              a(icol, icol) = pivinv 
              a(:icol-1, :) = a(:icol-1, :) - outerprod(dumc(:icol-1), a(icol,:))
              b(:icol-1, :) = b(:icol-1, :) - outerprod(dumc(:icol-1), b(icol,:))
              a(icol+1:, :) = a(icol+1:, :) - outerprod(dumc(icol+1:), a(icol,:))
              b(icol+1:, :) = b(icol+1:, :) - outerprod(dumc(icol+1:), b(icol,:))
            enddo   
            do l=n,1,-1
              call swap(a(:, indxr(l)), a(:, indxc(l))) 
            enddo 
          end subroutine 

          function outerprod(a, b)
            real*8, dimension(:), intent(in)    :: a, b 
            real*8, dimension(size(a), size(b)) :: outerprod
            outerprod = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a)) 
          end function 

          function outerand(a, b)
            logical, dimension(:), intent(in)    :: a, b 
            logical, dimension(size(a), size(b)) :: outerand
            outerand = spread(a, dim=2, ncopies=size(b)) .and. spread(b, dim=1, ncopies=size(a)) 
          end function

          subroutine swap(a, b)
            real*8, dimension(:), intent(inout) :: a, b 
            real*8, dimension(size(a))          :: dum 
            dum = a 
            a = b 
            b = dum 
          end subroutine 

          subroutine swap0(a, b)
            real*8, intent(inout) :: a, b 
            real*8                :: dum
            dum = a 
            a = b 
            b = dum  
          end subroutine 

          subroutine eigsrt(n, d, v)
            integer, intent(in)                    :: n 
            real*8, dimension(n), intent(inout)    :: d
            real*8, dimension(n, n), intent(inout) :: v
            do i=1, n-1
              j = imaxloc(d(i:n)) + i - 1 
              if(j.ne.i) then
                call swap0(d(i), d(j))
                call swap(v(:, i), v(:, j)) 
              endif  
            enddo   
          end subroutine 

          function imaxloc(arr)
            real*8, intent(in), dimension(:) :: arr
            integer imaxloc 
            integer, dimension(1) :: imax 
            imax = maxloc(arr)
            imaxloc = imax(1)
          end function 

          recursive subroutine fft(n, x)
            integer*4, intent(in)                   :: n 
            complex*16, intent(inout), dimension(n) :: x
            complex*16, allocatable, dimension(:)   :: even, odd
            complex*16                              :: t 
            real*8                                  :: pi 
            if(n.le.1) return 
            allocate(odd((n+1)/2))
            allocate(even(n/2))
            pi = dacos(-1.d0)
            odd = x(1:n:2)
            even = x(2:n:2)
            call fft((n+1)/2, odd)
            call fft(n/2, even)
            do i=1, n/2
              t = exp(dcmplx(.0d0, -2.d0*pi*dble(i-1)/dble(n)))*even(i)
              x(i) = odd(i) + t
              x(i+n/2) = odd(i) - t  
            enddo 
            deallocate(odd)
            deallocate(even)
          end subroutine 

          recursive function deter(n, mat) result(accum)
            integer :: n 
            real*8, dimension(n, n)     :: mat 
            real*8                      :: accum 
            real*8, dimension(n-1, n-1) :: submat
            integer                     :: sgn
            if(n.eq.1) then
              accum = mat(1,1) 
            else 
              accum = .0d0 
              sgn = 1 
              do i=1, n
                submat(:, :i-1) = mat(2:, :i-1)
                submat(:, i:) = mat(2:, i+1:)
                accum = accum + sgn*mat(1, i)*deter(n-1, submat)
                sgn = -sgn   
              enddo 
            endif   
          end function 

          subroutine integr_method(getfunval, n, x0, x1, outputreal)
            procedure(getfunval_proc), pointer, intent(in) :: getfunval
            integer, intent(in)                :: n 
            real*8, intent(in)                 :: x0, x1 
            real*8, intent(out)                 :: outputreal
            real*8                             :: s, h, xi 
            integer*4                          :: m             
            m = n 
            if((m/2)*2 .ne. m) m = m + 1 
            s = .0d0 
            h = (x1-x0) / dble(m)
            do i=2, m-2, 2
              xi = x0 + dble(i)*h 
              s = s + 2.d0 * getfunval(xi) + 4.d0 * getfunval(xi+h)
            enddo 
            outputreal = (s + getfunval(x0) + getfunval(x1) + 4.d0*getfunval(x0+h)) * h/3.d0            
          end subroutine 

          SUBROUTINE tred2(n,a,d,e)
            IMPLICIT NONE
            INTEGER*4, INTENT(IN) :: n 
            REAL*8, DIMENSION(:,:), INTENT(INOUT) :: a
            REAL*8, DIMENSION(:), INTENT(OUT) :: d,e
            INTEGER  :: i,j,l
            REAL*8 :: f,g,h,hh,scale
            REAL*8, DIMENSION(size(a,1)) :: gg
            LOGICAL :: yesvec=.true.
            do i=n,2,-1
              l=i-1
              h=0.0
              if (l > 1) then
                scale=sum(abs(a(i,1:l)))
                if (scale == 0.0) then
                  e(i)=a(i,l)
                else
                  a(i,1:l)=a(i,1:l)/scale
                  h=sum(a(i,1:l)**2)
                  f=a(i,l)
                  g=-sign(sqrt(h),f)
                  e(i)=scale*g
                  h=h-f*g
                  a(i,l)=f-g
                  if (yesvec) a(1:l,i)=a(i,1:l)/h
                  do j=1,l
                    e(j)=(dot_product(a(j,1:j),a(i,1:j))+dot_product(a(j+1:l,j),a(i,j+1:l)))/h
                  end do
                  f=dot_product(e(1:l),a(i,1:l))
                  hh=f/(h+h)
                  e(1:l)=e(1:l)-hh*a(i,1:l)
                  do j=1,l
                    a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
                  end do
                end if
              else
                e(i)=a(i,l)
              end if
              d(i)=h
            end do
            if (yesvec) d(1)=0.0
            e(1)=0.0
            do i=1,n
              if (yesvec) then
                l=i-1
                if (d(i) /= 0.0) then
                  gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
                  a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
                end if
                d(i)=a(i,i)
                a(i,i)=1.0
                a(i,1:l)=0.0
                a(1:l,i)=0.0
              else
                d(i)=a(i,i)
              end if
            end do
          END SUBROUTINE

          SUBROUTINE tqli(n,d,e,z)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n 
            REAL*8, DIMENSION(:), INTENT(INOUT) :: d,e
            REAL*8, DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
            INTEGER :: i,iter,l,m,ndum
            REAL*8 :: b,c,dd,f,g,p,r,s
            REAL*8, DIMENSION(size(e)) :: ff
            e(:)=eoshift(e(:),1)
            do l=1,n
              iter=0
              iterate: do
                do m=l,n-1
                  dd=abs(d(m))+abs(d(m+1))
                  if (abs(e(m))+dd == dd) exit
                end do
                if (m == l) exit iterate
                iter=iter+1
                g=(d(l+1)-d(l))/(2.0d0*e(l))
                r=dsqrt(1.d0 + g**2)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                  f=s*e(i)
                  b=c*e(i)
                  r=dsqrt(f**2 + g**2)
                  e(i+1)=r
                  if (r == 0.0) then
                    d(i+1)=d(i+1)-p
                    e(m)=0.0
                    cycle iterate
                  end if
                  s=f/r
                  c=g/r
                  g=d(i+1)-p
                  r=(d(i)-g)*s+2.0d0*c*b
                  p=s*r
                  d(i+1)=g+p
                  g=c*r-b
                  if (present(z)) then
                    ff(1:n)=z(1:n,i+1)
                    z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                    z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
                  end if
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
              end do iterate
            end do
          END SUBROUTINE tqli

        end module