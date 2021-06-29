!======================================================================!
 subroutine calc_distance( na, m, nb, xa, xb, Dmat, same, dist_measure)!
!                                                                      !
!This subroutine computes the distance matrix between two datasets     !
!----------------------------------------------------------------------!
!INPUTS:                                                               !
!           na: number of rows in matrix xa                            !
!           nb: number of rows in matrix xb                            !
!           m : number of columns in matrices xa and xb                !
!           xa: data matrix a                                          !
!           xb: data matrix b                                          !
!         same: xa == xb?                                              !
! dist_measure: "cb" for cityblock or "eu" for Euclidean               !
!OUTPUT:                                                               !
!       Dmat: (na, nb) distance matrix                                 !
!----------------------------------------------------------------------!

    implicit none

    integer, intent( in)                :: na, nb, m
    real*8 , intent( in)                :: xa( na,m), xb( nb,m)
    logical, intent( in)                :: same
    character( len=2), intent( in)      :: dist_measure
    real*8 , intent(out)                :: Dmat( na,nb)

    integer     :: i, j
    real*8      :: dvec( m)

    interface 
        function pair_distance( n, a, b, dist_measure) result( d)          !
            integer, intent( in)                :: n
            real*8 , intent( in)                :: a( n), b( n)
            character( len=2), intent( in)      :: dist_measure
            real*8                              :: d
        end function pair_distance
    end interface

! interface for f2py
!f2py intent( in)       :: na, nb, m, xa, xb, same
!f2py intent(out)       :: Dmat

    Dmat(:,:) = 0.d0
    if ( same ) then ! all( xa==xb)
        do i=1  ,na
        do j=i+1,nb
            !dvec      = xa(i,:) - xa(j,:)
            !Dmat(i,j) = dsqrt( dot_product( dvec, dvec))
            Dmat(i,j) = pair_distance( m, xa(i,:), xb(j,:), dist_measure)
            Dmat(j,i) = Dmat(i,j)
        enddo
        enddo
    else
        do i=1,na
        do j=1,nb
            !dvec      = xa(i,:) - xb(j,:)
            !Dmat(i,j) = dsqrt( dot_product( dvec, dvec))
            Dmat(i,j) = pair_distance( m, xa(i,:), xb(j,:), dist_measure)
        enddo
        enddo
    endif

 end subroutine calc_distance
!-------------------------------------------------------------------!


!===================================================================!
 function pair_distance( n, a, b, dist_measure) result( d)          !
!                                                                   !
!This function computes the distance matrix between two points      !
!-------------------------------------------------------------------!
!INPUTS:                                                            !
!       n    : number of elements in array a                        !
!       a    : first  array of coordinates                          !
!       b    : second array of coordinates                          !
!dist_measure: name of distance type                                !
!OUTPUT:                                                            !
!       d : distance                                                !
!-------------------------------------------------------------------!
    integer, intent( in)                :: n
    real*8 , intent( in)                :: a( n), b( n)
    character( len=2), intent( in)     :: dist_measure
    real*8                              :: d

    select case ( dist_measure )
        case ('cb')
            d = sum( dabs( a-b))
        case default
            d = dsqrt( dot_product( a-b, a-b))
    end select

 end function pair_distance


program test_dist

integer         :: i,j,k
integer, parameter      :: na=3, nb=4, m=2
real*8          :: xa(na, m), xb(nb, m), Dmat(na,nb)

do j=1,size( xa,2)
    do i=1,size( xa,1)
        call random_number( xa(i,j))
    enddo

    do i=1,size( xb,1)
        call random_number( xb(i,j))
    enddo
enddo

do i=1,size(xa,1)
    print*,xa(i,:)
enddo
print*
do i=1,size(xb,1)
    print*,xb(i,:)
enddo

call calc_distance( na, m, nb, xa, xb, Dmat, .false., 'cb') 

print*
do i=1,na
    print*,Dmat(i,:)
enddo

end program test_dist

