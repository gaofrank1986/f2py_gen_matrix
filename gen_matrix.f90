
module gen_matrix
    
    use mesh
    implicit none
    
    real(8),allocatable :: amata(:,:,:),bmata(:,:)
    real(8),parameter :: ex(4) = (/1.,1.,-1.,-1./)
    real(8),parameter :: ey(4) = (/1.,-1.,-1.,1./)
contains
    
    subroutine test()
            print *, ex
    end subroutine
    
    subroutine calc_matrix()
        implicit none
       
        integer :: inode,ii,ielem,i
        real(8) :: xp,yp,zp,s_angle 
        real(8) :: angle(nnode)
    
        allocate(amata(nnode,nnode,nsys))
        allocate(bmata(nnode,nsys))
        
        amata = 0
        bmata = 0        
        
        do inode = 1,nnf
            XP=XYZ(1,INODE)
            YP=XYZ(2,INODE)
            ZP=XYZ(3,INODE)
            CALL SOLIDANGLE(INODE,NNODE,NELEM,&
            &NCN,NCON,NODQUA,H,XYZ,DXYZE,S_ANGLE)    
            S_ANGLE=1.0d0-S_ANGLE
            angle(inode) = s_angle
            amata(inode,inode,1:nsys) = s_angle

            do ielem = 1,nelemf

                II=0   
                DO I=1, NODNOE(INODE)
                    IF(IELEM .EQ. NODELE(INODE,I)) THEN
                        II=II+1
                    ENDIF
                ENDDO

 
            end do! ielem
        end do! inode   
    end subroutine calc_matrix        
end module
