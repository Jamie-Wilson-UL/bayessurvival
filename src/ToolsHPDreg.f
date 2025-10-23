c=======================================================================
      subroutine hpddensregmf(nsave,npred,alpha,tint,
     &                        worksam,llower,lupper)
c=======================================================================
      implicit none 
      integer tint
      integer nsave,npred
      double precision alpha
      double precision worksam(nsave)
      double precision llower(npred)
      double precision lupper(npred)
      integer maxnsave,maxnpred
      parameter(maxnsave=30000,maxnpred=500)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxnpred)
      integer i,ii,j

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpddensreg'")
      end if   
      if(maxnpred.lt.npred)then
         call rexit("Increase 'maxnpred' in 'hpddensreg'")
      end if   

      open(unit=2,file='dppackage2.out',status='old',
     &     form='unformatted')
      do i=1,nsave
         read(2) (workm(i,j),j=1,npred)
      end do

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            worksam(i)=workm(i,ii)
         end do  
         call hpd(nsave,alpha,worksam,alow,aupp)
         if(tint.eq.1)then
            llower(ii)=alow(1)
            lupper(ii)=aupp(1)
         else
            llower(ii)=alow(2)
            lupper(ii)=aupp(2)
         end if
      end do
      close(unit=2)
      return
      end

c=======================================================================
      subroutine hpddensreg(nsave,npred,ngrid,alpha,tint,
     &                      workv1,fs,llower,lupper)
c=======================================================================
      implicit none 
      integer tint
      integer nsave,npred,ngrid 
      double precision alpha
      double precision fs(ngrid)
      double precision workv1(nsave)
      double precision llower(npred,ngrid)
      double precision lupper(npred,ngrid)
      integer maxnsave,maxngrid
      parameter(maxnsave=30000,maxngrid=500)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxngrid)
      integer i,ii,j,l

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpddensreg'")
      end if   
      if(maxngrid.lt.ngrid)then
         call rexit("Increase 'maxngrid' in 'hpddensreg'")
      end if   

      open(unit=1,file='dppackage1.out',status='old',
     &     form='unformatted')
      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               read(1) (fs(l),l=1,ngrid)
               if(ii.eq.j)then
                  do l=1,ngrid
                     workm(i,l)=fs(l)
                  end do
               end if
            end do
         end do  
         rewind(unit=1)
         do i=1,ngrid
            do j=1,nsave
               workv1(j)=workm(j,i)
            end do
            call hpd(nsave,alpha,workv1,alow,aupp)
            if(tint.eq.1)then
               llower(ii,i)=alow(1)
               lupper(ii,i)=aupp(1)
            else
               llower(ii,i)=alow(2)
               lupper(ii,i)=aupp(2)
            end if
         end do
      end do
      close(unit=1)
      return
      end 