      program memspec
      integer n,nfdt,m,nmax,mmax
      parameter (nmax=15000,mmax=300)
      integer i,j,nobs,time,year,first,last,step,beg,end,ii
      real evlmem,fdt(mmax),mse,cof(Mmax),fpe,aic,spec,lspec
      real date(nmax),data(Nmax) 
      real freq(mmax),period(mmax),x(nmax)
      character*32 infile,outfile
c
      open(3,file='\fire\mem\memspec.inf')
      read (3,900) infile
      read (3,900) outfile 
  900 format (a32)     
      open(1,file=infile)
      open(2,file=outfile) 
c 
      read (3,*) nobs,m,nfdt,delta,n,first,last,step
c   
	ii=0   
      do 10 i=nobs,1,-1
		ii=ii+1
          read(1,*) date(i),data(i) 
cc          data(i)=log(data(i)+0.001)
		date(i)=22+ii
   10 continue
c
      do 11 j=1,nfdt
          fdt(j)=0.5*j/nfdt
          freq(j)=fdt(j)/delta
          period(j)=1.0/freq(j)
   11 continue         
c     
c     
      do 18 time=first,last,step
      year=date(time)
      beg=time-(n/2)
      end=beg+n                  
      ii=0 
      do 19 i=beg,end 
          ii=ii+1
          x(ii)=data(i)  
   19 continue    
c
c      do 19 i=1,n
c          x(i)=data(i)
c   19 continue 
c 
c      do 17 m=1,50,1
c      
      call memcof(x,N,M,mse,cof) 
      fpe=float(n+m)/float(n-m)*mse
      aic=n*alog(mse)+2.0*m
c
      do 16 j=1,nfdt 
          spec=evlmem(fdt(j),cof,m,mse)   
          lspec=alog(spec)
cc          write(*,903) freq(j),period(j),spec  
cc          write(2,903) freq(j),period(j),spec   
cc  903     format (f8.4,f8.4,f8.4)
          write(*,903) year,freq(j),period(j),spec,lspec  
          write(2,903) year,freq(j),period(j),spec,lspec
  903     format (i6,f8.4,f9.2,2f9.3)     
16    continue 
c      write (*,904) year,m,aic                                                        
c      write (2,904) year,m,aic
c  904 format (i6,i4,f11.4)   
c   17 continue 
c
   18 continue
c
      END
