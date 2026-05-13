      program influx
c
c         reads level, age, and charcoal concentration, and interpolates
c         to annual values, and then averages concentrations over a
c         particular interval to produce concentration values at even
c         age intervals
c
c         input data:
c             depth in meters
c             age in years
c             Tilia concentration (char/cm3)
c             Tilia influx (char/cm2/yr)--ignored
c             Tilia deposition time (yr/cm)--ignored
c
c         info file contents:
c             line 1:  input file name
c             line 2:  output file name
c             line 3:  first year, last year, output averaging interval
c             line 4:  missing concentration data code
c                         (i.e. ignore input if conc = missing)
c
      implicit none
      integer maxlvls,maxyrs
      parameter (maxlvls=5000,maxyrs=45000)
      real dp,ag,cn,tf,sp
      real depth(maxlvls),age(maxlvls),conc(maxlvls)
      real tflux(maxlvls),span(maxlvls),cmiss
      real target,intconc(0:maxyrs),intspan(0:maxyrs)
      real aveconc,avespan,intage(0:maxyrs),aveflux(0:maxyrs),intflux
      integer i,j,l,firstyr,lastyr
      integer nlvls,nmiss,step
      character*64 infile,outfile
c
c         open info file and read file names
c
      open (1,file='d:\fire\influx.inf')
c
      read (1,'(a64)') infile
      write (*,901) infile
  901 format (' Input file: ',a64)
      read (1,'(a64)') outfile
      write (*,902) outfile
  902 format (' Output file: ',a64)
c
c          read first year, last year and
c          averaging interval for influx calculations
c
      read (1,*) firstyr,lastyr,step
      write (*,903) firstyr,lastyr,step
  903 format (' First year, last year, averaging interval: ',3i6)
c
c         read missing data code
c
      read (1,*) cmiss
      write (*,904) cmiss
  904 format (' Observations with conc = ',f6.1,'  will be ignored')
c
c         open input and output files
c
      open (2,file=infile)
      open (3,file=outfile)
c
      nmiss=0
      i=0
      do while (.not.eof(2))
          read (2,*) dp,ag,cn,tf,sp
c
c         check for missing concentration value
c
          if (cn.ne.cmiss) then
              i=i+1
              depth(i)=dp
              age(i)=ag
              conc(i)=cn
              tflux(i)=tf
              span(i)=sp
          else
              nmiss=nmiss+1
          endif
c
      end do
      nlvls=i
      write (*,905) nmiss
  905 format (' Number of levels with missing data: ',i5)
      write (*,906) nlvls
  906 format (' Number of usable levels in data set: ',i5)
c
c         interpolation to produce pseudo-annual values
c         (uses hunt.for from Numerical Recipes)
c
      j=1
      do 20 i=firstyr,lastyr+step
          target=i
          call hunt(age,nlvls,target,j)
c          write (*,800) i,target,j,age(j),age(j+1)
          if (j.le.1) then
              intconc(i)=conc(1)
              intspan(i)=span(1)
          else if (j.eq.nlvls) then
              intconc(i)=conc(nlvls)
              intspan(i)=span(nlvls)
          else
              intconc(i)=conc(j)
     +         + ((target-age(j))/(age(j+1)-age(j)))*(conc(j+1)-conc(j))
              intspan(i)=span(j)
     +         + ((target-age(j))/(age(j+1)-age(j)))*(span(j+1)-span(j))
          endif
          intflux=intconc(i)/intspan(i)
c          write (*,800) i,target,j,intconc(i),intspan(i),intflux
c          write (3,800) i,target,j,intconc(i),intspan(i),intflux
  800     format (i6,f8.0,i6,3f8.3)
   20 continue
c
c         calcuate average concentrations and influxes
c             over averaging intervals
c
      do 30 l=firstyr,lastyr,step
          aveconc=0.0
          avespan=0.0
          do 31 i=l,l+step-1
              aveconc=aveconc+intconc(i)
              avespan=avespan+intspan(i)
   31     continue
          aveconc=aveconc/step
          avespan=avespan/step
c
          aveflux(l)=aveconc/avespan
          intage(l)=l
c
          write (3,911) intage(l),aveflux(l)
  911     format (f10.2,3f10.3)
   30 continue
c
      stop
      end
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D04'%1'i0(95.)61j'1..