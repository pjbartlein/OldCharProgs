      program statst  
c
c         calculates a weighted mean, standard deviation and
c             autocorrelation of charcoal influx data, and 
c             filters and passes peaks in the influx data 
c         this version assumes equally spaced observations 
c         
c         input (output from influx.for):
c             age (years)
c             average influx (char/cm2/yr)
c
c         info file contents:
c             line 1:  input filename
c             line 2:  output filename
c             line 3:  window width (number of samples, not years)
c             line 4:  peak-to-mean ratio threshold value
c             line 5:  transformation parameter
c             line 6:  window width for smoothing peak mag and freq
c             line 7:  interval for calculating frequency
c             line 8:  missing data code 
c
c         output  
c             col 1:  age(yrs)
c             col 2:  influx (char/cm2/yr)
c             col 3:  mean (background) influx (char/cm2/yr)
c             col 4:  standard deviation of influx (char/cm2/yr)
c             col 5:  ar1 (lag-1 autocorrelation)
c             col 6:  difference (in transformed units if a transformation
c                         was requested (between influx and mean (char/cm2/yr)
c             col 7:  ratio (in transformed units if a transformation was
c                         requested) between influx and mean
c             col 8:  peak (1 if under peak, 0 otherwise)
c             col 9:  start of peak (1 if at beginning of peak, 0 otherwise)
c             col 10: total influx under peak (char/cm2)
c             col 11: smoothed value of frequency of peaks (scaled to no. of
c                         peaks per interval)
c             col 12: smoothed value of total influxes under peaks
c             col 13: influx, deviations about mean 
c
      implicit none
      integer maxobs,wmax
      parameter (maxobs=4400,wmax=3000)
      real age(maxobs),influx(maxobs),tflux(maxobs),miss
      real thresh,trans,const,dist,weight(wmax)
      real w(wmax),x(wmax),y(wmax),incl(wmax)
      real wsum,sum,tmean(maxobs),ssq,denom,tvar,tstdv,scp,tratio
      real mean(maxobs),var(maxobs),stddev(maxobs),ar1(maxobs)
      real diff(maxobs),ratio(maxobs),ptotal,pflux(maxobs),dev(maxobs)
      real interval,freq(maxobs),mag(maxobs),duration,logten
      integer peak(maxobs),pstart(maxobs),use(maxobs),nuse
      integer nh,nobs,i,t,in,ii,id,width,smooth,begt,endt
      character*64 infile,outfile                                            
c 
c         open info file and read file names
c
      open (1,file='d:\fire\statst.inf')
c
      read (1,'(a64)') infile
      write (*,901) infile
  901 format (' Input file: ',a32)
      read (1,'(a64)') outfile
      write (*,902) outfile
  902 format (' Output file: ',a32)
c     
c         open input and output files
c      
      open (2,file=infile)    
      open (3,file=outfile)
c     
      i=0
      do while (.not.eof(2))
          i=i+1
          read (2,*) age(i),influx(i)
      end do
      nobs=i
      write (*,903) nobs
  903 format (' Number of observations in data set: ',i5)
      age(nobs+1)=age(nobs)+(age(nobs)-age(nobs-1))
      influx(nobs+1)=influx(nobs)
c
c          read window width, peak-to-mean ratio threshold 
c             transformation exponent, mag and freq smoother
c             width,  interval for freq calculations, and
c             missing data code (skip if missing)
c
      read (1,*) width
      write (*,904) width
  904 format (' Window width for smoothing: ',i5)
      read (1,*) thresh
      write (*,905) thresh
  905 format (' Peak-to-mean ratio threshold ',f6.2)
      read (1,*) trans
      write (*,906) trans
  906 format (' Transformation exponent: ',f6.2)
      read (1,*) smooth
      write (*,907) smooth
  907 format (' Window width for smoothing peak mag and freq: ',i5)
      read (1,*) interval
      write (*,908) interval
  908 format (' Interval for frequency calculations: ',f7.1)
      read (1,*) miss
      write (*,909) miss
  909 format (' Missing data code:  ',f10.4)      
c
c         check for missing data, and transform data if necessary
c
      const=1.
      do 10 t=1,nobs
          use(t)=0
          if (influx(t).ne.miss) use(t)=1
          if (use(t).ne.0) then
              if (trans.eq.0.0) then
                  tflux(t)=log(influx(t)+const)
              elseif (trans.eq.1.0) then
                  tflux(t)=influx(t)
              else    
                  tflux(t)=influx(t)**trans
              endif
          endif    
   10 continue 
      logten=log(10.)                     
c
c         calculate window half-width, and weights
c
      nh=int(width/2)+1
      do 11 i=1,width+1 
          dist=float((abs(nh-i)))/float(nh)
          weight(i)=(1.0-dist**3.0)**3.0  
   11 continue  
c 
c         loop over data
c
      do 20 t=1,nobs
      if (use(t).ne.0) then
c          
          in=0
          do 21 ii=1,width+1
              id=(t-nh)+ii  
              if (id.ge.1.and.id.le.nobs) then
                  in=in+1           
                  x(in)=tflux(id)
                  w(in)=weight(ii)
                  incl(in)=use(id)
              endif
   21     continue
c
c         calculate statistics
c
c         weighted mean
c
      wsum=0.0
      sum=0.0
      tmean(t)=0.0
      do 30 i=1,in
          if (incl(i).eq.1) then
              sum=sum+x(i)*w(i)
              wsum=wsum+w(i)
          endif    
   30 continue
      if (wsum.ne.0) tmean(t)=sum/wsum
c
c         variance and standard deviation
c      
      ssq=0.0
      tvar=0.0
      nuse=0
      do 32 i=1,in
          if (incl(i).eq.1) then
              nuse=nuse+1
              ssq=ssq+w(i)*(x(i)-tmean(t))*(x(i)-tmean(t))
          endif    
   32 continue
      denom=0.0
      if (nuse.ne.0) denom=(float(nuse-1)*wsum/float(nuse))
      if (denom.ne.0.0) tvar=ssq/denom
      tstdv=sqrt(tvar)
c
c         first-order autocorrelation coefficient
c     
      scp=0.0
      ar1(t)=0.0
      nuse=0
      do 33 i=1,in-1
          if (incl(i).eq.1) then
              nuse=nuse+1
              scp=scp+w(i)*(x(i+1)-tmean(t))*(x(i)-tmean(t))
          endif    
   33 continue
      denom=0.0
      if (nuse.ne.0) denom=tvar*(nuse-1) 
      if (denom.ne.0.0) ar1(t)=scp/denom 
c
c         untransform data if necessary
c
      
      if (trans.eq.0.0) then
          mean(t)=exp(tmean(t)+(tvar/2.0))-const
          var(t)=(mean(t)**2.0)*(exp(tvar)-1.0)
          stddev(t)=sqrt(var(t))
      elseif (trans.eq.1.0) then
          mean(t)=tmean(t)
          var(t)=tvar
          stddev(t)=tstdv
      else
          mean(t)=tmean(t)**(1.0/trans)
          var(t)=tvar
          stddev(t)=tstdv**(1.0/trans)
      endif
c
c         deviations about the mean
c     
      dev(t)=0.0
c      if (stddev(t).ne.0.0) then                            
c          dev(t)=(influx(t)-mean(t))/stddev(t)
c      if (tstdv.ne.0.0) then
c          dev(t)=(tflux(t)-tmean(t))/tstdv
c      endif
      if (trans.eq.0.0) then
          dev(t)=log10(influx(t))-log10(mean(t))
c          dev(t)=influx(t)-mean(t)
      else    
          dev(t)=influx(t)-mean(t)
      endif    
c     
      endif 
		if (t.eq.410) then
			do i=1,in
				write (3,800) t,i,x(i),w(i)
				write (*,800) t,i,x(i),w(i)
  800				format (2i6,2f9.4)
			end do
		end if
c
   20 continue
c
c         peak statistics
c     
      diff(nobs+1)=0.0
      ratio(nobs+1)=thresh
      do 41 t=nobs,1,-1
      if (use(t).ne.0) then
c      
      diff(t)=0.0         
      if (tflux(t).gt.tmean(t)) diff(t)=tflux(t)-tmean(t)
c
c         define peak using difference between transformed (if requested)
c         value of influx and the mean
c
c      peak(t)=0
c      if (diff(t).gt.0.0) peak(t)=1
c         
c         define peak using ratio between transformed (if requested)
c         value of influx and the mean
c
      ratio(t)=thresh
      tratio=0.0
      peak(t)=0
      if (tflux(t).gt.tmean(t)) then 
          tratio=tflux(t)/tmean(t) 
          if (tratio.gt.thresh) then 
              ratio(t)=tratio
              peak(t)=1
          endif    
      endif
c
c         note starting observation of a peak
c  
      pstart(t)=0
      if ((ratio(t+1).le.thresh).and.(ratio(t).gt.thresh)) pstart(t)=1
c      if ((diff(t+1).le.0.0).and.(diff(t).gt.0.0)) pstart(t)=1
c     
      endif
   41 continue 
c
c         calculate total influx under peak (magnitude)
c
      ptotal=0.0
      do 42 t=nobs,1,-1
      if (use(t).ne.0) then    
          pflux(t)=0.0
          if (peak(t).eq.1) then
              ptotal=ptotal+influx(t)*(age(t+1)-age(t))
          else
              pflux(t+1)=ptotal
              ptotal=0.0
          endif
      endif    
   42 continue 
      pflux(1)=0.0
      if (peak(1).eq.1) pflux(1)=ptotal   
c
c         smooth frequency and magnitudes of peaks
c
c         calculate window half-width, and weights 
c
      nh=int(smooth/2)+1
      do 51 i=1,smooth+1 
          dist=float((abs(nh-i)))/float(nh)
          weight(i)=(1.0-dist**3.0)**3.0  
   51 continue  
c 
c         loop over data
c
      do 52 t=1,nobs
      if (use(t).ne.0) then
c          
          in=0
          do 53 ii=1,smooth+1
              id=(t-nh)+ii  
              if (id.ge.1.and.id.le.nobs) then
                  in=in+1           
                  x(in)=pstart(id)
                  y(in)=pflux(id)
                  w(in)=weight(ii)
                  incl(in)=use(id)
              endif
   53     continue
          begt=(t-nh)+1
          if (begt.lt.1) begt=1
          endt=(t-nh)+smooth+1
          if (endt.gt.nobs) endt=nobs
c
c         calculate statistics
c
c         weighted proportion of obs with peaks (frequency)
c
      wsum=0.0
      sum=0.0
      freq(t)=0.0
      do 54 i=1,in
          if (incl(i).ne.0) then
              sum=sum+x(i)*w(i)
              wsum=wsum+w(i)
          endif    
   54 continue
      duration=age(endt)-age(begt)
      if (wsum.ne.0) freq(t)=in*(sum/wsum)*(interval/duration)
c
c         weighted total influx (magnitude)
c
      wsum=0.0
      sum=0.0
      mag(t)=0.0
      do 55 i=1,in
          if (incl(i).ne.0) then
              if (y(i).gt.0.0) then
                  sum=sum+y(i)*w(i)
                  wsum=wsum+w(i)
              endif
          endif        
   55 continue
      if (wsum.ne.0) mag(t)=sum/wsum  
c     
      endif
   52 continue      
c 
c         write out the results     
c      
      do 60 t=1,nobs
      if (use(t).eq.1) then
          write (3,910) age(t),influx(t),mean(t),stddev(t),ar1(t),
     +    diff(t),ratio(t),peak(t),pstart(t),pflux(t),freq(t),
     +    mag(t),dev(t)
  910     format (f8.1,3f8.3,3f8.3,2i3,f10.3,f7.1,f10.3,f8.3)
      else
          write (3,911) age(t)
  911     format (f8.1,3(7x,'*'),3(7x,'*'),2(2x,'*'),9x,'*',
     +    6x,'*',9x,'*',7x,'*')
      endif
   60 continue
c                    
      end       