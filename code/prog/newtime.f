      program calcnewdate
C     ===================

      implicit none

      integer   date1(5),date2(5)
      integer	iargc
      real      diff
      character*(80) arg,yychar,cdat
      character*(2)  cdate(5)
      integer	flag,nc

c     check for sufficient requested arguments
      if (iargc().ne.2) then
        print*,'USAGE: newtime date (format (YY)YYMMDD_HH(MM)) timestep'
        call exit(1)
      endif

c     read and transform input
      call getarg(1,arg)
      call lenchar(arg,nc)
      if (nc.eq.9) then
        yychar=''
        read(arg(1:2),'(i2)',err=120)date1(1)
        read(arg(3:4),'(i2)',err=120)date1(2)
        read(arg(5:6),'(i2)',err=120)date1(3)
        read(arg(8:9),'(i2)',err=120)date1(4)
        date1(5)=99
      else if (nc.eq.11) then
        yychar=arg(1:2)
        read(arg(3:4),'(i2)',err=120)date1(1)
        read(arg(5:6),'(i2)',err=120)date1(2)
        read(arg(7:8),'(i2)',err=120)date1(3)
        read(arg(10:11),'(i2)',err=120)date1(4)
        date1(5)=99
      else if (nc.eq.13) then
        yychar=arg(1:2)
        read(arg(3:4),'(i2)',err=120)date1(1)
        read(arg(5:6),'(i2)',err=120)date1(2)
        read(arg(7:8),'(i2)',err=120)date1(3)
        read(arg(10:11),'(i2)',err=120)date1(4)
        read(arg(12:13),'(i2)',err=120)date1(5)
      else
        print*,'USAGE: newtime date (format (YY)YYMMDD_HH(MM)) timestep'
        call exit(1)
      endif

      call getarg(2,arg)
      call checkchar(arg,".",flag)
      if (flag.eq.0) arg=trim(arg)//"."
      read(arg,*) diff

      call newdate_test(date1,diff,date2)

      if ((date2(1).lt.date1(1)).and.
     >    (diff.gt.0.).and.
     >    (yychar.eq.'19')) yychar='20'

      if (date2(1).lt.0) date2(1)=date2(1)+100

      if ((date2(1).gt.date1(1)).and.
     >    (diff.lt.0.).and.
     >    (yychar.eq.'20')) yychar='19'

      if (date2(1).lt.10) then
        write(cdate(1),'(a,i1)')'0',date2(1)
      else
        write(cdate(1),'(i2)')date2(1)
      endif
      if (date2(2).lt.10) then
        write(cdate(2),'(a,i1)')'0',date2(2)
      else
        write(cdate(2),'(i2)')date2(2)
      endif
      if (date2(3).lt.10) then
        write(cdate(3),'(a,i1)')'0',date2(3)
      else
        write(cdate(3),'(i2)')date2(3)
      endif
      if (date2(4).lt.10) then
        write(cdate(4),'(a,i1)')'0',date2(4)
      else
        write(cdate(4),'(i2)')date2(4)
      endif
      if (date1(5).eq.99) then
        cdate(5)=''
      else
        if (date2(5).lt.10) then
          write(cdate(5),'(a,i1)')'0',date2(5)
        else
          write(cdate(5),'(i2)')date2(5)
        endif
      endif

      cdat=trim(yychar)//cdate(1)//cdate(2)//cdate(3)//
     >     '_'//cdate(4)//cdate(5)
      write(*,'(a)')trim(cdat)

      goto 200

 120  write(*,*)
     > "*** error: date must be in format (YY)YYMMDD_HH(MM) ***"
 
 200  continue

      end

      subroutine checkchar(string,char,flag)
C     ======================================

      character*(*)	string
      character*(1)	char
      integer	n,flag

      flag=0
      do n=1,len(string)
        if (string(n:n).eq.char) then
          flag=n
          return
        endif
      enddo
      end

      subroutine lenchar(string,lstr)
C     ===============================
 
      character*(*)     string
      integer   n,lstr
 
      do n=1,len(string)
        if (string(n:n).eq."") then
          lstr=n-1
          goto 100
        endif
      enddo
 100  continue
      end

      subroutine newdate_test(date1,diff,date2)
C     ====================================
C
C     Routine calculates the new date when diff (in hours) is added to
C     date1.
C     date1	int	input	array contains a date in the form
C				year,month,day,hour,minute
C     diff	real	input	timestep to go from date1; if date1(5)
C                               is equal to 99, diff is assumed to be
C                               given in hours, otherwise in minutes
C     date2	int	output	array contains new date in the same form

      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real	diff
      logical	yearchange

      data idays/31,28,31,30,31,30,31,31,30,31,30,31/

      yearchange=.false.

      if ((mod(date1(1),4).eq.0).and.(date1(2).le.2)) idays(2)=29

      date2(1)=date1(1)
      date2(2)=date1(2)
      date2(3)=date1(3)
      if (date1(5).eq.99) then
        date2(4)=date1(4)+int(diff)
        date2(5)=0
      else
        date2(4)=date1(4)
        date2(5)=date1(5)+int(diff)
      endif

      if (date2(5).ge.60) then
        date2(4)=date2(4)+int(date2(5)/60)
        date2(5)=date2(5)-int(date2(5)/60)*60
      endif
      if (date2(5).lt.0) then
        if (mod(date2(5),60).eq.0) then
          date2(4)=date2(4)-int(abs(date2(5))/60)
          date2(5)=date2(5)+int(abs(date2(5))/60)*60
        else
          date2(4)=date2(4)-(1+int(abs(date2(5))/60))
          date2(5)=date2(5)+(1+int(abs(date2(5))/60))*60
        endif
      endif

      if (date2(4).ge.24) then
        date2(3)=date2(3)+int(date2(4)/24)
        date2(4)=date2(4)-int(date2(4)/24)*24
      endif
      if (date2(4).lt.0) then
        if (mod(date2(4),24).eq.0) then
          date2(3)=date2(3)-int(abs(date2(4))/24)
          date2(4)=date2(4)+int(abs(date2(4))/24)*24
        else
          date2(3)=date2(3)-(1+int(abs(date2(4))/24))
          date2(4)=date2(4)+(1+int(abs(date2(4))/24))*24
        endif
      endif

  100 if (date2(3).gt.idays(date2(2))) then
        if ((date2(2).eq.2).and.(mod(date2(1),4).eq.0)) idays(2)=29
        date2(3)=date2(3)-idays(date2(2))
        if (idays(2).eq.29) idays(2)=28
        date2(2)=date2(2)+1
        if (date2(2).gt.12) then
*         date2(1)=date2(1)+int(date2(2)/12)
*         date2(2)=date2(2)-int(date2(2)/12)*12
          date2(1)=date2(1)+1
          date2(2)=date2(2)-12
        endif
        if (date2(2).lt.1) then
          date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
          date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
        endif
        goto 100
      endif     
  200 if (date2(3).lt.1) then
        date2(2)=date2(2)-1
        if (date2(2).gt.12) then
          date2(1)=date2(1)+int(date2(2)/12)
          date2(2)=date2(2)-int(date2(2)/12)*12
        endif
        if (date2(2).lt.1) then
          date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
          date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
        endif
        if ((date2(2).eq.2).and.(mod(date2(1),4).eq.0)) idays(2)=29
        date2(3)=date2(3)+idays(date2(2))
        if (idays(2).eq.29) idays(2)=28
        goto 200
      endif

      if (date2(2).gt.12) then
        date2(1)=date2(1)+int(date2(2)/12)
        date2(2)=date2(2)-int(date2(2)/12)*12
      endif
      if (date2(2).lt.1) then
        date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
        date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
      endif

      if (date2(1).lt.1000) then
      if (date2(1).ge.100) date2(1)=date2(1)-100
      endif

      return
      end
