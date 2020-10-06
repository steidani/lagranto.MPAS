F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= newtime.o ${LAGRANTO}/lib/times.a
LIBS	= 

.f.o: $*.f $(INCS)
	${F77} -c ${FFLAGS} $*.f

newtime:	$(OBJS)
	${F77} -o newtime $(OBJS) $(LIBS)
