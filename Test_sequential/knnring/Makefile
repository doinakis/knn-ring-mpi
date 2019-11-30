#
#Implementation of Vantage point tree pthreads
#*Doinakis Michail
#*e-mail: doinakis@ece.auth.gr
#
SHELL	:=	/bin/bash
CC	=	gcc-7
CFLAGS	=	-O3	-Wall -g -Werror
INCLUDES	=
LDFLAGS	=
LIBS	= -lm	-lopenblas
SRC	=	./src/knnring
TYPES	=	sequential

MAIN	=	./src/main
.PRECIOUS:	%.a

all:	$(addprefix	$(MAIN)_,	$(TYPES))
		mv $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES))) ./lib/


lib:	$(addsuffix	.a,	$(addprefix	$(SRC)_,$(TYPES)))
		mv $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES))) ./lib/

$(MAIN)_%:	$(MAIN).c	$(SRC)_%.a
	$(CC)	$(CFLAGS)	$(INCLUDES)	-o	$@	$^	$(LDFLAGS)	$(LIBS)

.o.a:
	ar	rcs	$@	$<
.c.o:
	$(CC)	$(CFLAGS)	$(INCLUDES)	-o	$@	-c	$<

clean:
	$(RM)	*.o	*~	$(addprefix	$(MAIN)_,	$(TYPES))	knnring_*.a
