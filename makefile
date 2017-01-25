CC = gcc
TYPE = %.c
CFLAGS = -Wall
LDFLAGS = -lm 
PROGNAME = testQR
OPTIONS = $(OPT)
CFILES = testQR.c easyQR.c
OFILES = $(CFILES:%.c=%.o)

ifeq ($(strip $(OPT)),debug)
	
	CFLAGS = -fopenmp -Wall
	OPTIONS = -g -pg	

endif

ifeq ($(strip $(OPT)),)
		
	CFLAGS = -fopenmp
	OPTIONS = -march=native

endif

ifeq ($(strip $(OPT)),O3)
	
	CFLAGS = -fopenmp
	OPTIONS = -O3

endif

ifeq ($(strip $(OPT)),serial)

	CFLAGS = -Wall
	OPTIONS = -g -pg

endif




compile: $(OFILES)

	$(CC) $(CFLAGS) $(OPTIONS) -o $(PROGNAME) $(OFILES) $(LDFLAGS)

%.o : $(TYPE)
	$(CC) $(CFLAGS) $(OPTIONS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o $(PROGNAME)
