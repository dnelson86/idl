
# in ~/.bashrc include the line:
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/n/sw/idl-8.1/idl/idl81/external/include
# comment out INT64_PRECISION OPT for 32bit version

OPT = -DUINT64_PRECISION
#OPT = -DINT64_PRECISION

EXEC = CalcMatch_uint64.so

SRCS   = main.c fhtr.c

OBJS = $(SRCS:.c=.o)
INCL    = -L/usr/lib64 -L/usr/local/lib64 -shared

CFLAGS =  -fPIC -m64 -O3 -Wall

LIBS   =  -lm

CC     =  cc 

$(EXEC): $(OBJS)
	$(CC) $(INCL) $(CFLAGS) $(OBJS) $(LIBS) $(OPT) -o $(EXEC)

$(OBJS):
	$(CC) -c $(SRCS) $(CFLAGS) $(OPT)

clean:
	rm -f $(OBJS) *.so
