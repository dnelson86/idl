
# in ~/.bashrc include the line:
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/n/sw/idl-8.1/idl/idl81/external/include
# comment out TWODIMS OPT for 3D version

#OPT = -DTWODIMS
EXEC = CalcSphMap_3D.so

SRCS   = main.c

OBJS = $(SRCS:.c=.o)
INCL    = -L/usr/lib64 -L/usr/local/lib64 -shared

CFLAGS =  -fPIC -m64 -O3 -g -Wall

LIBS   =  -lm

CC     =  cc 

$(EXEC): $(OBJS)
	$(CC) $(INCL) $(CFLAGS) $(OBJS) $(LIBS) $(OPT) -o $(EXEC)

$(OBJS):
	$(CC) -c $(SRCS) $(CFLAGS) $(OPT)

clean:
	rm -f $(OBJS)
