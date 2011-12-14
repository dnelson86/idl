
# in ~/.bashrc include the line:
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/n/sw/idl-8.1/idl/idl81/external/include

EXEC   = CalcHSML.so

SRCS   = main.c tree.c

OBJS   = $(SRCS:.c=.o)
INCL   = 

CFLAGS =  -fPIC -m64 -O3 -g -Wall

LIBS   =  -lm

CC     =  cc 

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib64 -L/usr/local/lib64 -shared  $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm $(OBJS) $(EXEC)

