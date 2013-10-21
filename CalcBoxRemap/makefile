
# in ~/.bashrc include the line:
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/n/sw/idl-8.1/idl/idl81/external/include

OPT = -DSKIPZ_STRIDE2
EXEC = CalcBoxRemap_skipZ.so

SRCS   = main.cpp cuboid.cpp

OBJS = $(SRCS:.cpp=.o)
INCL    = -L/usr/lib64 -L/usr/local/lib64 -shared

CFLAGS =  -fPIC -m64 -O3 -g -Wall #-march=bdver1 -mprefer-avx128 #ITC/holynoke nodes only

LIBS   =  -lm

CC     =  g++

$(EXEC): $(OBJS)
	$(CC) $(INCL) $(CFLAGS) $(OBJS) $(LIBS) $(OPT) -o $(EXEC)

$(OBJS):
	$(CC) -c $(SRCS) $(CFLAGS) $(OPT)

clean:
	rm -f $(OBJS)
