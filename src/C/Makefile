CC		= mpicc
CPPFLAGS	= -MMD -MP
CFLAGS		= -pedantic -Wall -std=gnu99
INCLUDES	= -I/usr/include/atlas
LIBS		= -latlas -llapack -lblas -lm

SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c,%.o,$(SRCS))
DEPS = $(patsubst %.c,%.d,$(SRCS))

TGT = 1D-pure

all: $(TGT)

$(TGT): $(OBJS)
	$(LINK.o) -o $@ $(OBJS) $(LIBS) 

%.o: %.c
	$(COMPILE.c) -o $@ $(INCLUDES) $<

-include $(DEPS)

clean:
	@$(RM) $(TGT) $(OBJS) $(DEPS)
