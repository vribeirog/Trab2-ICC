CC = gcc

# likwid home
    LIKWID_HOME = /home/soft/likwid

CFLAGS = -g -O3 -march=native -mavx -fopt-info-vec -I${LIKWID_HOME}/include -DLIKWID_PERFMON
LFLAGS = -lm -L${LIKWID_HOME}/lib -llikwid

PROG = cgSolver
MODULES = utils cgSolver matvet sislin gradconj $(PROG)
OBJS = $(addsuffix .o,$(MODULES))
SRCS = $(addsuffix .c,$(MODULES)) $(addsuffix .h,$(MODULES))

# Lista de arquivos para distribuição
DISTFILES = *.c *.h Makefile README.md
DISTDIR = ibs20-vrg20

.PHONY: clean purge dist all

all: $(PROG)

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

debug:   CFLAGS+=-D__DEBUG__
debug: $(PROG)

clean:
	@echo "Limpando sujeira ....."
	@rm -rf core *~ *.bak
	@rm -f a.out *.o $(PROG)

purge: clean
	@echo "Fazendo a faxina ....."
	@rm -f a.out *.o $(PROG)


dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
