SAM_DIR=

CC=gcc

# Our headers are in include/, and we also need the SAM headers, which
# are in the root directory of the samtools distribution.
CFLAGS=-Iinclude -I$(SAM_DIR) -Wall

# The location of the compiled sam library
SAM_LIB=$(SAM_DIR)/libbam.a

# SAM requires libz
LDFLAGS=-lz

bins=bin/findexons bin/dumptranscripts

all : check_sam $(bins) 

.PHONY : check_sam

check_sam :
ifndef SAM_DIR
	$(error Pre-pade requires the samtools library. Please add SAM_DIR=DIR to your command, where DIR is the directory containing the samtool headers and compiled libraries.)
endif

%.o : src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

bin/dumptranscripts : geneindex.o dumptranscripts.o $(SAM_LIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

bin/findexons : findexons.o geneindex.o $(SAM_LIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

bin/testgeneindex : testgeneindex.o geneindex.o testutils.o $(SAM_LIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

bin/testsamutils : testsamutils.o geneindex.o testutils.o $(SAM_LIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

test_geneindex : bin/testgeneindex
	bin/testgeneindex

test : bin/testgeneindex bin/testsamutils
	bin/testgeneindex
	bin/testsamutils
	@echo
	@echo "It looks like all tests passed."

clean :
	rm -f *.o $(bins)
	rm -f `find . -name \*~`

