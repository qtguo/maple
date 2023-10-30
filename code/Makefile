CFLAGS += -Wall -std=c++14 -O3
LDFLAGS +=

all: baim format_graph gene_real

baim: *.cpp *.h
		g++ ${CFLAGS}  baim.cpp  sfmt/SFMT.c  -o $@ ${LDFLAGS}  -lpthread

format_graph: utils/format_graph.cpp
	g++ ${CFLAGS} $^ -o $@

gene_real: utils/generate_realization.cpp sfmt/SFMT.c
	g++ ${CFLAGS} $^ -o $@

clean: 
	rm baim format_graph gene_real
