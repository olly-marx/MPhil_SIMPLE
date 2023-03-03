CC = c++
CFLAGS = -std=c++17 -Wpedantic -Wall -Wextra -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -O2
LIBFLAGS = -larmadillo -lconfig++
DBGFLAGS = -g -O0 -fsanitize=address -fsanitize=bounds -lubsan

BIN = ./bin
SRC = ./src
OBJ = ./obj
SET = ./stg
DAT = ./dat

TARGET = $(BIN)/fvMesh

INCS = $(wildcard $(SRC)/*.H)
SRCS = $(wildcard $(SRC)/*.C)
OBJS = $(patsubst $(SRC)/%.C, $(OBJ)/%.o, $(SRCS))
INCDIRS = -I./ $(addprefix -I, $(SRC))

PHONY := $(TARGET)
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBFLAGS)

$(OBJ)/%.o: $(SRC)/%.C $(INCS)
	$(CC) $(CFLAGS) -o $@ -c $< $(INC_DIRS)

PHONY += debug
debug: 
	CFLAGS += $(DBGFLAGS)
	$(TARGET)

PHONY += release
release: 
	CFLAGS += $(OPTFLAG)
	$(TARGET)

PHONY += refactor
refactor: 
	CFLAGS += $(REFFLAGS)
	$(TARGET)

PHONY += plotGifs
FINDDATS = $(wildcard $(DAT)/*.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plotGifs:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plotGif.plt ; \
	done

PHONY += plot2D
FINDDATS = $(wildcard $(DAT)/*.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plot2D:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plot2D.plt ; \
	done

PHONY += clean
clean:
	rm $(OBJ)/*.o
	rm $(BIN)/*

PHONY += cleanresults
cleanresults:
	rm $(DAT)/*.dat
	rm *.gif

.PHONY: $(PHONY)
