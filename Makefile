CC = c++
CFLAGS = -std=c++17 -O2
LIBFLAGS = -larmadillo -lconfig++
DBGFLAGS = -g -O0 -fsanitize=address -fsanitize=bounds -lubsan

BIN = ./bin
SRC = ./src
OBJ = ./obj
SET = ./stg
DAT = ./dat

TARGET = $(BIN)/SIMPLE

INCS = $(wildcard $(SRC)/*.H)
SRCS = $(wildcard $(SRC)/*.C)
OBJS = $(patsubst $(SRC)/%.C, $(OBJ)/%.o, $(SRCS))
INCDIRS = -I./ $(addprefix -I, $(SRC))

PHONY := $(TARGET)
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBFLAGS)

$(OBJ)/%.o: $(SRC)/%.C $(INCS)
	$(CC) $(CFLAGS) -o $@ -c $< $(INC_DIRS)

PHONY += test
FINDDATS = $(wildcard $(DAT)/*.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
test: 
	$(TARGET) 0 1
	$(TARGET) 0 5
	$(TARGET) 0 10
	$(TARGET) 0 20
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plot2D.plt ; \
	done

.PHONY: $(PHONY)
