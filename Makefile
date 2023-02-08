CC = c++
CFLAGS = -std=c++17 -Wpedantic -Wall -Wextra -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast -Wsign-conversion \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2
DBGFLAGS = -g -O0 -fsanitize=address -fsanitize=bounds -lubsan
OPTFLAG = -O2

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

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

$(OBJ)/%.o: $(SRC)/%.C $(INCS)
	$(CC) $(CFLAGS) -o $@ -c $< $(INC_DIRS)

$(OBJ)/%.o: $(SRC)/%.C $(INCS)
	$(CC) $(CFLAGS) -o $@ -c $< $(INC_DIRS)

debug: CFLAGS += $(DBGFLAGS)
debug : $(TARGET)

release: CFLAGS += $(OPTFLAG)
release: $(TARGET)

refactor: CFLAGS += $(REFFLAGS)
refactor: $(TARGET)

run: $(TARGET)

.PHONY: clean
clean:
	rm $(OBJ)/*.o
	rm $(BIN)/*

