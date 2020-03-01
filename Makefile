CC=g++
EXE=rt
INC_DIR=include
SRC_DIR=src
OBJ_DIR=obj
SRC=$(wildcard $(SRC_DIR)/*.cpp)
OBJ=$(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
CPPFLAGS=-Iinclude
CFLAGS=-Wall
LDFLAGS=
LDLIBS=

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir $@

clean:
	$(RM) $(OBJ)
