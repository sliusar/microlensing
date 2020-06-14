CC          = nvcc
LD          = nvcc
CFLAG       = 
PROG_NAME   = gravlens

SRC_DIR     = src
OBJECTS_DIR   = bin/objects
BIN_DIR     = bin

SRC_LIST = $(wildcard $(SRC_DIR)/*.cu)
OBJ_LIST = $(OBJECTS_DIR)/$(notdir $(SRC_LIST:.cu=.o))

.PHONY: all clean $(PROG_NAME) compile

all: $(PROG_NAME)

compile: 
	$(CC) -c $(CFLAG) $(SRC_LIST) -o $(OBJ_LIST)

$(PROG_NAME): compile
	$(LD) $(OBJ_LIST) -o $(BIN_DIR)/$@

clean:
	rm -f $(BIN_DIR)/$(PROG_NAME) $(OBJECTS_DIR)/*.o