TARGET  = Acou2D
#CC      = gcc
#CFLAGS  = -Wall -Wextra -g -pedantic -std=c99 -fopenmp 
#OPT     = -O3
CC      = mpicc
CFLAGS  = -Wall -g -openmp 
OPT     = -O3
LIBS    = -lm

DEFINES = -DVERB -DTEST_VEL

INCLUDE = -I./include

BIN_PTH = ./bin
SRC_PTH = ./source
SOURCE  = $(SRC_PTH)/Test_2DAcou.c \
          $(SRC_PTH)/init.c $(SRC_PTH)/output.c\
          $(SRC_PTH)/FM_acoustic.c $(SRC_PTH)/params.c \
	  $(SRC_PTH)/FFT.c $(SRC_PTH)/operation_complex.c\
	  $(SRC_PTH)/kernel_fd_acou.c 
	 
OBJ     = $(SOURCE:.c=.o) 

## Default rule executed
all: $(TARGET)
	@true

$(TARGET): $(OBJ)
	@echo
	@echo "=> Linking the target $@"
	@$(CC) $(CFLAGS) $(OPT) -o $(BIN_PTH)/$@ $^ $(LIBS)
	@echo '=> Link done'
	@-rm -f $(OBJ)

%.o: %.c
	@mkdir -p $(dir $@)
	@echo
	@echo "=> Compiling $<"
	@$(CC) $(CFLAGS) $(OPT) $(DEFINES) -c $< -o $@ $(INCLUDE)

# the rule to clean binaries
clean:
	@-rm -f bin/$(TARGET) $(OBJ)
	@echo "=> Clean done"
