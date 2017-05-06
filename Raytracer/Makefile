FILE=raytracer
FILEEXT=raytracer_ext
########
#   Directories
S_DIR=src
B_DIR=build
EXT_DIR=bin

########
#   Output
EXEC=$(B_DIR)/$(FILE)

# default build settings
CC_OPTS=-c -pipe -Wall -Wno-switch -ggdb -g3 -O3
LN_OPTS=-lX11
CC=c++ -fopenmp

########
#       SDL options
SDL_CFLAGS := $(shell sdl-config --cflags)
GLM_CFLAGS := -I$(GLMDIR)
SDL_LDFLAGS := $(shell sdl-config --libs)
#######

######
#	Openmp compilation
CPP = /usr/local/opt/llvm/bin/clang++ -fopenmp
CPPFLAGS = -I/usr/local/opt/llvm/include
LDFLAGS = -L/usr/local/opt/llvm/lib

########
#   This is the default action
all:Build


########
#   Object list
#
OBJ = $(B_DIR)/$(FILE).o


########
#   Objects
$(B_DIR)/$(FILE).o : $(S_DIR)/$(FILE).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

########
#   Main build rule
Build : $(OBJ) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS) $(LDFLAGS) $(CPPFLAGS)

########
# Extension stuff

# Extension object list
OBJEXT = $(B_DIR)/$(FILEEXT).o

# Extension objects
$(B_DIR)/$(FILEEXT).o : $(S_DIR)/$(FILEEXT).cpp $(S_DIR)/SDLauxiliary.h $(EXT_DIR)/TestModel2.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILEEXT).o $(S_DIR)/$(FILEEXT).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

EXECEXT=$(B_DIR)/$(FILEEXT)

extension: $(OBJEXT) Makefile
	$(CC) $(LN_OPTS) -o $(EXECEXT) $(OBJEXT) $(SDL_LDFLAGS) $(LDFLAGS) $(CPPFLAGS)

########


clean:
	rm -f $(B_DIR)/*
