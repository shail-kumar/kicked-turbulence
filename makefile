#~ all: shell
INCLUDE = -I/home/shailen/local/include/
OBJ = shell.o input_functions.o evolution.o processing_functions.o output_functions.o
COMPILER = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

shell: $(OBJ)
	$(COMPILER) $(DEBUG) $(OBJ) -o shell $(INCLUDE)
	
shell.o: shell.h shell.cc
	$(COMPILER) $(CFLAGS) shell.cc $(INCLUDE)
	
input_functions.o: shell.h input_functions.cc
	$(COMPILER) $(CFLAGS) input_functions.cc $(INCLUDE)
	
evolution.o: evolution.h evolution.cc shell.h
	$(COMPILER) $(CFLAGS) evolution.cc $(INCLUDE)
	
processing_functions.o: shell.h evolution.h processing_functions.cc
	$(COMPILER) $(CFLAGS) processing_functions.cc $(INCLUDE)
	
output_functions.o: shell.h output_functions.cc
	$(COMPILER) $(CFLAGS) output_functions.cc $(INCLUDE)
clean:
	rm *o shell
	
tar:
	tar cvf shell.tar shell.h shell.cc input_functions.cc \
	evolution.h evolution.cc processing_functions.cc \
	output_functions.cc makefile parameters.d initial_field.d
	
run:
	@./shell
