OPENCL_INCLUDE_PATH    = /home/_dev/OpenCL/SDK/include/
OPENCL_LIB_PATH        = /home/_dev/OpenCL/SDK/lib/x86_64/
OPENCL_LIB             = OpenCL

CC                     = g++
WARN                   = -Wall

INCLUDE                += -I. -I/usr/include -I/usr/local/include -I$(OPENCL_INCLUDE_PATH)
LIBPATH                += -L. -L/usr/lib -L$(OPENCL_LIB_PATH)
LIBRARIES              += -l$(OPENCL_LIB)

CFLAGS                 = -c -g -O3 $(INCLUDE) $(WARN)

OBJECTS                = gravlens.o
TARGET                 = gravlens

RM                     = rm -f

$(TARGET): $(OBJECTS)
	@echo "Linking '$@'"
	$(CC) -o $@ $(OBJECTS) $(LIBPATH) $(LIBRARIES)

%.o: %.cpp
	@echo "Building '$@'"
	$(CC) $(CFLAGS) -c $<

gravlens.o: gravlens.cpp gravlens.hpp

clean:
	@echo "Cleaning $(OBJECTS) $(TARGET)"
	$(RM) $(OBJECTS) $(TARGET)

