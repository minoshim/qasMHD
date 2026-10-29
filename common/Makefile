DIR0 = ..
include $(DIR0)/Makefile.inc
TARGET = $(DIR0)/lib$(LIBNAME).a
OBJS = common_func.o common_interp.o mhd_eigen.o mhd_flux.o mhd_func.o
DEPS = $(OBJS:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(DEPS)

-include $(DEPS)
