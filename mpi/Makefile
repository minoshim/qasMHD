DIR0 = ..
include $(DIR0)/Makefile.inc
TARGET = $(DIR0)/lib$(LIBMPI).a
OBJS = common_mpi.o
DEPS = $(OBJS:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(DEPS)

-include $(DEPS)
