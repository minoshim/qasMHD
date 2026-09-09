DIR0 = ../../..
include $(DIR0)/Makefile.inc
DIRI = $(DIR0)/common
DIRL = $(DIR0)
CPPFLAGS += -I$(DIRI)
LDFLAGS += -L$(DIRL)
LDLIBS += -lm -l$(LIBNAME)
LIBRARIES = $(DIRL)/lib$(LIBNAME).a
TARGET = a.out
OBJS = main.o mhd2d_class.o mhd2d_paras.o mhd2d_solve.o gmhd2d_class.o gmhd2d_init_.o gmhd2d_solve.o
DEPS = $(OBJS:.o=.d)
DATADIR = dat

.PHONY: all clean cdata

all: $(TARGET)

$(TARGET): $(OBJS) $(LIBRARIES)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(DIRL)/lib$(LIBNAME).a:
	$(MAKE) -C $(DIRL)/common

clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

cdata:
	$(RM) $(DATADIR)/*.dat

-include $(DEPS)
