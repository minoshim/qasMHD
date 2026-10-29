ifndef DIR0
$(error DIR0 must point to the repository root before including case.mk)
endif

ifndef CASE_SRCS
$(error CASE_SRCS must list the case-specific source files)
endif

include $(DIR0)/Makefile.inc

SERIAL1D_COMMON := ../common
ROOT_COMMON := $(DIR0)/common
DIRL := $(DIR0)
OBJDIR ?= build
TARGET ?= a.out
DATADIR ?= dat
HALL ?= 0

ifneq ($(filter $(HALL),0 1),$(HALL))
$(error HALL must be either 0 or 1)
endif

COMMON_NAMES := mhd1d_class mhd1d_solve
ifeq ($(HALL),1)
COMMON_NAMES += hmhd1d_class hmhd1d_solve
endif

CASE_OBJS := $(addprefix $(OBJDIR)/case_,$(CASE_SRCS:.cpp=.o))
COMMON_OBJS := $(addprefix $(OBJDIR)/common_,$(addsuffix .o,$(COMMON_NAMES)))
OBJS := $(CASE_OBJS) $(COMMON_OBJS)
DEPS := $(OBJS:.o=.d)

CPPFLAGS += -I. -I$(SERIAL1D_COMMON) -I$(ROOT_COMMON)
LDFLAGS += -L$(DIRL)
LDLIBS += -lm -l$(LIBNAME)
LIBRARIES := $(DIRL)/lib$(LIBNAME).a

.PHONY: all clean cdata

all: $(TARGET)

$(TARGET): $(OBJS) $(LIBRARIES)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

$(OBJDIR):
	mkdir -p $@

$(OBJDIR)/case_%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/common_%.o: $(SERIAL1D_COMMON)/%.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(DIRL)/lib$(LIBNAME).a:
	$(MAKE) -C $(DIRL)/common

clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

cdata:
	$(RM) $(DATADIR)/*.dat

-include $(DEPS)
