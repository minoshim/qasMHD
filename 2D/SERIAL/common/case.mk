ifndef DIR0
$(error DIR0 must point to the repository root before including case.mk)
endif

ifndef CASE_SRCS
$(error CASE_SRCS must list the case-specific source files)
endif

include $(DIR0)/Makefile.inc

SERIAL2D_COMMON := ../common
ROOT_COMMON := $(DIR0)/common
DIRL := $(DIR0)
OBJDIR ?= build
TARGET ?= a.out
DATADIR ?= dat
DISSIPATION ?= 0
GRAVITY ?= 0

ifneq ($(DISSIPATION),0)
ifneq ($(DISSIPATION),1)
$(error DISSIPATION must be either 0 or 1)
endif
endif
ifneq ($(GRAVITY),0)
ifneq ($(GRAVITY),1)
$(error GRAVITY must be either 0 or 1)
endif
endif

COMMON_NAMES := mhd2d_class mhd2d_solve
ifeq ($(DISSIPATION),1)
COMMON_NAMES += dmhd2d_class dmhd2d_solve
endif
ifeq ($(GRAVITY),1)
COMMON_NAMES += gmhd2d_class gmhd2d_solve
endif

CASE_OBJS := $(addprefix $(OBJDIR)/case_,$(CASE_SRCS:.cpp=.o))
COMMON_OBJS := $(addprefix $(OBJDIR)/common_,$(addsuffix .o,$(COMMON_NAMES)))
OBJS := $(CASE_OBJS) $(COMMON_OBJS)
DEPS := $(OBJS:.o=.d)

# Compile shared sources separately for each case, using that case's mymacros.hpp.
CPPFLAGS += -I. -I$(SERIAL2D_COMMON) -I$(ROOT_COMMON)
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

$(OBJDIR)/common_%.o: $(SERIAL2D_COMMON)/%.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(DIRL)/lib$(LIBNAME).a:
	$(MAKE) -C $(DIRL)/common

clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

cdata:
	$(RM) $(DATADIR)/*.dat

-include $(DEPS)
