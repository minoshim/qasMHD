DIR0 := ../../..
include $(DIR0)/Makefile.inc

TARGET ?= merge.out
OBJDIR ?= build
OBJECT := $(OBJDIR)/merge.o

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJECT)
	$(CXX) $(LDFLAGS) $^ -o $@

$(OBJECT): merge.cpp Makefile $(DIR0)/Makefile.inc | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $@

clean:
	$(RM) $(TARGET) $(OBJECT) $(OBJECT:.o=.d)

-include $(OBJECT:.o=.d)
