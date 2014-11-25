BUILD = release


########## Available targets

cxxflags.debug =    -DTRACE_ON
linkflags.debug =

cxxflags.trace =    -DNDEBUG -DTRACE_ON -O3
linkflags.trace =

cxxflags.release =  -DNDEBUG -O3
linkflags.release =

cxxflags.gprof =    -DNDEBUG -O3 -g -pg
linkflags.gprof =   -g -pg

cxxflags.stats =    -DNDEBUG -DSTATS_ON -O3
linkflags.stats =


########## Shared options

SOURCE_DIR = src
BUILD_DIR = build/$(BUILD)

BINARY = $(BUILD_DIR)/aspino
GCC = g++
CXX = $(GCC)
CXXFLAGS = $(cxxflags.$(BUILD)) -Isrc/glucose-syrup -Wall -Wextra
LINK = $(GCC)
LINKFLAGS = $(linkflags.$(BUILD)) -lm -lz -lgflags -lpthread

SRCS = $(shell find $(SOURCE_DIR) -name '*.cc')

OBJS = $(patsubst $(SOURCE_DIR)%.cc,$(BUILD_DIR)%.o, $(SRCS))
DEPS = $(patsubst $(SOURCE_DIR)%.cc,$(BUILD_DIR)%.d, $(SRCS))

all: $(BINARY)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.d: $(SOURCE_DIR)/%.cc
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MM -MT '$(@:.d=.o)' $< -MF $@
	
$(BINARY): $(OBJS) $(DEPS)
	$(LINK) $(LINKFLAGS) $(LIBS) $(OBJS) -o $(BINARY)

static: $(OBJS) $(DEPS)
	$(LINK) $(LINKFLAGS) $(LIBS) $(OBJS) -static -o $(BINARY)

run: $(BINARY)
	./$(BINARY)

########## Tests
include tests/Makefile.tests.inc


########## Clean

clean-dep:
	rm -f $(DEPS)
clean: clean-dep
	rm -f $(OBJS)

distclean: clean
	rm -fr $(BUILD_DIR)

-include $(DEPS)
