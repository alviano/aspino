BUILD = release


########## Available targets

cxxflags.debug =    -DTRACE_ON -DSAFE_EXIT -g
linkflags.debug =   -g

cxxflags.trace =    -DNDEBUG -DTRACE_ON -O3
linkflags.trace =

cxxflags.release =  -DNDEBUG -O3
linkflags.release =

cxxflags.gprof =    -DNDEBUG -DSAFE_EXIT -O3 -g -pg
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
LINKFLAGS = $(linkflags.$(BUILD))
LIBS = -lm -lz

SRCS = $(shell find $(SOURCE_DIR) -name '*.cc')

OBJS = $(patsubst $(SOURCE_DIR)%.cc,$(BUILD_DIR)%.o, $(SRCS))
DEPS = $(patsubst $(SOURCE_DIR)%.cc,$(BUILD_DIR)%.d, $(SRCS))

EXTERN = glucose-syrup

.PHONY: $(EXTERN)

all: $(EXTERN) $(BINARY)

glucose-syrup:
	@if [ ! -d src/glucose-syrup ]; then \
	    echo "************************************************************"; \
	    echo "* Hey! Directory src/glucose-syrup is missing.             *"; \
	    echo "* Did you run bootstrap.sh?                                *"; \
	    echo "************************************************************"; \
	    exit 1; \
    fi

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.d: $(SOURCE_DIR)/%.cc
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MM -MT '$(@:.d=.o)' $< -MF $@
	
$(BINARY): $(OBJS) $(DEPS)
	$(LINK) $(OBJS) -o $(BINARY) $(LINKFLAGS) $(LIBS)

static: $(OBJS) $(DEPS)
	$(LINK) $(OBJS) -static -o $(BINARY) $(LINKFLAGS) $(LIBS)

run: $(BINARY)
	./$(BINARY)

########## Tests
include tests/Makefile.tests.inc


########## Clean
.PHONY: clean-dep clean distclean

clean-dep:
	rm -f $(DEPS)
clean: clean-dep
	rm -f $(OBJS)

distclean: clean
	rm -fr $(BUILD_DIR)

-include $(DEPS)
