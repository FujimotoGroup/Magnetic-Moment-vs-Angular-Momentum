#COMPILER  = icpc
COMPILER  = g++
#CFLAGS    = -std=c++11 -qmkl -g -MMD -MP -Wall -Wextra -Winit-self -Wno-missing-field-initializers
CFLAGS    = -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm
CFLAGS   += -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/local/include
CFLAGS   += -L/usr/local/lib -lgts -lglib-2.0 -lm
ifeq "$(shell getconf LONG_BIT)" "64"
	LDFLAGS =
else
	LDFLAGS =
endif
HEADDIR   = ./header
SRCDIR    = ./src
OBJDIR    = ./obj
EXEDIR    = ./bin
LIBDIR    = ./lib
HEADS     = parameters.hpp
LIBS      = initialize.cpp \
			hamiltonian.cpp \
			dispersion.cpp
INCLUDE   = -I$(LIBDIR) -I$(HEADDIR)
TARGET  = $(notdir $(SOURCES:%.cpp=%))
SOURCES   = $(wildcard $(SRCDIR)/*.cpp)
HEADERS = $(addprefix $(HEADDIR)/, $(notdir $(HEADS)))
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(LIBS:%.cpp=%.o)))
DEPENDS   = $(OBJECTS:.o=.d)
###############################################################################
.SUFFIXES: .cpp

.PHONY : clean all
###############################################################################

$(TARGET): $(OBJECTS)
	$(COMPILER) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(SRCDIR)/$@.cpp $(CFLAGS)

$(OBJDIR)/%.o: $(LIBDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

all: clean $(TARGET)

clean:
	-rm -f $(OBJECTS) $(DEPENDS) $(TARGET)

-include $(DEPENDS)

%.o: %.cpp
	@:
