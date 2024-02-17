VAPOR_HOME = ../vapor
GUI_HOME = .

EXE += blast_gui
SRC += \
	$(GUI_HOME)/imgui/backends/imgui_impl_metal.mm \
	$(GUI_HOME)/imgui/backends/imgui_impl_sdl2.cpp \
	$(GUI_HOME)/imgui/imgui.cpp \
	$(GUI_HOME)/imgui/imgui_demo.cpp \
	$(GUI_HOME)/imgui/imgui_draw.cpp \
	$(GUI_HOME)/imgui/imgui_tables.cpp \
	$(GUI_HOME)/imgui/imgui_widgets.cpp \
	$(GUI_HOME)/implot/implot.cpp \
	$(GUI_HOME)/implot/implot_demo.cpp \
	$(GUI_HOME)/implot/implot_items.cpp \
	blast_gui.mm
OBJ += $(addsuffix .o, $(basename $(SRC)))
LIB += -framework Metal -framework Cocoa
LIB += -lhdf5
LIB += $(shell sdl2-config --libs)

CXXFLAGS += -std=c++17
CXXFLAGS += -Wall
CXXFLAGS += -O3
CXXFLAGS += $(shell sdl2-config --cflags)
CXXFLAGS += -I$(GUI_HOME)/imgui
CXXFLAGS += -I$(GUI_HOME)/imgui/backends
CXXFLAGS += -I$(GUI_HOME)/implot
CXXFLAGS += -I$(GUI_HOME)/imgui-filebrowser
CXXFLAGS += -I$(VAPOR_HOME)/include

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
%.o : %.mm
	$(CXX) $(CXXFLAGS) -ObjC++ -fobjc-weak -fobjc-arc -c -o $@ $<
$(EXE) : $(OBJ)
	$(CXX) -o $@ $^ $(LIB)
show:
	echo $(LIB)
clean:
	rm -f $(EXE) $(OBJ)
