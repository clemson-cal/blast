VAPOR_HOME = ../vapor

EXE = blast-ui
IMGUI_DIR = imgui
IMPLOT_DIR = implot
IMGUI_FILEBROWSER = imgui-filebrowser

SOURCES = main.mm
SOURCES += \
	$(IMGUI_DIR)/imgui.cpp \
	$(IMGUI_DIR)/imgui_demo.cpp \
	$(IMGUI_DIR)/imgui_draw.cpp \
	$(IMGUI_DIR)/imgui_tables.cpp \
	$(IMGUI_DIR)/imgui_widgets.cpp \

SOURCES += \
	$(IMGUI_DIR)/backends/imgui_impl_sdl2.cpp \
	$(IMGUI_DIR)/backends/imgui_impl_metal.mm \

SOURCES += \
	$(IMPLOT_DIR)/implot.cpp \
	$(IMPLOT_DIR)/implot_items.cpp \
	$(IMPLOT_DIR)/implot_demo.cpp \

OBJS = $(addsuffix .o, $(basename $(notdir $(SOURCES))))

LIBS += -framework Metal -framework MetalKit -framework Cocoa -framework IOKit -framework CoreVideo -framework QuartzCore
LIBS += `sdl2-config --libs`
LIBS += -L/usr/local/lib
LIBS += -lhdf5

CXXFLAGS = -std=c++17 -I/usr/local/include
CXXFLAGS += `sdl2-config --cflags`
CXXFLAGS += -Wall -Wformat
CXXFLAGS += -O0
CXXFLAGS += -I$(IMGUI_DIR)
CXXFLAGS += -I$(IMGUI_DIR)/backends
CXXFLAGS += -I$(IMPLOT_DIR)
CXXFLAGS += -I$(IMGUI_FILEBROWSER)
CXXFLAGS += -I$(VAPOR_HOME)/include

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:$(IMGUI_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:$(IMGUI_DIR)/backends/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:$(IMPLOT_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:%.mm
	$(CXX) $(CXXFLAGS) -ObjC++ -fobjc-weak -fobjc-arc -c -o $@ $<

%.o:$(IMGUI_DIR)/backends/%.mm
	$(CXX) $(CXXFLAGS) -ObjC++ -fobjc-weak -fobjc-arc -c -o $@ $<

all: $(EXE)
	@echo Build complete

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(EXE) $(OBJS)
