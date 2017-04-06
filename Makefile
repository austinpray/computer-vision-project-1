.PHONY: all

CPPFLAGS=-I/usr/local/include/opencv -I/usr/local/include/opencv2 -I./eigen -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_ml -lopencv_video -lopencv_features2d -lopencv_calib3d -lopencv_objdetect -lopencv_stitching -lopencv_imgcodecs

p1: p1.cpp
	llvm-g++ -g -o p1 p1.cpp $(CPPFLAGS)

p2: p2.cpp
	llvm-g++ -g -o p2 p2.cpp $(CPPFLAGS)

p3: p3.cpp
	llvm-g++ -g -o p3 p3.cpp $(CPPFLAGS)

p4: p4.cpp
	llvm-g++ -g -o p4 p4.cpp $(CPPFLAGS)

all: p1 p2 p3 p4;
