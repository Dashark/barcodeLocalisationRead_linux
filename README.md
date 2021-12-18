# barcodeLocalisationRead_linux
Barcode reader that crops  barcode from image and then reads the barcode

This is an implementation of the algorithm described in the paper "Robust Recognition of 1-D Barcodes Using Camera Phones Steffen Wachenfeld, Sebastian Terlunen, Xiaoyi Jiang Computer Vision and Pattern Recognition Group, Department of Computer Science, University of Munster, Germany"

This uses opencv libraries for the initial processing of the input jpg files. The opencv libraries are used to extract teh RGB datat from the image. If the RGB data is already known, this can be removed.

To compile, you will have to include opencv c++ libraries. Most Linux and UNIX OS's will have a package manager from which opencv libraries can be installed with a single command. For example on Open Suse (I used WSL,Opensuse Leap from Microsoft Store on Windows 10)

sudo zypper install opencv

sudo zypper install opencv-devel

to compile g++ -g main.cpp -lopencv_core -lopencv_imgproc -lopencv_imgcodecs -lpthread -o main

or you can use CMake, my CMakeLists.txt is

cmake_minimum_required(VERSION 3.10)

project(Test)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(${DIR}/headers /opt/local/include) link_directories (${DIR}/src/utils /opt/local/lib)

ADD_EXECUTABLE(main main.cpp) TARGET_LINK_LIBRARIES(main opencv_core opencv_imgproc -lopencv_imgcodecs ${EXTERNAL_LIBS})

to run

./main barcode_1.jpg true

the first parameter is the file path of the jpeg image and the second parameter turns debug output on or off

Ten sample images have been uploaded, of which all are correctly read.

If you are using Windows 10 or later, you can make use of the Windows Subsystem Linux to install various flavours of virtual Linux machines to compile and run the code.
