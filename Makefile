# $<: 所有依赖的目标集
# $@: 目标集
# $*  目标模式中“%”及其之前的部分
# %.o:%.c  模式规则
# @echo: 显示正在编译的模块

#---------------------------------
# Set correct variables and paths
#---------------------------------
CIMG_VERSION = 1.6.3_pre042715
CC           = g++
TARGET       = Paint
OBJ          =  Point.o  Color.o  Utility.o MaskedImage.o NNF.o FeaturePaint.o main.o
X11Path      = /usr/X11
ARMAPath     = /usr/local/Cellar
CFLAGS       =  -Wall -g -O2 -std=c++11 
LIBS         = -lm -lpthread -lblas  -llapack  -larmadillo  
CONF_LIBS    = -L$(X11Path)/lib -lpthread -lX11 -L$(ARMAPath)/lib
DIR         = -I$(X11Path)/include  -I$(ARMAPath)/include
PNG          =  $(CIMG_PNG_CFLAGS) $(CIMG_PNG_LIBS)
JPEG         = $(CIMG_JPEG_CFLAGS) $( -ljpeg)


#png
CIMG_PNG_CFLAGS  = -Dcimg_use_png
CIMG_PNG_LIBS    = -lpng -lz

#jpeg
CIMG_JPEG_CFLAGS = -Dcimg_use_jpeg
CIMG_JPEG_LIBS = -ljpeg

$(TARGET):$(OBJ)
	$(CC)  $(CONF_LIBS) $(LIBS) $(CIMG_PNG_LIBS) -o $@ $(OBJ)
%.o:%.cpp
	$(CC) $(CFLAGS)  $(DIR) $(CIMG_PNG_CFLAGS) -o $@ -c $<


.PHONY:clean

clean:
	rm -fr *.o $(TARGET)
