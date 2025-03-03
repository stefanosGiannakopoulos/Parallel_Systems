.KEEP_STATE:

CPP          = g++
NVCC        = /usr/local/cuda-12.1/bin/nvcc

CFLAGS  = -O3
NVCCFLAGS   =$(CFLAGS) --ptxas-options=-v -arch=sm_60\
#-gencode=arch=compute_35,code=sm_35\
#-gencode=arch=compute_52,code=sm_52\
-gencode=arch=compute_60,code=sm_60\
#-gencode=arch=compute_61,code=sm_61\
#-gencode=arch=compute_75,code=sm_75\
#-gencode=arch=compute_80,code=sm_80\
#-gencode=arch=compute_86,code=sm_86\
#-gencode=arch=compute_87,code=sm_87\
#-gencode=arch=compute_86,code=compute_86


all: cuda_properties

cuda_properties: cuda_properties.o
	$(NVCC) $(CFLAGS) $^ -o $@

cuda_properties.o: cuda_properties.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f cuda_properties.o cuda_properties

