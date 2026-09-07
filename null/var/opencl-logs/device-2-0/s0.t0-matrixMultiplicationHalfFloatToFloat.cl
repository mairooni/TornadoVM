#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matrixMultiplicationHalfFloatToFloat(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *A, __global uchar *B, __global uchar *C, __private int size)
{
  long l_7, l_6, l_13, l_14; 
  half half_9, half_10; 
  float f_16; 
  short sh_11; 
  ulong ul_1, ul_8, ul_15, ul_0; 
  int i_4, i_3, i_2, i_17, i_5, i_12; 

  // BLOCK 0
  ul_0  =  (ulong) A;
  ul_1  =  (ulong) C;
  i_2  =  get_global_size(0);
  i_3  =  get_global_id(0);
  // BLOCK 1 MERGES [0 2 ]
  i_4  =  i_3;
  for(;i_4 < 256;)
  {
    // BLOCK 2
    i_5  =  i_4 + 12;
    l_6  =  (long) i_5;
    l_7  =  l_6 << 1;
    ul_8  =  ul_0 + l_7;
    half_9  =  *((__global half *) ul_8);
    half_10  =  half_9 + 16L;
    sh_11  =  *((__global short *) half_10);
    i_12  =  i_4 + 6;
    l_13  =  (long) i_12;
    l_14  =  l_13 << 2;
    ul_15  =  ul_1 + l_14;
    f_16  =  convert_float(sh_11);
    *((__global float *) ul_15)  =  f_16;
    i_17  =  i_2 + i_4;
    i_4  =  i_17;
  }  // B2
  // BLOCK 3
  return;
}  //  kernel
