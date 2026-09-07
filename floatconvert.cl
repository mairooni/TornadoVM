#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matrixMultiplicationHalfFloatToFloat(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *A, __global uchar *B, __global uchar *C, __private int size)
{
  bool b_24, b_23, b_16, b_12; 
  int i_33, i_2, i_3, i_4, i_29, i_5; 
  long l_6, l_7, l_30, l_31; 
  half half_10, half_9; 
  float f_13, f_18, f_17, f_20, f_19, f_26, f_25, f_28, f_27; 
  ulong ul_1, ul_0, ul_32, ul_8; 
  short sh_11, sh_15, sh_14, sh_22, sh_21; 

  // BLOCK 0
  ul_0  =  (ulong) A;
  ul_1  =  (ulong) C;
  i_2  =  get_global_size(0);
  i_3  =  get_global_id(0);
  // BLOCK 1 MERGES [0 12 ]
  i_4  =  i_3;
  for(;i_4 < 4;)
  {
    // BLOCK 2
    i_5  =  i_4 + 12;
    l_6  =  (long) i_5;
    l_7  =  l_6 << 1;
    ul_8  =  ul_0 + l_7;
    half_9  =  *((__global half *) ul_8);
    f_20 = convert_float((float) half_9);
    // BLOCK 12 MERGES [6 9 11 10 ]
    i_29  =  i_4 + 6;
    l_30  =  (long) i_29;
    l_31  =  l_30 << 2;
    ul_32  =  ul_1 + l_31;
    *((__global float *) ul_32)  =  f_20;
    i_33  =  i_2 + i_4;
    i_4  =  i_33;
  }  // B12
  // BLOCK 13
  return;
}  //  kernel
