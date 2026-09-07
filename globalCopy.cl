#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void floatGlobalCopy(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *a, __global uchar *b, __global uchar *c)
{
  long l_9, l_10; 
  ulong ul_2, ul_1, ul_0, ul_15, ul_5, ul_13, ul_11; 
  float f_12, f_14, f_17, f_16; 
  int i_4, i_3, i_18, i_8, i_6; 

  // BLOCK 0
  ul_0  =  (ulong) a;
  ul_1  =  (ulong) b;
  ul_2  =  (ulong) c;
  i_3  =  get_global_size(0);
  i_4  =  get_global_id(0);
  // BLOCK 1 MERGES [0 2 ]
  ul_5  =  ul_0;
  i_6  =  i_4;
  for(;i_6 < 16;)
  {
    // BLOCK 2
    ulong ul_7 = ul_5;
    i_8  =  i_6 + 6;
    l_9  =  (long) i_8;
    l_10  =  l_9 << 2;
    ul_11  =  ul_7 + l_10;
    f_12  =  *((__global float *) ul_11); // value1 -> base is ul_5 
    ul_13  =  ul_1 + l_10;
    f_14  =  *((__global float *) ul_13); 
    ul_15  =  ul_2 + l_10;
    f_16  =  f_14 + 3.0F; // value2 
    f_17  =  f_16 + f_12;
    *((__global float *) ul_15)  =  f_17; // value2 -> base is ul_1
    i_18  =  i_3 + i_6;
    ul_5  =  ul_1; // <- race condition
    i_6  =  i_18;
  }  // B2
  // BLOCK 3
  return;
}  //  kernel
