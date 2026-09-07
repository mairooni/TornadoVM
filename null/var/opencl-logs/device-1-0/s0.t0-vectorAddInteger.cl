#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void vectorAddInteger(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *a, __global uchar *b, __global uchar *c)
{
  int i_8, i_9, i_6, i_7, i_4, i_18, i_17, i_15, i_13; 
  long l_11, l_10; 
  ulong ul_12, ul_14, ul_1, ul_0, ul_16, ul_3, ul_2; 
  bool b_5; 

  // BLOCK 0
  ul_0  =  (ulong) a;
  ul_1  =  (ulong) b;
  ul_2  =  (ulong) c;
  ul_3  =  ul_0 + 24L;
  i_4  =  *((__global int *) ul_3);
  b_5  =  i_4 < 13;
  if(b_5)
  {
    // BLOCK 1
    i_6  =  get_global_size(0);
    i_7  =  get_global_id(0);
    // BLOCK 2 MERGES [1 3 ]
    i_8  =  i_7;
    for(;i_8 < 4096;)
    {
      // BLOCK 3
      i_9  =  i_8 + 6;
      l_10  =  (long) i_9;
      l_11  =  l_10 << 2;
      ul_12  =  ul_0 + l_11;
      i_13  =  *((__global int *) ul_12);
      ul_14  =  ul_1 + l_11;
      i_15  =  *((__global int *) ul_14);
      ul_16  =  ul_2 + l_11;
      i_17  =  i_13 + i_15;
      *((__global int *) ul_16)  =  i_17;
      i_18  =  i_6 + i_8;
      i_8  =  i_18;
    }  // B3
    // BLOCK 4
    return;
    else
    {
      // BLOCK 5
      return;
    }  // B5
  }  //  kernel
