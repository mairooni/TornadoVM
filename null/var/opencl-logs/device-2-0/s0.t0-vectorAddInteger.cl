#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void vectorAddInteger(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *a, __global uchar *b, __global uchar *c)
{
  int i_8, i_6, i_7, i_4, i_5, i_16, i_17, i_14, i_12; 
  long l_9, l_10; 
  ulong ul_1, ul_2, ul_3, ul_13, ul_15, ul_0, ul_11; 

  // BLOCK 0
  ul_0  =  (ulong) a;
  ul_1  =  (ulong) b;
  ul_2  =  (ulong) c;
  ul_3  =  ul_0 + 24L;
  i_4  =  *((__global int *) ul_3);
  i_5  =  get_local_size(0);
  i_6  =  get_local_id(0);
  // BLOCK 1 MERGES [0 2 ]
  i_7  =  i_6;
  for(;i_4 >= i_7;)
  {
    // BLOCK 3
    return;
  }  // B3
  // BLOCK 2
  i_8  =  i_7 + 6;
  l_9  =  (long) i_8;
  l_10  =  l_9 << 2;
  ul_11  =  ul_0 + l_10;
  i_12  =  *((__global int *) ul_11);
  ul_13  =  ul_1 + l_10;
  i_14  =  *((__global int *) ul_13);
  ul_15  =  ul_2 + l_10;
  i_16  =  i_12 + i_14;
  *((__global int *) ul_15)  =  i_16;
  i_17  =  i_5 + i_7;
  i_7  =  i_17;
}  // B2
}  //  kernel
