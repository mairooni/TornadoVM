#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void testInfinity(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *fin, __global uchar *fout)
{
  int i_5, i_4, i_10, i_13, i_12, i_3, i_2; 
  ulong ul_8, ul_1, ul_0, ul_11; 
  float f_9; 
  long l_6, l_7; 

  // BLOCK 0
  ul_0  =  (ulong) fin;
  ul_1  =  (ulong) fout;
  i_2  =  get_global_size(0);
  i_3  =  get_global_id(0);
  // BLOCK 1 MERGES [0 5 ]
  i_4  =  i_3;
  for(;i_4 < 128;)
  {
    // BLOCK 2
    i_5  =  i_4 + 6;
    l_6  =  (long) i_5;
    l_7  =  l_6 << 2;
    ul_8  =  ul_0 + l_7;
    f_9  =  *((__global float *) ul_8);
    i_10  =  i_2 + i_4;
    ul_11  =  ul_1 + l_7;
    i_12  =  isless(64.0F, f_9);
    if(i_12 == 1)
    {
      // BLOCK 3
      *((__global float *) ul_11)  =  InfinityF;
    }  // B3
    else
    {
      // BLOCK 4
      *((__global float *) ul_11)  =  -InfinityF;
    }  // B4
    // BLOCK 5 MERGES [3 4 ]
    i_13  =  i_10;
    i_4  =  i_13;
  }  // B5
  // BLOCK 6
  return;
}  //  kernel
