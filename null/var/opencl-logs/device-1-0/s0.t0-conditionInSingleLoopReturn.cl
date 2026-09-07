#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void conditionInSingleLoopReturn(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *a, __global uchar *b, __global uchar *c)
{
  ulong ul_5, ul_15, ul_0, ul_1, ul_17, ul_2, ul_13, ul_7; 
  int i_10, i_14, i_16, i_18, i_3, i_19, i_4, i_20, i_6, i_8; 
  long l_12, l_11; 
  bool b_9; 

  // BLOCK 0
  ul_0  =  (ulong) a;
  ul_1  =  (ulong) b;
  ul_2  =  (ulong) c;
  i_3  =  get_global_id(0);
  // BLOCK 1 MERGES [0 3 ]
  i_4  =  i_3;
  for(;i_4 < 2048;)
  {
    // BLOCK 2
    ul_5  =  ul_0 + 24L;
    i_6  =  *((__global int *) ul_5);
    ul_7  =  ul_1 + 24L;
    i_8  =  *((__global int *) ul_7);
    b_9  =  i_8 < i_6;
    if(b_9)
    {
      // BLOCK 4
      return;
    }  // B4
    else
    {
      // BLOCK 3
      i_10  =  i_4 + 6;
      l_11  =  (long) i_10;
      l_12  =  l_11 << 2;
      ul_13  =  ul_0 + l_12;
      i_14  =  *((__global int *) ul_13);
      ul_15  =  ul_1 + l_12;
      i_16  =  *((__global int *) ul_15);
      ul_17  =  ul_2 + l_12;
      i_18  =  i_14 + i_16;
      *((__global int *) ul_17)  =  i_18;
      i_19  =  get_global_size(0);
      i_20  =  i_19 + i_4;
      i_4  =  i_20;
    }  // B3
  }  // B3
  // BLOCK 5
  return;
}  // B5
}  //  kernel
