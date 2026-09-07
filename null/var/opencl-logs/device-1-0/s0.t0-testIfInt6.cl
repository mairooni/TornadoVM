#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void testIfInt6(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *a)
{
  int i_4, i_3, i_2, i_1, i_11, i_13, i_8, i_9; 
  ulong ul_7, ul_0; 
  bool b_12, b_10; 
  long l_5, l_6; 

  // BLOCK 0
  ul_0  =  (ulong) a;
  i_1  =  get_global_size(0);
  i_2  =  get_global_id(0);
  // BLOCK 1 MERGES [0 8 ]
  i_3  =  i_2;
  for(;i_3 < 256;)
  {
    // BLOCK 2
    i_4  =  i_3 + 6;
    l_5  =  (long) i_4;
    l_6  =  l_5 << 2;
    ul_7  =  ul_0 + l_6;
    i_8  =  *((__global int *) ul_7);
    i_9  =  i_1 + i_3;
    b_10  =  i_8 < 0;
    if(b_10)
    {
      // BLOCK 3
    }  // B3
    else
    {
      // BLOCK 4
      i_11  =  *((__global int *) ul_7);
      b_12  =  i_11 < 2;
      if(b_12)
      {
        // BLOCK 5
        *((__global int *) ul_7)  =  100;
      }  // B5
      else
      {
        // BLOCK 6
      }  // B6
    }  // B4
    // BLOCK 7 MERGES [3 6 ]
    *((__global int *) ul_7)  =  200;
  }  // B7
  // BLOCK 8 MERGES [5 7 ]
  i_13  =  i_9;
  i_3  =  i_13;
}  // B8
// BLOCK 9
return;
}  //  kernel
