#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void badCascadeKernel2(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics)
{
  int i_11, i_10, i_9, i_0, i_4, i_3, i_2, i_1, i_8, i_7, i_6; 
  bool b_5; 

  // BLOCK 0
  i_0  =  get_global_size(0);
  i_1  =  get_global_id(0);
  // BLOCK 1 MERGES [0 11 ]
  i_2  =  i_1;
  for(;i_2 < 100;)
  {
    // BLOCK 2
    // BLOCK 3 MERGES [2 10 ]
    i_3  =  1;
    i_4  =  0;
    for(;i_3 == 0;)
    {
      // BLOCK 6
    }  // B6
    // BLOCK 7 MERGES [6 5 ]
    // BLOCK 8 MERGES [7 9 ]
    i_6  =  i_3;
    i_7  =  0;
    for(;i_7 < i_2;)
    {
      // BLOCK 9
      i_8  =  i_7 + 1;
      i_9  =  (i_7 == 0) ? 1 : 0;
      i_6  =  i_9;
      i_7  =  i_8;
    }  // B9
    // BLOCK 10
    i_10  =  i_4 + 1;
    i_3  =  i_6;
    i_4  =  i_10;
    // BLOCK 4
    b_5  =  i_4 < 100;
    if(b_5)
    {
      // BLOCK 5
    }  // B5
    else
    {
      // BLOCK 11
      i_11  =  i_0 + i_2;
      i_2  =  i_11;
      break;
    }  // B11
  }  // B11
  // BLOCK 12
  return;
}  // B12
}  //  kernel
