#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void computeInit(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *this)
{
  ulong ul_5, ul_0, ul_2, ul_11, ul_7; 
  long l_10, l_9; 
  int i_12, i_8, i_1, i_3, i_4; 
  uint ui_6; 

  // BLOCK 0
  ul_0  =  (ulong) this;
  i_1  =  get_global_size(0);
  ul_2  =  ul_0 + 12L;
  i_3  =  get_global_id(0);
  // BLOCK 1 MERGES [0 2 ]
  i_4  =  i_3;
  for(;i_4 < 1024;)
  {
    // BLOCK 2
    ul_5  =  ul_0 + 12L;
    ui_6  =  *((__global uint *) ul_5);
    ul_7  =  ul_0 + ((ulong) ui_6 << 3);
    i_8  =  i_4 + 4;
    l_9  =  (long) i_8;
    l_10  =  l_9 << 2;
    ul_11  =  ul_7 + l_10;
    *((__global int *) ul_11)  =  100;
    i_12  =  i_1 + i_4;
    i_4  =  i_12;
  }  // B2
  // BLOCK 3
  return;
}  //  kernel
